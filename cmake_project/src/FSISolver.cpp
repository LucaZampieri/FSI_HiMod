#include "../include/FSISolver.hpp"

#include <cmath>
#include <fstream>

namespace LifeV
{

  FSISolver::FSISolver( const assembler_ptrType& HM_,
             const GetPot& dataFile_,
                   ReferenceMap* refMap_ ) :
             HM( HM_ ),
             DataVelFESpace( HM->velocityFespace()->dof() ), DataPressFESpace( HM->pressureFespace()->dof() ),
             data( new FSIData( dataFile_ ) ),
             refMap( refMap_ ),
             uMeshSize( data->L() / (2*data->Nelements()) ),
             udof( DataVelFESpace.numTotalDof() ), pdof( DataPressFESpace.numTotalDof() ),
             nQuadRho( refMap->qrRho().nbQuadPt() ), nQuadTheta( refMap->qrTheta().nbQuadPt() ),
             numbStep( static_cast<int>( round( (data->T()-data->t0())/data->dt() ) ) ),
             Comm( new Epetra_SerialComm ),
             Map( MapEpetra( (data->mx() + data->mr() + data->mtheta()) * udof + data->mp() * pdof, Comm) ), Map_3D( MapEpetra( nQuadRho * nQuadTheta * ( 3 * udof + pdof ), Comm ) ),
             Map_ALE( udof * nQuadRho * nQuadTheta, Comm ), Map_Wall( MapEpetra( nQuadTheta * udof, Comm ) ), Map_Eta( MapEpetra( udof, Comm ) ),
             block_row( std::vector<UInt>( data->mx() + data->mr() + data->mtheta(), udof ) ), block_col( std::vector<UInt>( data->mx() + data->mr() + data->mtheta(), udof ) ),
             precPtr( new prec_Type ), belosList( Teuchos::rcp ( new Teuchos::ParameterList ) ),
             SystemMatrix( new matrix_Type( Map ) ), Rhs( new vector_Type( Map, Repeated ) ),
             Solution( new vector_Type( Map, Unique ) ),
             solution_3D( vector_Type( Map_3D, Unique ) ), solution_3DOld( vector_Type( Map_3D, Unique ) ),
             urWall( vector_Type( Map_Wall, Unique ) ), urWallOld( vector_Type( Map_Wall, Unique ) ),
             wr( vector_Type( Map_ALE, Unique ) ),
             etar( vector_Type( Map_Eta, Unique ) ), etarOld( vector_Type( Map_Eta, Unique ) ),
             exporterVtk( HiModExporterVtk( *(HM->modalspace()), uMeshSize, data->Nelements(), 1 ) )
  {
    refMap->evaluateMapFSI( data->t0() );

    block_row.resize( data->mx() + data->mr() + data->mtheta() + data->mp(), pdof );
    block_col.resize( data->mx() + data->mr() + data->mtheta() + data->mp(), pdof );

    SystemMatrix->setBlockStructure( block_row, block_col );
    SystemMatrix->zero();

    *Rhs *= 0.0;
    Rhs->setBlockStructure( block_row );

    precPtr->setDataFromGetPot( dataFile_, "prec" );
    belosList = Teuchos::getParametersFromXmlFile( "data/SolverParamList.xml" );
    linearSolver.setCommunicator( Comm );
    linearSolver.setParameters( *belosList );
    linearSolver.setPreconditioner( precPtr );

    setGrid();
    setContainers3D();

    wr *= 0;

    etar *= 0;
    etarOld *= 0;

    *Solution *= 0;

    solution_3DOld = HM->evaluateForce3DGridFSI( data->ux0(), data->ur0(), data->utheta0(), 0, grid );
    urWallOld = HM->evaluateInitialVelocityWallFSI( data->ur0(), grid, data->R() );
  }

  void FSISolver::solveSystem(const Real& t)
  {
    vector_Type f = HM->evaluateForce3DGridFSI( data->fx(), data->fr(), data->ftheta(), t, grid );

    HM->addStokesProblemFSI( SystemMatrix, data->mu(), data->rho_f(), data->dt(), *refMap, t );
    HM->addALEProblemFSI( SystemMatrix, data->rho_f(), wr, *refMap, t );
    HM->addWallProblemFSI( SystemMatrix, data->rho_s(), data->h_s(), data->dt(), data->E(), data->csi(), data->R(), *refMap, t );
    HM->addrhsFSI( Rhs, f, solution_3DOld, data->rho_f(), data->dt() );
    HM->addWallrhsFSI( Rhs, urWallOld, etarOld, data->rho_s(), data->h_s(), data->dt(), data->E(), data->csi(), data->R() );
    HM->addBCxInOutFSI( Rhs, data->p1()(t), data->p2()(t) );
    HM->addBCrInOutFSI( SystemMatrix, Rhs );
    HM->addBCthetaInOutFSI( SystemMatrix, Rhs );
    SystemMatrix->globalAssemble();
    Rhs->globalAssemble();

    linearSolver.setOperator( SystemMatrix );
    linearSolver.setRightHandSide( Rhs );
    linearSolver.solve( Solution );

    solution_3D = HM->evaluateBase3DGridFSI( *Solution );
    urWall = HM->evaluateBaseWallGridFSI( *Solution );

  }

  void FSISolver::expandSolution ()
  {
    for ( UInt s( 0 ); s != udof; ++s ) // on FE nodes
    {
      for ( UInt ntheta( 0 ); ntheta != nQuadTheta; ++ntheta ) // on quadrature node nquadTheta
      {
        for ( UInt nrho( 0 ); nrho != nQuadRho; ++nrho ) // on quadrature node nqadRho
        {
          ux[s][nrho][ntheta] = solution_3D[ xcoord2index( s, nrho, ntheta ) ];
          ur[s][nrho][ntheta] = solution_3D[ rcoord2index( s, nrho, ntheta ) ];
          utheta[s][nrho][ntheta] = solution_3D[ thetacoord2index( s, nrho, ntheta ) ];
          uy[s][nrho][ntheta] = ur[s][nrho][ntheta]*std::cos(std::get<2>(grid[s][nrho][ntheta])) - utheta[s][nrho][ntheta]*std::sin(std::get<2>(grid[s][nrho][ntheta]));
          uz[s][nrho][ntheta] = ur[s][nrho][ntheta]*std::sin(std::get<2>(grid[s][nrho][ntheta])) + utheta[s][nrho][ntheta]*std::cos(std::get<2>(grid[s][nrho][ntheta]));
          p[s][nrho][ntheta] = solution_3D[ pcoord2index( s, nrho, ntheta ) ];
        }
      }
    }
  }

  void FSISolver::computeALEVelocity()
  {
    for (UInt i = 0; i < udof; i++)
      for (UInt j = 0; j < nQuadRho; j++)
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          wr[coord2index(i,j,k)] = ( std::get<1>(grid[i][j][k]) - std::get<1>(gridOld[i][j][k]) ) / data->dt();
        }
  }

  void FSISolver::computeDisplacement()
  {
    Real sum;
    for (UInt i = 0; i < udof; i++)
    {
      sum = 0;
      etarOld[i] = etar[i];
      for (UInt k = 0; k < nQuadTheta; k++)
      {
        sum = etarOld[i] + data->dt() * urWall[coord2indexWall(i,k)];
      }
      etar[i] = sum / nQuadTheta;
    }
  }

  void FSISolver::updateMap( const Real& t )
  {
    std::vector<Real> r( udof, data->R() );
    for (UInt i = 0; i < udof; i++)
    {
      r[i] += etar[i];
    }
    refMap->setRadius( r );
    refMap->evaluateMapFSI( t );
  }

  void FSISolver::updateGrid()
  {
    for (UInt i = 0; i < udof; i++)
      for (UInt j = 0; j < nQuadRho; j++)
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          gridOld[i][j][k] = grid[i][j][k];
          std::get<1>(grid[i][j][k]) = std::get<1>(gridOld[i][j][k]) * (( data->R() + etar[i] ) / (data->R() + etarOld[i] ));
        }

    setCartesianGrid();
  }

  void FSISolver::setGrid()
  {
    gridOld.resize(udof);
    cartesianGrid.resize(udof);
    for (UInt i = 0; i < udof; i++)
    {
      gridOld[i].resize(nQuadRho);
      cartesianGrid[i].resize(nQuadRho);
      for (UInt j = 0; j < nQuadRho; j++)
      {
        gridOld[i][j].resize(nQuadTheta);
        cartesianGrid[i][j].resize(nQuadTheta);
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          std::get<0>(gridOld[i][j][k]) = refMap->nodes()[i];
          std::get<1>(gridOld[i][j][k]) = refMap->qrRho().quadPointCoor(j,0) * data->R();
          std::get<2>(gridOld[i][j][k]) = refMap->qrTheta().quadPointCoor(k,0) * data->theta();
        }
      }
    }

    grid = gridOld;
    setCartesianGrid();
  }

  void FSISolver::setCartesianGrid()
  {
    for (UInt i = 0; i < udof; i++)
    {
      for (UInt j = 0; j < nQuadRho; j++)
      {
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          std::get<0>(cartesianGrid[i][j][k]) = std::get<0>(grid[i][j][k]);
          std::get<1>(cartesianGrid[i][j][k]) = std::get<1>(grid[i][j][k]) * std::cos(std::get<2>(grid[i][j][k]));
          std::get<2>(cartesianGrid[i][j][k]) = std::get<1>(grid[i][j][k]) * std::sin(std::get<2>(grid[i][j][k]));
        }
      }
    }
  }

  void FSISolver::setContainers3D()
  {
    ux.resize(udof);
    ur.resize(udof);
    utheta.resize(udof);
    uy.resize(udof);
    uz.resize(udof);
    p.resize(udof);
    for (UInt i = 0; i < udof; i++)
    {
      ux[i].resize(nQuadRho);
      ur[i].resize(nQuadRho);
      utheta[i].resize(nQuadRho);
      uy[i].resize(nQuadRho);
      uz[i].resize(nQuadRho);
      p[i].resize(nQuadRho);
      for (UInt j = 0; j < nQuadRho; j++)
      {
        ux[i][j].resize(nQuadTheta);
        ur[i][j].resize(nQuadTheta);
        utheta[i][j].resize(nQuadTheta);
        uy[i][j].resize(nQuadTheta);
        uz[i][j].resize(nQuadTheta);
        p[i][j].resize(nQuadTheta);
      }
    }
  }

  void FSISolver::solve()
  {
    Real dt = data->dt();
    Real t = data->t0() - dt;
    for (UInt iter = 1; iter <= numbStep; ++iter)
    {
      t += dt;
      solveSystem( t );
      expandSolution();
      computeDisplacement();
      updateGrid();
      updateMap( t );
      computeALEVelocity();
      exporterVtk.writeSolution( std::string( "output/Solution" ), solution_3D, iter, 1 );
      save( t, iter );
    }
  }

  void FSISolver::save( const Real& t, const UInt& iter ) const
  {
    std::ofstream solution( "output/Iteration" + std::to_string(iter) + ".txt" );

    /*
    solution <<"SOLUTION OF THE FLUID STRUCTURE INTERACTION PROBLEM: t = " <<t <<" s\n" <<std::endl;

    for (UInt i = 0; i < udof; i++)
    {
      for (UInt j = 0; j < nQuadRho; j++)
      {
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          solution <<"Coordinates: x = " <<std::get<0>(cartesianGrid[i][j][k])
                              <<", y = " <<std::get<1>(cartesianGrid[i][j][k])
                              <<", z = " <<std::get<2>(cartesianGrid[i][j][k])
                             <<", (r = " <<std::get<1>(grid[i][j][k])
                          <<", theta = " <<std::get<2>(grid[i][j][k]) <<")" <<std::endl;
          solution <<"Solution: ux = " <<ux[i][j][k]
                           <<", uy = " <<uy[i][j][k]
                           <<", uz = " <<uz[i][j][k]
                           <<", (ur = " <<ur[i][j][k]
                        <<", utheta = " <<utheta[i][j][k] <<")"
                             <<", p = " <<p[i][j][k] <<std::endl;
          solution <<std::endl;
        }
      }
    }

    solution <<"Radial displacemet of the structure:\n" <<std::endl;
    for (UInt i = 0; i < udof; i++)
    {
      solution <<"x = " <<std::get<0>(grid[i][0][0]) <<"\t" <<"etar = " <<etar[i] <<std::endl;
    }
	*/

	for (UInt i = 0; i < udof; i++)
    {
      for (UInt j = 0; j < nQuadRho; j++)
      {
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          solution <<std::get<0>(cartesianGrid[i][j][k]) <<" "
        		   <<std::get<1>(cartesianGrid[i][j][k]) <<" "
                   <<std::get<2>(cartesianGrid[i][j][k]) <<" "
                   <<std::get<1>(grid[i][j][k]) <<" "
                   <<std::get<2>(grid[i][j][k]) <<" "
                   <<ux[i][j][k] <<" "
                   <<uy[i][j][k] <<" "
                   <<uz[i][j][k] <<" "
                   <<ur[i][j][k] <<" "
                   <<utheta[i][j][k] <<" "
                   <<p[i][j][k]  <<" "
                   <<std::endl;
        }
      }
    }
    for (UInt i = 0; i < udof; i++)
    {
      solution <<std::get<0>(cartesianGrid[i][0][0]) <<" "
	           <<etar[i] <<" "
	           <<std::endl;
    }

    solution.close();
  }

}
