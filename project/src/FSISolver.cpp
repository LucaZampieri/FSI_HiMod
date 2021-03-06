#include "include/FSISolver.hpp"

#include <cmath>
#include <fstream>

#define DEBUG 1

#ifndef FSI_ENABLED
  #define FSI_ENABLED
#endif

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
             udof( DataVelFESpace.numTotalDof() ),
             pdof( DataPressFESpace.numTotalDof() ),
             nQuadRho( refMap->qrRho().nbQuadPt() ),
             nQuadTheta( refMap->qrTheta().nbQuadPt() ),
             numbStep( static_cast<int>( round( (data->T()-data->t0())/data->dt() ) ) ),
             Comm( new Epetra_SerialComm ),
             Map( MapEpetra( (data->mx() + data->mr() + data->mtheta()) * udof + data->mp() * pdof, Comm) ),
             Map_3D( MapEpetra( nQuadRho * nQuadTheta * ( 3 * udof + pdof ), Comm ) ),
             //Map_ALE( udof * nQuadRho * nQuadTheta, Comm ),
             Map_Wall( MapEpetra( nQuadTheta * udof, Comm ) ),
             Map_Eta( MapEpetra( udof, Comm ) ),
             block_row( std::vector<UInt>( data->mx() + data->mr() + data->mtheta(), udof ) ),
             block_col( std::vector<UInt>( data->mx() + data->mr() + data->mtheta(), udof ) ),
             precPtr( new prec_Type ),
             belosList( Teuchos::rcp ( new Teuchos::ParameterList ) ),
             SystemMatrix( new matrix_Type( Map ) ),
             Rhs( new vector_Type( Map, Repeated ) ),
             Solution( new vector_Type( Map, Unique ) ),
             solution_3D( vector_Type( Map_3D, Unique ) ),
             solution_3DOld( vector_Type( Map_3D, Unique ) ),
             urWall( vector_Type( Map_Wall, Unique ) ),
             urWallOld( vector_Type( Map_Wall, Unique ) ),
             //wr( vector_Type( Map_ALE, Unique ) ),
             etar( vector_Type( Map_Eta, Unique ) ),
             etarOld( vector_Type( Map_Eta, Unique ) ),
             exporterVtk( HiModExporterVtk( *(HM->modalspace()), uMeshSize, data->Nelements(), 1 ) )
  {
    //data -> printAll(); // print all the values found in datafile
    etar *= 0;
    etarOld *= 0;
    *Solution *= 0;
    *Rhs *= 0.0;


    refMap->evaluateAxialMap( data->L()/(pdof - 1) /* h */, udof, data->Radius(), data->dRadius());

    block_row.resize( data->mx() + data->mr() + data->mtheta() + data->mp(), pdof );
    block_col.resize( data->mx() + data->mr() + data->mtheta() + data->mp(), pdof );

    std::cout << "count of the number of nodes : " << (data->mx() + data->mr() + data->mtheta()) * udof + data->mp() * pdof << std::endl;

    SystemMatrix->setBlockStructure( block_row, block_col );
    SystemMatrix->zero();


    Rhs->setBlockStructure( block_row );

    precPtr->setDataFromGetPot( dataFile_, "prec" );
    belosList = Teuchos::getParametersFromXmlFile( "data/SolverParamList.xml" );
    linearSolver.setCommunicator( Comm );
    linearSolver.setParameters( *belosList );
    linearSolver.setPreconditioner( precPtr );

    setGrid();
    setContainers3D();

    solution_3DOld = HM->evaluateForce3DGridFSI( data->ux0(), data->ur0(), data->utheta0(), 0, grid );
    urWallOld = HM->evaluateInitialVelocityWallFSI( data->ur0(), grid, data->Radius() );
    #ifdef FSI_ENABLED
      std::cout << "\nFSI ENABLED in fsisolver.cpp. \n";
    #endif

  }

  void FSISolver::solveSystem(const Real& t)
  {
    Rhs->setBlockStructure( block_row );
    *Rhs *= 0.0;
    *Solution *= 0.0;
    SystemMatrix->setBlockStructure( block_row, block_col );
    SystemMatrix->zero();


    vector_Type f = HM->evaluateForce3DGridFSI( data->fx(), data->fr(), data->ftheta(), t, grid );

    HM->addStokesProblemFSI( SystemMatrix, data->nu(), data->rho_s(), data->h_s(), data->e(), *refMap, t, data->alpha());

    #ifdef FSI_ENABLED
      HM->addRhsFSI( Rhs, data->alpha(), data->rho_s(), data->h_s(), data->e(), f, solution_3DOld, urWallOld, etar);
    #endif
    HM->addBcFSI( SystemMatrix, Rhs, data->p1()(t), data->p2()(t), data->uIn(), data->utheta0() /* 0 */ ,t );

    SystemMatrix->globalAssemble();
    Rhs->globalAssemble();



    linearSolver.setOperator( SystemMatrix );
    linearSolver.setRightHandSide( Rhs );
    linearSolver.solve( Solution );

    solution_3D = HM->evaluateBase3DGrid( *Solution );
    //urWall = HM->evaluateBaseWallGridFSI( *Solution );
    solution_3DOld = solution_3D;
  }

  void FSISolver::solve()
  {
    Real dt = data->dt();
    Real t = data->t0() ; // -data->dt()
    Real finalTime = data->T();
    for (UInt iter = 1; iter <= numbStep; ++iter)
    {
      t += dt;
      std::cout << "\n----------------------------------------\n";
      std::cout << "Entering iteration: "+to_string(iter)+" / "+to_string(numbStep) << std::endl;
      std::cout << "Time in simulation: "+to_string(t) + " / "+ to_string(finalTime)   << std::endl;
      std::cout << "----------------------------------------\n\n";
      solveSystem( t );
      expandSolution();
      computeDisplacement();
      //updateGrid();
      //updateMap( t );
      //computeALEVelocity();
      exporterVtk.writeSolution( std::string( "output/Solution" ), solution_3D, iter, 1 );
      save( iter );

    }
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

    #ifdef DEBUG
      std::cout << "\n\n\nCheck: \n\n";

      std::cout << "\n\n\nPrinting ux ,ur, utheta: \n\n";
      std::cout << "Printing ux (horizontal): \n\n";
      std::cout << ux[0][0][16]   << "     "
                << ux[1][0][16]   << "     "
                << ux[3][0][16]   << "     "
                << ux[13][0][16]  << "     "
                << ux[25][0][16]  << "     "
                << ux[30][0][16]  << "     "
                <<  std::endl;


      std::cout << "\n\n Printing ux: \n\n";

      for (int j=0 ; j!= udof; ++j )
      {
        std::cout << ux[j][nQuadRho-1][16]   << "     ";
        std::cout << ur[j][nQuadRho-1][16];
        std::cout << std::endl;
        /*
        std::cout << ur[0][i][16]   << "     "
                  << ur[1][i][16]   << "     "
                  << ur[3][i][16]   << "     "
                  << ur[13][i][16]  << "     "
                  << ur[25][i][16]  << "     "
                  << ur[30][i][16]  << "     " */
      }

    #endif
    /*
    std::cout << "\n\n Printing utheta: \n\n\n";
    for (int i=0; i != nQuadTheta; ++i)
    {
      std::cout << utheta[16][16][i] << std::endl;
    }*/

  }

  /*
  void FSISolver::computeALEVelocity()
  {
    for (UInt i = 0; i < udof; i++)
      for (UInt j = 0; j < nQuadRho; j++)
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          wr[coord2index(i,j,k)] = ( std::get<1>(grid[i][j][k]) - std::get<1>(gridOld[i][j][k]) ) / data->dt();
        }
  }
  */ // end comment of computeALEVelocity
  void FSISolver::computeDisplacement()
  {
    Real sum;
    std::cout << "printing etar\n";
    for (UInt i = 0; i < udof; i++)
    {
      sum = 0;
      // etarOld[i] = etar[i];

      /* // useless if axysymmetric
      for (UInt k = 0; k < nQuadTheta; k++)
      {
        // sum = etarOld[i] + data->dt() * urWall[coord2indexWall(i,k)];
        sum +=  ur[i][nQuadRho-1][k];
      }
      etar[i] =  etar[i] + data->dt() * sum / nQuadTheta;
      */
      etar[i] +=  data->dt() * ur[i][nQuadRho-1][16]; // 16 -->
      std::cout << "etar: " <<etar[i] << std::endl;
    }
  }

  /*
  void FSISolver::updateMap( const Real& t )
  {
    std::vector<Real> r( udof, data->Radius() );
    for (UInt i = 0; i < udof; i++)
    {
      r[i] += etar[i];
    }
    refMap->setRadius( r );
    refMap->evaluateAxialMap( t );
  }*/

  /*
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
  */
  void FSISolver::setGrid()
  {
    //gridOld.resize(udof);
    grid.resize(udof);
    cartesianGrid.resize(udof);
    for (UInt i = 0; i < udof; i++)
    {
      grid[i].resize(nQuadRho);
      cartesianGrid[i].resize(nQuadRho);
      for (UInt j = 0; j < nQuadRho; j++)
      {
        grid[i][j].resize(nQuadTheta);
        cartesianGrid[i][j].resize(nQuadTheta);
        for (UInt k = 0; k < nQuadTheta; k++)
        {
          Real x_tmp = data->L()/(2*(pdof - 1)) /* h */ *i; // temporary x to update grid
          std::get<0>(grid[i][j][k]) = x_tmp;
          std::get<1>(grid[i][j][k]) = refMap->qrRho().quadPointCoor(j,0) * data->Radius()(0,x_tmp,0,0,0);
          std::get<2>(grid[i][j][k]) = refMap->qrTheta().quadPointCoor(k,0) * data->theta();
        }
      }
    }

    //grid = gridOld;
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
  { //TODO We could use fixed-size arrays instead of vectors
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

  //void FSISolver::save( const Real& t, const UInt& iter ) const
  void FSISolver::save( const UInt& iter ) const
  {

      std::ofstream solution( "output/Iteration" + std::to_string(iter) + ".txt" );

      for (UInt i = 0; i < udof; i++)
      {
          Real x_tmp = std::get<0>(cartesianGrid[i][0][0]);
          solution <<x_tmp <<" "
                   <<data->Radius()(0,x_tmp,0,0,0) <<" "
	                 <<etar[i] <<" "
	                 <<std::endl;
      }

      solution.close();
  }

}
