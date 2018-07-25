//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file main.cpp
    @brief Tutorial introducing how to use the himod package

    @author Sofia Guzzetti <sofia.guzzetti@gmail.com>
    @date 03-2014

    In this tutorial we explain how to solve Navier-Stokes Equations
    with the Hierarchical Model Reduction method, with
    lateral homogeneous Dirichlet boundary conditions.
 */

// ----------------------------------------------------//
// There is only a serial version of the code so far   //
// ----------------------------------------------------//

#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/lexical_cast.hpp>
#include <string.h>
#include <string>
#include <stdlib.h>

// -----------------------------------------------------//
// These headers are needed for the solution of the     //
// linear system                                        //
// -----------------------------------------------------//
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wwrite-strings"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/filter/GetPot.hpp>
/*
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"
*/
// ------------------------------------------------------//
// Even himod is a 3D solver only a monodimensional grid //
// is needed, we use a structured one                    //
// ------------------------------------------------------//
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

// ------------------------------------------------------ //
// these headers are needed for the data structures       //
// ------------------------------------------------------ //
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

// ------------------------------------------------------ //
// We have to include the headers of the FESpace for      //
// the monodimensional supporting fiber                   //
// ModalSpace deals with all the aspects reagarding the   //
// transversal fiber like the creation of the modal       //
// basis                                                  //
// ------------------------------------------------------ //
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/himod/modalbasis/NSModalSpaceCircular.hpp>
#include <lifev/himod/modalbasis/NSHiModAssembler.hpp>
#include <lifev/himod/tools/ReferenceMap.hpp>
#include <lifev/himod/tools/BCstructure.hpp>
#include <lifev/himod/tools/utilityFunctionsCone.hpp>
#include <lifev/himod/tools/HiModExporterVtk.hpp>

#include <fstream>
#include <sstream>

// We use the LifeV namespace ...
using namespace LifeV;

// and we make some usefull typedef for the mesh ...
typedef RegionMesh<LinearLine>              mesh_Type;

// ... the data structures ...
typedef MatrixEpetraStructured<Real>        matrix_Type;
typedef VectorEpetraStructured              vector_Type;
typedef VectorSmall<3>                      TreDvector_type;

// ... for the preconditioner
typedef LifeV::Preconditioner                      basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
typedef PreconditionerIfpack                prec_Type;
typedef boost::shared_ptr<prec_Type>        precPtr_Type;

// Now we are ready to start ...
int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif
    {

        // needed to properly destroy all objects inside before mpi finalize
#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
        ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "-  Resolution of Navier-Stokes Equations using a Hierarchical Model Reduction approach   -" << std::endl;
        std::cout << "-           				Master Thesis by Sofia Guzzetti           						-" << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;


        // we use getpot to obtain some data given by the user

        // ********** GetPot **********
        GetPot command_line( argc, argv );
        GetPot dataFile( "data" );
        UInt verbose = dataFile("set/verbose", 0);

        // total number of modes used in the transversal fiber
        UInt mx = dataFile( "himod/mx", 10 );
        UInt mr = dataFile( "himod/mr", 10 );
        UInt mtheta = dataFile( "himod/mtheta", 10 );
        UInt mp = dataFile( "himod/mp", 8 );
        
        Real lx( dataFile ("mesh/lx", 1.) );
        const UInt Nelements( dataFile( "mesh/num_elements", 10 ) );
        Real h = lx / ( 2 * Nelements );
		
		// Cynematic viscosity
		Real nu = dataFile( "physics/nu", 0.1 );
		
		// time discretization
		Real t = dataFile( "time/t_In", 0 );
		Real Ntimesteps = dataFile( "time/Nsteps", 100 );
		Real dt = dataFile( "time/dt", 0.01 );
		Real alpha = 0;
		
		// Fixed Point Method ( uncomment in case of Navier-Stokes )
		Real tol = dataFile( "fixedPoint/tol", 1e-4 );
		Real PFmax = dataFile( "fixedPoint/PF_max", 1000 );
	
        //-----------------------------------------------------------------------
        //          SPACES AND PROBLEM
        //----------------------------------------------------------------------

        // Mesh definition, we used a structured mesh for the dominant, supporting fiber
        // the number of elements is user defined with getpot
        if (verbose)
        {
            std::cout << "Mesh definition" << std::endl;
        }

        boost::shared_ptr< mesh_Type > fullMeshPtr( new mesh_Type );
        regularMesh1D( *fullMeshPtr, 0, Nelements, false, dataFile( "mesh/lx", 1. ), 0.0 );
        Real uMeshSize = lx / ( 2 * Nelements );


        if (verbose)
        {
            std::cout << "FESpace definition" << std::endl;
        }

        // We define here the FE spaces for the scalar component of the velocity on the supporting fiber...
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
        						( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, dataFile( "mesh/Poly_type_velocity", "P2" ), 1, Comm) );
        						
		// ...and for the pressure
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > pSpace
        						( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, dataFile( "mesh/Poly_type_pressure", "P1" ), 1, Comm) );

        // display the number of FEM dof
        DOF DataVelFESpace( uSpace->dof() );
        UInt udof = DataVelFESpace.numTotalDof();
        
        DOF DataPressFESpace( pSpace->dof() );
        UInt pdof = DataPressFESpace.numTotalDof();
        
        UInt numdof = ( mx + mr + mtheta ) * udof + mp * pdof;
        if (verbose)
        {
            std::cout << "Velocity (scalar) degrees of freedom FEM = " << udof << std::endl;
            std::cout << "Pressure degrees of freedom FEM = " << pdof << std::endl;
            std::cout << "Total degrees of freedom FEM = " << numdof << std::endl;
            std::cout << "-- done! --" << std::endl;
        }

        // Read the dimension of the transversal fiber from getpot
        // Real R = dataFile( "himod/rho", 1. );
        // Real dR = dataFile( "himod/drho", 0 );
        Real Theta = dataFile( "himod/theta", 2. * M_PI );

		if (verbose)
        {
            std::cout << "Map definition" << std::endl;
        }

		ReferenceMap	refMap( Jr, Jtheta, Dr, Dtheta, Drtheta, Dthetar, Jacobian, inverseRhat, uSpace->map() );
		
        if (verbose)
        {
	        std::cout << "-- done! --" << std::endl;
            std::cout << "ModalBasis definition" << std::endl;
        }

        // creation of the modal space. R and dR are defined in utilityFunctionsCone.hpp
        QuadratureRule	quadRho = refMap.qrRho();
        QuadratureRule	quadTheta = refMap.qrTheta();
        boost::shared_ptr< NSModalSpaceCircular > MBsin( new NSModalSpaceCircular( R, dR, Theta, mx, mr, mtheta, mp, &refMap, &quadRho, &quadTheta ) );
        boost::shared_ptr< NSModalSpaceCircular > MBcos( new NSModalSpaceCircular( R, dR, Theta, mx, mr, mtheta, mp, &refMap, &quadRho, &quadTheta ) );

        if (verbose)
        {
            std::cout << "Creating Uneducated Basis on the transversal section" << std::endl;
        }
		
        // actually compute the basis
        MBcos->evaluateBasis( "ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShiftNatural",
        						cosine );
        						
		MBsin->evaluateBasis( "ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShift",
        						"ChebyshevQuadraticShiftNatural",
        						sine );
               
        if (verbose)
        {
            // Some info, from the ModalSpace
            MBsin->showMe();
	        std::cout << "-- done! --" << std::endl;
            std::cout << "NSHiModAssembler definition" << std::endl;
        }
        // we mix the FEM space with the ModalSpace, creating a HiModSpace
        NSHiModAssembler< mesh_Type , matrix_Type , vector_Type, Circular > HMsin( uSpace, pSpace, MBsin, Comm );
        NSHiModAssembler< mesh_Type , matrix_Type , vector_Type, Circular > HMcos( uSpace, pSpace, MBcos, Comm );

        if (verbose)
        {
            std::cout << "Creation of mapepetra for vector structured" << std::endl;
        }

        // we define a Map that will define the vector on the himod space
        // mx + mr + mtheta blocks of monodimensional fem u-functions + mp blocks of monodimensional fem p-functions
        MapEpetra Map( ( mx + mr + mtheta ) * DataVelFESpace.numTotalDof() + mp * DataPressFESpace.numTotalDof(), Comm );

        //Using the map it is possible to define the system matrix
        boost::shared_ptr<matrix_Type> sinSystemMatrix( new matrix_Type( Map ) );
        boost::shared_ptr<matrix_Type> cosSystemMatrix( new matrix_Type( Map ) );
        // boost::shared_ptr<matrix_Type> systemMatrix( new matrix_Type( uSpace->map() | uSpace->map() | uSpace->map() | pSpace->map() ) );

        //We set here the block structure of the matrix according to the theory
        std::vector<UInt> block_row( mx + mr + mtheta, DataVelFESpace.numTotalDof() );
        block_row.resize( mx + mr + mtheta + mp, DataPressFESpace.numTotalDof() );

        std::vector<UInt> block_col( mx + mr + mtheta, DataVelFESpace.numTotalDof() );
        block_col.resize( mx + mr + mtheta + mp, DataPressFESpace.numTotalDof() );

        sinSystemMatrix->setBlockStructure( block_row, block_col );
        *sinSystemMatrix *= 0.0;
        
        cosSystemMatrix->setBlockStructure( block_row, block_col );
        *cosSystemMatrix *= 0.0;

		// also the rhs should have a block structure
		boost::shared_ptr<vector_Type> sinRhs( new vector_Type( Map, Repeated ) );
	   	*sinRhs *= 0.0;
   		sinRhs->setBlockStructure( block_row );
   		
   		boost::shared_ptr<vector_Type> cosRhs( new vector_Type( Map, Repeated ) );
	   	*cosRhs *= 0.0;
   		cosRhs->setBlockStructure( block_row );
   	    
   		// here we create some vectors to build the force term projected on the modal space and then interpolated on the 1D fem structure
   		boost::shared_ptr<vector_Type> fSin_interpolated( new vector_Type( Map, Repeated ) );
   		*fSin_interpolated *= 0.0;
   		fSin_interpolated->setBlockStructure( block_row );
   		
   		boost::shared_ptr<vector_Type> fCos_interpolated( new vector_Type( Map, Repeated ) );
   		*fCos_interpolated *= 0.0;
   		fCos_interpolated->setBlockStructure( block_row );
   		   
		boost::shared_ptr<vector_Type> u_old3Dsin_interpolated( new vector_Type( Map, Repeated ) );
   		*u_old3Dsin_interpolated *= 0.0;
   		u_old3Dsin_interpolated->setBlockStructure( block_row );

		boost::shared_ptr<vector_Type> u_old3Dcos_interpolated( new vector_Type( Map, Repeated ) );
   		*u_old3Dcos_interpolated *= 0.0;
   		u_old3Dcos_interpolated->setBlockStructure( block_row );

	
        // We want to evaluate the solution on a 3D grid made by the FEM grid times the quadrature grid        
        vector_Type	u_oldSin( Map, Unique );
		u_oldSin *= 0.0;
		vector_Type	u_oldCos( Map, Unique );
		u_oldCos *= 0.0;
	    
   		boost::shared_ptr<vector_Type> sinSolution;
   		sinSolution.reset ( new vector_Type ( Map, Unique ) );
   		boost::shared_ptr<vector_Type> cosSolution;
   		cosSolution.reset ( new vector_Type ( Map, Unique ) );

        UInt nquadRho = MBsin->qrRho()->nbQuadPt();
        UInt nquadTheta = MBsin->qrTheta()->nbQuadPt();

        MapEpetra Map_3D( nquadRho * nquadTheta * ( 3 * udof + pdof ), Comm );
        vector_Type u_old3D( Map_3D, Unique );
        
        if( verbose )
        {
	        std::cout << "-- done! --" << std::endl;
	        std::cout<< " Defining inflow and outflow BC " <<std::endl;
        }
        // Boundary conditions
        BCScalarData	xInflow( dir, xInf );
        BCScalarData	rInflow( dir, rInf );
        BCScalarData	thetaInflow( dir, thetaInf );
        
        BCScalarData	xOutflow( neu, xOut );
        BCScalarData	rOutflow( dir, rOut );
        BCScalarData	thetaOutflow( dir, thetaOut );
        
        BCdata		BC( xInflow, rInflow, thetaInflow,
			    xOutflow, rOutflow, thetaOutflow );
        					
        // Solver initialization
        if( verbose )
    	{
	        std::cout << "-- done! --" << std::endl;
   	        std::cout << std::endl << "[Solvers initialization]" << std::endl;
   		}
	   	boost::shared_ptr<prec_Type> precPtr( new prec_Type );
   		precPtr->setDataFromGetPot( dataFile, "prec" );
	
   		if( verbose )
   		{
   		    std::cout << "Setting up LinearSolver (Belos)... " << std::flush<<std::endl;
   		}
   		Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
   		belosList = Teuchos::getParametersFromXmlFile( "SolverParamList.xml" );
	
		LinearSolver sinLinearSolver;
		sinLinearSolver.setCommunicator( Comm );
   		sinLinearSolver.setParameters( *belosList );
   		sinLinearSolver.setPreconditioner( precPtr );
   		
   		LinearSolver cosLinearSolver;
		cosLinearSolver.setCommunicator( Comm );
   		cosLinearSolver.setParameters( *belosList );
   		cosLinearSolver.setPreconditioner( precPtr );
   		
   		if( verbose )
   		{
		        std::cout << "-- done! --" << std::endl;
		        std::cout<< " -- Ready to go! --" << std::endl;
   		}
	
   		for( UInt t_iter = 0; t_iter != Ntimesteps; ++t_iter )
		{
		 	t = t + dt;
		 	std::cout << "-------------------------------------------" << std::endl;
		    std::cout << "Time step " << t_iter + 1 << " out of " << Ntimesteps << std::endl;
		    std::cout << "t = " << t << std::endl;
		
		    // NOT ANYMORE because the map is time independent now
		    //MBsin->map()->evaluateMap( t );
		    //MBcos->map()->evaluateMap( t );
            MBsin->map()->evaluateAxialMap( h, udof, R, dR );
            MBcos->map()->evaluateAxialMap( h, udof, R, dR );
		    
		    if( verbose )
		    {
			std::cout<<"Assembling rhs"<<std::endl;
		    }
		    // Notes that this is not the only way to assemble the rhs, you can use the functor, for more details look test_addrhs.
    		    LifeChrono rhschrono;
    		    rhschrono.start();
    		
		    // Reconstruct the 3d (old) solution and integrate it in order to update the rhs (unsteady case)
		    u_old3D *= alpha;

			*u_old3Dsin_interpolated *= 0;
			*u_old3Dcos_interpolated *= 0;
			*fSin_interpolated *= 0;
			*fCos_interpolated *= 0;
    	    HMsin.interpolate( u_old3D, u_old3Dsin_interpolated, 1 );
    	    HMsin.interpolate( fxrho, frrho, fthetarho, fSin_interpolated, t, 1 );
    	    HMcos.interpolate( u_old3D, u_old3Dcos_interpolated, 1 );
    	    HMcos.interpolate( fxrho, frrho, fthetarho, fCos_interpolated, t, 1 );
    	    
    	    // now that f is actually defined on the himodspace it is possible to add it to the rhs
    	    HMsin.addrhs( sinRhs, fSin_interpolated, u_old3Dsin_interpolated );
    	    HMcos.addrhs( cosRhs, fCos_interpolated, u_old3Dcos_interpolated );
    	    rhschrono.stop();
		
    	    if( verbose )
    	    {
    	        std::cout << "Time of assemble rhs : " << rhschrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
				std::cout << "Assembling System Matrix" << std::endl;
    	    }
    		    
		    // we start a chrono to understand the time needed for the assembling
		    LifeChrono himodchrono;
		    himodchrono.start();
		    // Fill the matrix in for non constant radius
    	    HMsin.addStokesProblem( sinSystemMatrix, nu, refMap, t, alpha, 1 );
    	    HMcos.addStokesProblem( cosSystemMatrix, nu, refMap, t, alpha, 1 );
    	    himodchrono.stop();
	
    	    if (verbose)
    	    {
    	        std::cout << "Time of assemble : " << himodchrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
    	        std::cout << "Adding inflow and outflow BC" << std::endl;
    	    }
   		    // we add dirichlet bc on the inflow boundary, via an algebraic penalization technique
   		    HMcos.addBC( cosSystemMatrix, cosRhs, BC, t, cosine, 1 );
   		    HMsin.addBC( sinSystemMatrix, sinRhs, BC, t, sine, 1 );
			       
   		    if( verbose )
   		    {
		        std::cout << "-- done! --" << std::endl;
   		        std::cout << "Closing rhs and the matrix" << std::endl;
   		    }
   		    sinRhs->globalAssemble();
   		    cosRhs->globalAssemble();
   		    sinSystemMatrix->globalAssemble();
   		    cosSystemMatrix->globalAssemble();
				
				
   		    // Some export for post processing, might be interesting to explore the sparsity pattern of the matrix
   		    std::stringstream fileIter;
		    std::string filename;
		    fileIter << t_iter;

		    filename = "sinSystemMatrix" + fileIter.str();
   		    sinSystemMatrix->spy( filename );
   		    filename.erase();

		    filename = "cosSystemMatrix" + fileIter.str();
   		    cosSystemMatrix->spy( filename );
   		    filename.erase();

		    filename = "sinRhs" + fileIter.str();
   		    sinRhs->spy( filename );
   		    filename.erase();
	
		    filename = "cosRhs" + fileIter.str();
   		    cosRhs->spy( filename );
   		    filename.erase();
	
   		    //-----------------------------------------------------------------------
   		    //              SOLVE
   		    //-----------------------------------------------------------------------
	
   		    LifeChrono solvechrono;
   		    solvechrono.start();
	
   		    if( verbose )
   		    {
		        std::cout << "-- done! --" << std::endl;
   		        std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
   		    }
    		    
   		    sinLinearSolver.setOperator( sinSystemMatrix );
   		    sinLinearSolver.setRightHandSide( sinRhs );
   		    sinLinearSolver.solve( sinSolution );
	
   		    cosLinearSolver.setOperator( cosSystemMatrix );
   		    cosLinearSolver.setRightHandSide( cosRhs );
   		    cosLinearSolver.solve( cosSolution );
	
   		    solvechrono.stop();
	
   		    if( verbose )
   		    {
   		        std::cout << "Time to solve the problem: " << solvechrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
   		    }
			
   		    if( verbose )
   		    {
   		        std::cout << "Starting post processing on the solution" << std::endl;
   		    }

			filename = "SinModalSolution" + fileIter.str();
			sinSolution->spy( filename );
			filename.erase();

			filename = "CosModalSolution" + fileIter.str();
			cosSolution->spy( filename );
			filename.erase();

		    // postprocessing on the solution
		    // evaluation of the solution on the 3D grid through modal expansion
		    vector_Type solution_3D( Map_3D, Unique );
		    solution_3D = HMsin.evaluateBase3DGrid( *sinSolution );
		    solution_3D += HMcos.evaluateBase3DGrid( *cosSolution );

			u_old3D = solution_3D;

		    filename = "Solution3D" + fileIter.str();
			solution_3D.spy( filename );
			filename.erase();
			fileIter.str( "" );
	        
	        HiModExporterVtk exporterVtk( *MBsin,
		                              uMeshSize, Nelements, 1 );
  		    exporterVtk.writeSolution( std::string( "SteadyStokes" ), solution_3D, t_iter + 1, 1 );
  		    
  		    /*
		    std::ofstream velocityVtk("velocity3D.vtk");
		    std::ofstream pressureVtk("pressure3D.vtk");
			
			for( UInt ix = 0; ix!= udof; ++ix)
			{     
				Real Radius( R( 0, ix * h, 0, 0, 0 ) );
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						velocityVtk << ix * h <<"\t"<<
						MBsin -> qrRho() -> quadPointCoor (irho, 0) * Radius * cos (MBsin -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MBsin -> qrRho() -> quadPointCoor (irho, 0) * Radius * sin (MBsin -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						solution_3D[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						solution_3D[ (udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						solution_3D[ (2*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
											
					}
				}
			}
			
			h = lx / ( Nelements );
			
			for( UInt ix = 0; ix!= pdof; ++ix)
			{     
				Real Radius( R( 0, ix * h, 0, 0, 0 ) );
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						pressureVtk << ix * lx / Nelements <<"\t"<<
						MBsin -> qrRho() -> quadPointCoor (irho, 0) * Radius * cos (MBsin -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MBsin -> qrRho() -> quadPointCoor (irho, 0) * Radius * sin (MBsin -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						solution_3D[ (3*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
						
					}
				}
			}   
		
		velocityVtk.close();	
		pressureVtk.close();*/
		} // for loop on t
    } //unnamed namespace
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
