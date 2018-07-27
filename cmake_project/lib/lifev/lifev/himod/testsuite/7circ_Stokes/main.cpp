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
#include <lifev/himod/tools/utilityFunctions.hpp>
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
        const UInt Nelements( dataFile( "mesh/num_elements", 10 ) );
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
        Real R = dataFile( "himod/rho", 1. );
        Real dR = dataFile( "himod/drho", 0 );
        Real Theta = dataFile( "himod/theta", 2. * M_PI );

		if (verbose)
        {
            std::cout << "Map definition" << std::endl;
        }

		ReferenceMap	refMap( Jr, Jtheta, Dr, Dtheta, Drtheta, Dthetar, Jacobian, inverseRhat );
		
        if (verbose)
        {
	        std::cout << "-- done! --" << std::endl;
            std::cout << "ModalBasis definition" << std::endl;
        }

        // creation of the modal space
        QuadratureRule	quadRho = refMap.qrRho();
        QuadratureRule	quadTheta = refMap.qrTheta();
        boost::shared_ptr< NSModalSpaceCircular > MB( new NSModalSpaceCircular( R, dR, Theta, mx, mr, mtheta, mp, &refMap, &quadRho, &quadTheta ) );

        if (verbose)
        {
            std::cout << "Adding BC on the wall" << std::endl;
        }

        // applying lateral homogeneous boundary conditions (for the pressure Neumann conditions are enforced automatically)
		// mu * u' + chi * u = 0
        Real mux = 0;
        Real chix = 1;
        Real mur = 0;
        Real chir = 1;
        Real mutheta = 0;
        Real chitheta = 1;
        
		MB->addSliceBC( "dir", mux, chix,
						 "dir", mur, chir,
						 "dir", mutheta, chitheta );
				
        // actually compute the basis
        MB->evaluateBasis();
               
        if (verbose)
        {
            // Some info, from the ModalSpace
            MB->showMe();
	        std::cout << "-- done! --" << std::endl;
            std::cout << "NSHiModAssembler definition" << std::endl;
        }
        // we mix the FEM space with the ModalSpace, creating a HiModSpace
        NSHiModAssembler< mesh_Type , matrix_Type , vector_Type, Circular > HM( uSpace, pSpace, MB, Comm );

        if (verbose)
        {
            std::cout << "Creation of mapepetra for vector structured" << std::endl;
        }

        // we define a Map that will define the vector on the himod space
        // mx + mr + mtheta blocks of monodimensional fem u-functions + mp blocks of monodimensional fem p-functions
        MapEpetra Map( ( mx + mr + mtheta ) * DataVelFESpace.numTotalDof() + mp * DataPressFESpace.numTotalDof(), Comm );
        MapEpetra pressureMap( mp * DataPressFESpace.numTotalDof(), Comm );

        //Using the map it is possible to define the system matrix
        boost::shared_ptr<matrix_Type> systemMatrix( new matrix_Type( Map ) );
        boost::shared_ptr<matrix_Type> massPressMatrix( new matrix_Type( pressureMap ) );
        // boost::shared_ptr<matrix_Type> systemMatrix( new matrix_Type( uSpace->map() | uSpace->map() | uSpace->map() | pSpace->map() ) );

        std::vector<UInt> block_pressRow( mp, DataPressFESpace.numTotalDof() );
        std::vector<UInt> block_pressCol( mp, DataPressFESpace.numTotalDof() );
        massPressMatrix->setBlockStructure( block_pressRow, block_pressCol );
        *massPressMatrix *= 0.0;
        
        //We set here the block structure of the matrix according to the theory
        std::vector<UInt> block_row( MB->mx() + MB->mr() + MB->mtheta(), DataVelFESpace.numTotalDof() );
        block_row.resize( MB->mx() + MB->mr() + MB->mtheta() + MB->mp(), DataPressFESpace.numTotalDof() );

        std::vector<UInt> block_col( MB->mx() + MB->mr() + MB->mtheta(), DataVelFESpace.numTotalDof() );
        block_col.resize( MB->mx() + MB->mr() + MB->mtheta() + MB->mp(), DataPressFESpace.numTotalDof() );

        systemMatrix->setBlockStructure( block_row, block_col );
        *systemMatrix *= 0.0;

		// also the rhs should have a block structure
		boost::shared_ptr<vector_Type> rhs( new vector_Type( Map, Repeated ) );
	   	*rhs *= 0.0;
   		rhs->setBlockStructure( block_row );
   	    
   		// here we create some vectors to build the force term projected on the modal space and then interpolated on the 1D fem structure
   		boost::shared_ptr<vector_Type> f_interpolated( new vector_Type( Map, Repeated ) );
   		*f_interpolated *= 0.0;
   		f_interpolated->setBlockStructure( block_row );
   		   
		boost::shared_ptr<vector_Type> u_old3D_interpolated( new vector_Type( Map, Repeated ) );
   		*u_old3D_interpolated *= 0.0;
   		u_old3D_interpolated->setBlockStructure( block_row );
	
        // We want to evaluate the solution on a 3D grid made by the FEM grid times the quadrature grid        
        vector_Type	u_old( Map, Unique );
		u_old *= 0.0;
        vector_Type	u_kMinusOne( Map, Unique );
		u_kMinusOne *= 0.0;
	
	    
   		boost::shared_ptr<vector_Type> solution;
   		solution.reset ( new vector_Type ( Map, Unique ) );

        UInt nquadRho = MB->qrRho()->nbQuadPt();
        UInt nquadTheta = MB->qrTheta()->nbQuadPt();

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
	
		LinearSolver linearSolver;
		linearSolver.setCommunicator( Comm );
   		linearSolver.setParameters( *belosList );
   		linearSolver.setPreconditioner( precPtr );
   		
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
		
		    MB->map()->evaluateMap( t );
		    
		    if( verbose )
		    {
			std::cout<<"Assembling rhs"<<std::endl;
		    }
		    // Notes that this is not the only way to assemble the rhs, you can use the functor, for more details look test_addrhs.
    		    LifeChrono rhschrono;
    		    rhschrono.start();
    		
		    // Reconstruct the 3d (old) solution and integrate it in order to update the rhs (unsteady case)
		    u_old3D = HM.evaluateBase3DGrid( u_old );
		    u_old3D *= alpha;
	
			*u_old3D_interpolated *= 0;
			*f_interpolated *= 0;
    	    HM.interpolate( u_old3D, u_old3D_interpolated );
    	    HM.interpolate( fxrho, frrho, fthetarho, f_interpolated, t );
    	    
    	    // now that f is actually defined on the himodspace it is possible to add it to the rhs
    	    HM.addrhs( rhs, f_interpolated, u_old3D_interpolated );
    	    rhschrono.stop();
		
    	    if( verbose )
    	    {
    	        std::cout << "Time of assemble rhs : " << rhschrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
				std::cout << "Assembling System Matrix" << std::endl;
    	    }
    		    
    	    // Real test = 1;
    	    // UInt PFiter = 0;
    		    
    	    // while( test > tol && PFiter < PFmax )
		    // {
		    // ++PFiter
		    // evaluate u_kMinusOne
		
		    // we start a chrono to understand the time needed for the assembling
		    LifeChrono himodchrono;
		    himodchrono.start();
		    // Fill the matrix in
    	    HM.addStokesProblem( systemMatrix, nu, refMap, t, alpha );
    	    HM.pressureMassMatrix( massPressMatrix, refMap );
    	    massPressMatrix->spy( "massPressMatrix" );

    	    himodchrono.stop();
	
    	    if (verbose)
    	    {
    	        std::cout << "Time of assemble : " << himodchrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
    	        std::cout << "Adding inflow and outflow BC" << std::endl;
    	    }
   		    // we add dirichlet bc on the inflow boundary, via an algebraic penalization technique
   		    HM.addBC( systemMatrix, rhs, BC, t );
			       
   		    if( verbose )
   		    {
		        std::cout << "-- done! --" << std::endl;
   		        std::cout << "Closing rhs and the matrix" << std::endl;
   		    }
   		    rhs->globalAssemble();
   		    systemMatrix->globalAssemble();
				
   		    // Some export for post processing, might be interesting to explore the sparsity pattern of the matrix
   		    std::stringstream fileIter;
		    std::string filename;
		    fileIter << t_iter;

		    filename = "SystemMatrix" + fileIter.str();
   		    systemMatrix->spy( filename );
   		    filename.erase();

		    filename = "rhs" + fileIter.str();
   		    rhs->spy( filename );
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
    		    
   		    linearSolver.setOperator( systemMatrix );
   		    linearSolver.setRightHandSide( rhs );
   		    linearSolver.solve( solution );
	
   		    solvechrono.stop();
	
   		    if( verbose )
   		    {
   		        std::cout << "Time to solve the problem: " << solvechrono.diff() << std::endl;
		        std::cout << "-- done! --" << std::endl;
   		    }
				
			u_old = *solution;
			
   		    if( verbose )
   		    {
   		        std::cout << "Starting post processing on the solution" << std::endl;
   		    }

			filename = "ModalSolution" + fileIter.str();
			solution->spy( filename );
			filename.erase();

		    // postprocessing on the solution
		    // evaluation of the solution on the 3D grid through modal expansion
		    vector_Type solution_3D( Map_3D, Unique );
		    solution_3D = HM.evaluateBase3DGrid( *solution );

		    filename = "Solution3D" + fileIter.str();
			solution_3D.spy( filename );
			filename.erase();
			fileIter.str( "" );

		    vector_Type ues_3D( Map_3D, Unique );
		    //This is the same but from an analytical function
		    ues_3D = HM.evaluateBase3DGrid( x_ues, r_ues, theta_ues, p_es, t );
	
		    vector_Type err_3D( Map_3D, Unique );
		    err_3D = solution_3D;
		    err_3D += ues_3D * (-1);
			    
		    Real normrho = 1.0 / MB->Rho();
		    Real normtheta = 1.0 / sqrt(2. * M_PI );
	
		    // now that approximate solution and exact solution are evaluate on the grid it is easy to compute the L2 norm
		    Real normVelErr = HM.normL2( err_3D, "velocity" );
		    Real normVelSol = HM.normL2( ues_3D, "velocity" );
	
		    Real normPressErr = HM.normL2( err_3D, "pressure" );
		    Real normPressSol = HM.normL2( ues_3D, "pressure" );
	
		    // now we show how it is possible to export, in VTK format, the himod solution on a self-defined structured grid
	
		    std::cout << "---------------- Compute errors --------------------" << std::endl;
    	    std::cout << "L2 velocity error norm: " << normVelErr << std::endl;
    	    std::cout << "L2 velocity ues norm: " << normVelSol << std::endl;
  		    std::cout << "L2 velocity relative error norm: " << normVelErr / normVelSol << std::endl;
		    std::cout << "\t \t			 --------------------					" << std::endl;
    	    std::cout << "L2 pressure error norm: " << normPressErr << std::endl;
    	    std::cout << "L2 pressure ues norm: " << normPressSol << std::endl;
   		    std::cout << "L2 pressure relative error norm: " << normPressErr / normPressSol << std::endl;
		    std::cout << "-----------------------------------------------------	" << std::endl;
		    
		    HiModExporterVtk exporterVtk( *MB,
		                              uMeshSize, Nelements );
  		    exporterVtk.writeSolution( std::string( "SteadyStokes" ), solution_3D, t_iter );
  		    exporterVtk.writeSolution( std::string( "Poiseuille" ), ues_3D, t_iter );
		    
/*		    std::ofstream velocityVtk("velocity3D.vtk");
		    std::ofstream pressureVtk("pressure3D.vtk");
        	std::ofstream VexactVtk("ExactVelocity3D.vtk");
        	std::ofstream PexactVtk("ExactPressure3D.vtk");
		
			velocityVtk << "# vtk DataFile Version 3.1"<<endl;
			velocityVtk << "This file contains the evaluations of the velocity field for a Steady Stokes problem"<<endl;
			velocityVtk << "ASCII"<<endl;
			velocityVtk << "DATASET UNSTRUCTURED_GRID"<<endl;
			velocityVtk << endl;
			velocityVtk << "POINTS "<<nquadRho*nquadTheta*udof<<" FLOAT"<<endl;
			
			pressureVtk << "# vtk DataFile Version 3.1"<<endl;
			pressureVtk << "This file contains the evaluations of the pressure field for a Steady Stokes problem"<<endl;
			pressureVtk << "ASCII"<<endl;
			pressureVtk << "DATASET UNSTRUCTURED_GRID"<<endl;
			pressureVtk << endl;
			pressureVtk << "POINTS "<<nquadRho*nquadTheta*pdof<<" FLOAT"<<endl;
			
			VexactVtk << "# vtk DataFile Version 3.1"<<endl;
			VexactVtk << "This file contains the evaluations of the exact velocity field for a Steady Stokes problem"<<endl;
			VexactVtk << "ASCII"<<endl;
			VexactVtk << "DATASET UNSTRUCTURED_GRID"<<endl;
			VexactVtk << endl;
			VexactVtk << "POINTS "<<nquadRho*nquadTheta*udof<<" FLOAT"<<endl;
			
			PexactVtk << "# vtk DataFile Version 3.1"<<endl;
			PexactVtk << "This file contains the evaluations of the exact pressure field for a Steady Stokes problem"<<endl;
			PexactVtk << "ASCII"<<endl;
			PexactVtk << "DATASET UNSTRUCTURED_GRID"<<endl;
			PexactVtk << endl;
			PexactVtk << "POINTS "<<nquadRho*nquadTheta*pdof<<" FLOAT"<<endl;
			
			for( UInt ix = 0; ix!= udof; ++ix)
			{     
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						velocityVtk << ix * lx/2*1/Nelements <<"\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * cos (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * sin (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<< endl;
						
						VexactVtk << ix * lx/2*1/Nelements <<"\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * cos (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * sin (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<< endl;
											
					}
				}
			}
			
			velocityVtk << endl;
			velocityVtk << "POINT_DATA "<<nquadRho*nquadTheta*udof<<endl;
			velocityVtk << "VECTORS velocity FLOAT"<<endl;
			velocityVtk << "LOOKUP_TABLE default"<<endl;
			
			VexactVtk << endl;
			VexactVtk << "POINT_DATA "<<nquadRho*nquadTheta*udof<<endl;
			VexactVtk << "VECTORS velocity FLOAT"<<endl;
			VexactVtk << "LOOKUP_TABLE default"<<endl;
			
			for( UInt ix = 0; ix!= udof; ++ix)
			{     
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						velocityVtk << solution_3D[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						solution_3D[ (udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						solution_3D[ (2*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
						
						VexactVtk << ues_3D[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						ues_3D[ (udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << "\t" <<
						ues_3D[ (2*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
											
					}
				}
			}
			velocityVtk << endl;
			VexactVtk << endl;
			
			for( UInt ix = 0; ix!= pdof; ++ix)
			{     
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						pressureVtk << ix * lx/Nelements <<"\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * cos (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * sin (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<< endl;
						
						PexactVtk << ix * lx/Nelements <<"\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * cos (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<<
						MB -> qrRho() -> quadPointCoor (irho, 0) * R * sin (MB -> qrTheta() -> quadPointCoor (itheta, 0) * Theta ) << "\t"<< endl;
					}
				}
			}
			
			pressureVtk << endl;
			pressureVtk << "POINT_DATA "<<nquadRho*nquadTheta*udof<<endl;
			pressureVtk << "SCALARS pressure FLOAT"<<endl;
			pressureVtk << "LOOKUP_TABLE default"<<endl;
			
			PexactVtk << endl;
			PexactVtk << "POINT_DATA "<<nquadRho*nquadTheta*udof<<endl;
			PexactVtk << "SCALARS pressure FLOAT"<<endl;
			PexactVtk << "LOOKUP_TABLE default"<<endl;
			
			for( UInt ix = 0; ix!= pdof; ++ix)
			{     
				for( UInt itheta = 0; itheta!=nquadTheta; ++itheta) 
				{  
					for( UInt irho = 0; irho!=nquadRho; ++irho)
					{
						pressureVtk << solution_3D[ (3*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
						
						PexactVtk << ues_3D[ (3*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << endl;
					}
				}
			}
		
			pressureVtk << endl;
			PexactVtk << endl;
		
		velocityVtk.close();	
		VexactVtk.close();
		PexactVtk.close();
		pressureVtk.close();*/
		} // for loop on t
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
