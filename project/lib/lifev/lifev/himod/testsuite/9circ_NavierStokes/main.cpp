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
#include <lifev/himod/tools/utilityNSFunctions.hpp>

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
		
		// Cynematic viscosity
		Real nu = dataFile( "fluid/physics/nu", 0.1 );
		
		// time discretization
		Real t = dataFile( "time/t_In", 0 );
		Real Ntimesteps = dataFile( "time/Nsteps", 100 );
		Real dt = dataFile( "time/dt", 0.01 );
		Real alpha = 1./dt;
		
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

		// --------------------------------		ALGEBRAIC STRUCTURES	----------------------------------------------

        UInt nquadRho = MB->qrRho()->nbQuadPt();
        UInt nquadTheta = MB->qrTheta()->nbQuadPt();

        if (verbose)
        {
            std::cout << "Creation of mapepetra for vector structured" << std::endl;
        }

        // we define a Map that will define the vector on the himod space
        // mx + mr + mtheta blocks of monodimensional fem u-functions + mp blocks of monodimensional fem p-functions
        MapEpetra Map( ( mx + mr + mtheta ) * DataVelFESpace.numTotalDof() + mp * DataPressFESpace.numTotalDof(), Comm );

        // *******************		System matrix	***************************
        boost::shared_ptr<matrix_Type> systemMatrix( new matrix_Type( Map ) );

        //We set here the block structure of the matrix according to the theory
        std::vector<UInt> block_row( MB->mx() + MB->mr() + MB->mtheta(), DataVelFESpace.numTotalDof() );
        block_row.resize( MB->mx() + MB->mr() + MB->mtheta() + MB->mp(), DataPressFESpace.numTotalDof() );

        std::vector<UInt> block_col( MB->mx() + MB->mr() + MB->mtheta(), DataVelFESpace.numTotalDof() );
        block_col.resize( MB->mx() + MB->mr() + MB->mtheta() + MB->mp(), DataPressFESpace.numTotalDof() );

        systemMatrix->setBlockStructure( block_row, block_col );
        *systemMatrix *= 0.0;

		// *******************		rhs		***************************
		boost::shared_ptr<vector_Type> rhs( new vector_Type( Map, Repeated ) );
	   	*rhs *= 0.0;
   		rhs->setBlockStructure( block_row );
   	    
   		// *******************		modal coefficients of the force term * r	***************************
   		boost::shared_ptr<vector_Type> f_interpolated( new vector_Type( Map, Repeated ) );
   		*f_interpolated *= 0.0;
   		f_interpolated->setBlockStructure( block_row );
   		   
        // *******************		solution	***************************
        // update solution in TIME
        vector_Type	u_old( Map, Unique );
		u_old *= 0.0;
        // modal coefficients of u_old * r (to build the force term)
		boost::shared_ptr<vector_Type> u_old3D_interpolated( new vector_Type( Map, Repeated ) );
   		*u_old3D_interpolated *= 0.0;
   		u_old3D_interpolated->setBlockStructure( block_row );
		// 3D expansion of the solution from modal coefficients
		MapEpetra Map_3D( nquadRho * nquadTheta * ( 3 * udof + pdof ), Comm );
        vector_Type u_old3D( Map_3D, Unique );
        u_old *= 0.0;
        
        // update FIXED POINT solution
        vector_Type	u_kMinusOne( Map, Unique );
		u_kMinusOne *= 0.0;
		// 3D expansion of the fixed point (k-1)-solution from modal coefficients
		vector_Type u_kMinusOne3D( Map_3D, Unique );
   		u_kMinusOne3D *= 0.0;
	
	    // modal solution
   		boost::shared_ptr<vector_Type> solution;
   		solution.reset ( new vector_Type ( Map, Unique ) );
   		
   		vector_Type solution_3D( Map_3D, Unique );
   		solution_3D *= 0.0;
        
        // *******************		error	***************************
	    vector_Type err_3D( Map_3D, Unique );
	    err_3D *= 0.0;
        
		// ---------------------------------------------------------------------------------------------------------------
        
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
    	    HM.interpolate( u_old3D, u_old3D_interpolated ); // IMPORTANT: in this case u_old3D*r is interpolated implicitly
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
    		    
    	    Real test = 1;
    	    UInt PFiter = 0;
    	    
    	    std::stringstream fileIter;
			std::string filename;
    		
    		LifeChrono solvechrono;
	   		solvechrono.start();
	
    	    while( test > tol && PFiter < PFmax )
		    {
				++PFiter;
				u_kMinusOne3D = HM.evaluateBase3DGrid( u_kMinusOne );
		
				// we start a chrono to understand the time needed for the assembling
				LifeChrono himodchrono;
				himodchrono.start();
				// Fill the matrix in
			    HM.addNavierStokesProblem( systemMatrix, nu, refMap, t, alpha, u_kMinusOne3D );
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
	
	   		    if( verbose )
	   		    {
				    std::cout << "-- done! --" << std::endl;
	   		        std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
	   		    }
				    
	   		    linearSolver.setOperator( systemMatrix );
	   		    linearSolver.setRightHandSide( rhs );
	   		    linearSolver.solve( solution );
	
	   		    // compute test quantity
	   		    solution_3D = HM.evaluateBase3DGrid( *solution );
	   		    
			    err_3D = u_kMinusOne3D;
			    err_3D +=  solution_3D * (-1);
			    
			    test = HM.normL2( err_3D, "velocity" ) / HM.normL2( solution_3D, "velocity" );
			    
			    // update FIXED POINT solution
			    u_kMinusOne = *solution;
	   		    
   		    } // end while
    		
    		std::cout << "Fixed Point iterations: " << PFiter <<std::endl;
   		    if( PFiter == PFmax )
   		    {
   		    		std::cout << "WARNING: maximum number of Fixed Point iterations reached." << std::endl;
   		    		std::cout << "		   Test quantity: " << test << std::endl;
   		    }
    		
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

			std::cout << "---------------- Numerical solution norm(s) --------------------" << std::endl;
			Real normVel = HM.normL2( solution_3D, "velocity" );
			Real normPress = HM.normL2( solution_3D, "pressure" );
			std::cout << "L2 velocity norm: " << normVel << std::endl;
    	    std::cout << "L2 pressure norm: " << normPress << std::endl;
    	    std::cout << "----------------------------------------------------------------" << std::endl;
    	    
		} std::cout<<__LINE__<<std::endl;// for loop on t
    }std::cout<<__LINE__<<std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
