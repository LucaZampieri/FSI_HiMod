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

    @author Matteo Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi89@gmail.com>
    @date 01-09-2013

    In this tutorial we explain how to solve an ADR problem
    with the Hierarchical Model Reduction method, with
    lateral Dirichlet boundary conditions.
 */

// ----------------------------------------------------//
// For now there is only a serial version of the code  //
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

// -----------------------------------------------------//
// These headers are needed for the solution of the     //
// linear system                                        //
// -----------------------------------------------------//
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/filter/GetPot.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"

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
#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>
#include <lifev/himod/modalbasis/HiModAssembler.hpp>

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

// Now we define a set of data to assemble a problem with known solution
// To compute the analytical expression of the force term
// we used matlab symbolic toolbox

// The dirichlet inflow data
Real g (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real x = 0.0;
    Real Lx = 0.2;
    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}

// The force term
Real f (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    using std::exp;
    using std::pow;

    return     (2e7 * y * pow (Lx - x, 2) * (Ly - y) + 2e7 * z * pow (Lx - x, 2) * (Lz - z) + 4e7 * y * y * z * pow (Lx - x, 4) * (Ly - y) + 4e7 * y * z * z * pow (Lx - x, 4) * (Lz - z) - 4e7 * y * y * pow (Lx - x, 4) * (Ly - y) * (Lz - z) - 4e7 * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z) - 2e7 * y * z * (Ly - y) * (Lz - z) - 4e7 * y * y * z * z * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2) - 4e7 * y * z * z * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z) - 4e7 * y * y * y * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z) - 4e7 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) - 4e7 * y * y * y * z * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2) ) * exp (2 * y * z * pow (Lx - x, 2) );

}

// The exact solution
Real ues (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);


}

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
        std::cout << "----------------------------------------------------------------------------------" << std::endl;
        std::cout << "-  Resolution of an ADR problem using an Hierarchical Model Reduction approach   -" << std::endl;
        std::cout << "-               PACS project, by Matteo Aletti & Andrea Bortolossi               -" << std::endl;
        std::cout << "----------------------------------------------------------------------------------" << std::endl;


        // we use getpot to obtain some data given by the user

        // ********** GetPot **********
        GetPot command_line ( argc, argv );
        GetPot dataFile ( "data" );
        UInt verbose = dataFile ("set/verbose", 0);

        // total number of modes used in the transversal fiber
        UInt mtot = dataFile ("himod/m", 20);


        //-----------------------------------------------------------------------
        //          SPACES AND PROBLEM
        //----------------------------------------------------------------------

        // Mesh definition, we used a structured mesh for the dominant, supporting fiber
        // the number of elements is user defined with getpot
        if (verbose)
        {
            std::cout << "Mesh definition" << std::endl;
        }
        const UInt Nelements (dataFile ("mesh/num_elements", 10) );
        boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);
        regularMesh1D ( *fullMeshPtr, 0, Nelements, false, dataFile ("mesh/lx", 1.), 0.0);


        if (verbose)
        {
            std::cout << "FESpace definition" << std::endl;
        }

        // We define here the FE space on the supporting fiber, only P1 should be used
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
        ( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, "P1", 1, Comm) );

        // display the number of FEM dof
        DOF DataFESpace (uSpace->dof() );
        UInt numdof = DataFESpace.numTotalDof();
        if (verbose)
        {
            std::cout << "Degrees of freedom FEM = " << numdof << std::endl;
        }

        // we start a chrono to understand the time needed for the assembling
        LifeChrono himodchrono;
        himodchrono.start();

        // Read the dimension of the transveral rectangle from getpot
        Real Ly = dataFile ("himod/ly", 2.);
        Real Lz = dataFile ("himod/lz", 2.);

        if (verbose)
        {
            std::cout << "ModalBasis definition" << std::endl;
        }

        // creation of the modal space
        boost::shared_ptr<ModalSpaceRectangular> MB ( new ModalSpaceRectangular (Ly, Lz, mtot) );

        if (verbose)
        {
            std::cout << "Adding BC on the slices" << std::endl;
        }

        // applying lateral homogeneous boundary conditions
        MB->addSliceBCY ("dir", "dir");
        MB->addSliceBCZ ("dir", "dir");

        // actually compute the basis
        MB->evaluateBasis();

        if (verbose)
        {
            // Some infos, from the ModalSpace
            MB->showMe();
            std::cout << "HiModAssembler definition" << std::endl;
        }
        // we mix the FEM space with the ModalSpace, creating an HiModSpace
        HiModAssembler< mesh_Type , matrix_Type , vector_Type, Rectangular > HM (uSpace, MB, Comm);


        if (verbose)
        {
            std::cout << "Creation of mapepetra for vector structured" << std::endl;
        }

        // we define a Map that will define the vector on the himod space
        // mtot blocks of monodimensional fem functions
        MapEpetra Map (mtot * DataFESpace.numTotalDof(), Comm);


        if (verbose)
        {
            std::cout << "Assembling systemMatrix with an ADR problem" << std::endl;
        }
        //Using the map it is possible to define the system matrix
        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( Map ) );

        //We set here the block structure of the matrix following the theory
        std::vector<UInt> block_row (MB->mtot(), DataFESpace.numTotalDof() );
        std::vector<UInt> block_col (MB->mtot(), DataFESpace.numTotalDof() );
        systemMatrix->setBlockStructure (block_row, block_col);
        *systemMatrix *= 0.0;

        //Here we define the coeffcients of the bilinear form
        Real mu = 1.0;
        TreDvector_type beta;
        beta[0] = 0.0;
        beta[1] = 0.0;
        beta[2] = 0.0;
        Real sigma = 0.0;

        //With this command the matrix is filled
        HM.addADRProblem (systemMatrix, mu, beta, sigma);

        himodchrono.stop();

        if (verbose)
        {
            std::cout << "Time of assemble and modal basis creation : " << himodchrono.diff() << std::endl;
        }


        if (verbose)
        {
            std::cout << "Assembling rhs" << std::endl;
        }

        //Notes that this is not the only way to assemble the rhs, you can use the functor, for more details look test_addrhs.
        LifeChrono rhschrono;
        rhschrono.start();
        // also the rhs should have a block structure
        boost::shared_ptr<vector_Type> rhs (new vector_Type ( Map, Repeated ) );
        *rhs *= 0.0;
        rhs->setBlockStructure (block_row);

        // here we create a vector to contain the force term projected on the modal space and then interpolated on the 1D fem structure
        boost::shared_ptr<vector_Type> f_interpolated (new vector_Type ( Map, Repeated ) );
        *f_interpolated *= 0.0;
        f_interpolated->setBlockStructure (block_row);
        HM.interpolate (f, f_interpolated);

        // now that the f is actually defined on the himodspace it is possible to add it to the rhs
        HM.addrhs (rhs, f_interpolated);

        rhschrono.stop();

        if (verbose)
        {
            std::cout << "Time of assemble rhs : " << rhschrono.diff() << std::endl;
        }

        if (verbose)
        {
            std::cout << "Adding Dirichlet boundary conditions" << std::endl;
        }
        // we add dirichlet bc on the inflow boundary, via an algebraic penalization technique
        HM.addDirichletBC_In (systemMatrix, rhs, g);

        if (verbose)
        {
            std::cout << "Closing rhs and the matrix" << std::endl;
        }
        rhs->globalAssemble();
        systemMatrix->globalAssemble();

        // Some export for post processing, might be interesting to explore the sparsity pattern of the matrix
        // systemMatrix->spy ("SystemMatrix");
        // rhs->spy ("rhs");

        //-----------------------------------------------------------------------
        //              SOLVE
        //-----------------------------------------------------------------------

        LifeChrono solvechrono;
        solvechrono.start();

        if (verbose)
        {
            std::cout << std::endl << "[Solvers initialization]" << std::endl;
        }
        boost::shared_ptr<prec_Type> precPtr (new prec_Type);
        precPtr->setDataFromGetPot (dataFile, "prec");

        if (verbose)
        {
            std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
        }
        Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

        LinearSolver linearSolver;
        linearSolver.setCommunicator ( Comm );
        linearSolver.setParameters ( *belosList );
        linearSolver.setPreconditioner ( precPtr );

        if (verbose)
        {
            std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
        }
        boost::shared_ptr<vector_Type> solution;
        solution.reset ( new vector_Type ( Map, Unique ) );
        linearSolver.setOperator ( systemMatrix );
        linearSolver.setRightHandSide ( rhs );
        linearSolver.solve ( solution );

        solvechrono.stop();

        if (verbose)
        {
            std::cout << "Time for solve the problem: " << solvechrono.diff() << std::endl;
        }

        if (verbose)
        {
            std::cout << "Starting post processing on the solution" << std::endl;
        }

        //here we start the postprocess of the solution
        UInt nquadY = MB->qrY().nbQuadPt();
        UInt nquadZ = MB->qrZ().nbQuadPt();
        // We want to evaluate the solution on a 3D grid made by the FEM grid times the quadrature grid
        MapEpetra Map_3D (nquadY * nquadZ * numdof, Comm);
        vector_Type solution_3D (Map_3D, Unique);
        // This method does exactly this.
        solution_3D = HM.evaluateBase3DGrid (*solution);

        vector_Type ues_3D (Map_3D, Unique);
        //This is the same but from an analytical function

        ues_3D = HM.evaluateBase3DGrid (ues);

        vector_Type err_3D (Map_3D, Unique);
        err_3D = solution_3D;
        err_3D += ues_3D * (-1);

        // now that approximate solution and exact solution are evaluate on the grid it is easy to compute the L2 norm
        Real normL2 = HM.normL2 (err_3D);
        Real norm = HM.normL2 (ues_3D);

        // now we show how it is possible to export, in VTK format, the himod solution on a self-defined structured grid
        UInt n = 30;
        HM.exporterStructuredVTK (n, n, n, solution, dataFile, "DDDD", "solution");
        HM.exporterFunctionVTK (n, n, n, ues, dataFile, "DDDD", "forceterm");

        std::cout << "---------------- Compute errors --------------------" << std::endl;
        std::cout << "L2 error norm: " << normL2 << std::endl;
        std::cout << "L2 relative error norm: " << normL2 / norm << std::endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
