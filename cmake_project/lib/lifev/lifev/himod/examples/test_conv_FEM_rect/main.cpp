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
    @file
    @brief

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"

#include <lifev/core/array/MatrixEpetra.hpp>

//#include <lifev/core/filter/ExporterEnsight.hpp>

#pragma GCC diagnostic ignored "-Wconversion"
#include <lifev/core/filter/ExporterVTK.hpp>
#pragma GCC diagnostic warning "-Wconversion"

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/mesh/MeshData.hpp>

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


Real epsilon (1.0);

Real betaFct ( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    if (i == 0)
    {
        return 0;
    }
    if (i == 1)
    {
        return 0;
    }
    if (i == 2)
    {
        return 0;
    }
}

Real g_adr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real x = 0.0;
    Real Lx = 0.2;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}



Real fRhs (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    using std::exp;
    using std::pow;


    return (
               2e7 * y * pow (Lx - x, 2) * (Ly - y)
               + 2e7 * z * pow (Lx - x, 2) * (Lz - z)
               + 4e7 * y * y * z * pow (Lx - x, 4) * (Ly - y)
               + 4e7 * y * z * z * pow (Lx - x, 4) * (Lz - z)
               - 4e7 * y * y * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4e7 * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 2e7 * y * z * (Ly - y) * (Lz - z)
               - 1e7 * y * z * pow (Lx - x, 2) * (Ly - y)
               - 1e7 * y * z * pow (Lx - x, 2) * (Lz - z)
               + 1e7 * y * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               + 1e7 * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 5e7 * y * z * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
               + 3e7 * y * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * z * z * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               + 2e7 * y * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               + 2e7 * y * y * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4e7 * y * z * z * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * y * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * y * z * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               - 1e8 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
           ) * exp (2 * y * z * pow (Lx - x, 2) );
}

Real exactSolution (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}



typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);

    // Read first the data needed


    GetPot dataFile ( "data" );
    const UInt NelementsX (dataFile ("mesh/nelements", 10) );
    UInt NelementsY (2);
    UInt NelementsZ (2);

    const UInt maxit = 7;


    //-----------------------------------------------------------------------------------------------------
    boost::shared_ptr< mesh_Type > fullMeshPtraus ( new RegionMesh<LinearTetra> ( Comm ) );
    regularMesh3D ( *fullMeshPtraus, 1, 40, 20, 20, false,
                    0.2,   0.1,   0.1,
                    0.0,  0.0,  0.0);

    boost::shared_ptr< mesh_Type > meshPtraus;
    {
        MeshPartitioner< mesh_Type >   meshPartaus (fullMeshPtraus, Comm);
        meshPtraus = meshPartaus.meshPartition();
    }

    fullMeshPtraus.reset();

    std::string uOrd ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpaceaus ( new FESpace< mesh_Type, MapEpetra > (meshPtraus, uOrd, 1, Comm) );

    vector_Type vect (uFESpaceaus->map(), Unique);
    vect *= 0.0;
    uFESpaceaus->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), vect, 0.0 );

    Real l2uex (uFESpaceaus->l2Norm (vect) );
    //-----------------------------------------------------------------------------------------------------

    for (UInt i (0); i < maxit; ++i)
    {
        // Build and partition the mesh

        boost::shared_ptr< mesh_Type > fullMeshPtr ( new RegionMesh<LinearTetra> ( Comm ) );
        regularMesh3D ( *fullMeshPtr, 1, NelementsX, NelementsY, NelementsZ, false,
                        0.2,   0.1,   0.1,
                        0.0,  0.0,  0.0);

        boost::shared_ptr< mesh_Type > meshPtr;
        {
            MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
            meshPtr = meshPart.meshPartition();
        }

        fullMeshPtr.reset();

        // Build the FESpaces

        std::string uOrder ("P1");
        std::string bOrder ("P1");
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, Comm) );
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > betaFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, bOrder, 3, Comm) );

        // Build the assembler and the matrices

        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;

        adrAssembler.setup (uFESpace, betaFESpace);

        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );
        *systemMatrix *= 0.0;

        // Perform the assembly of the matrix

        adrAssembler.addDiffusion (systemMatrix, epsilon);

        // vector_Type beta (betaFESpace->map(), Repeated);
        // betaFESpace->interpolate (betaFct, beta, 0.0);
        // adrAssembler.addAdvection (systemMatrix, beta);

        //adrAssembler.addMass (systemMatrix, 0.3);

        systemMatrix->globalAssemble();

        // Definition and assembly of the RHS

        //vector_Type rhs(uFESpace->map(),Unique);
        vector_Type rhs (uFESpace->map(), Repeated);
        rhs *= 0.0;

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        fInterpolated *= 0.0;
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();

        // Definition and application of the BCs

        BCHandler bchandler;
        BCFunctionBase BCu ( exactSolution );

        bchandler.addBC ("Dirichlet", 1, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", 3, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", 4, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", 5, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", 6, Essential, Full, BCu, 1);

        bchandler.addBC ("Neumann", 2, Natural, Full, BCu, 1);

        bchandler.bcUpdate (*uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

        vector_Type rhsBC (rhs, Unique);
        bcManage (*systemMatrix, rhsBC, *uFESpace->mesh(), uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = rhsBC;

        // Definition of the solver

        SolverAztecOO linearSolver;

        linearSolver.setDataFromGetPot (dataFile, "solver");
        linearSolver.setupPreconditioner (dataFile, "prec");

        linearSolver.setMatrix (*systemMatrix);

        linearSolver.setCommunicator (Comm);

        // Definition of the solution

        vector_Type solution (uFESpace->map(), Unique);
        solution *= 0.0;

        // Solve the solution

        linearSolver.solveSystem (rhsBC, solution, systemMatrix);

        // Error computation

        vector_Type solutionErr (solution);
        solutionErr *= 0.0;
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), solutionErr, 0.0 );
        solutionErr -= solution;
        solutionErr.abs();

        Real l2error (uFESpace->l2Error (exactSolution, vector_Type (solution, Repeated), 0.0) );



        std::fstream fileOut;
        fileOut.open ("ConvFEM.txt", std::fstream::out | std::fstream::app);

        if (fileOut.is_open() == false)
        {
            std::cerr << "File not opened" << std::endl;
            exit (1);
        }


        fileOut << NelementsY* NelementsZ << '\t' << l2error / l2uex << std::endl;

        fileOut.close();

        NelementsY += i;
        NelementsZ += i;
    }
    /*
        if (l2error > 0.0055)
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }
        if (linferror > 0.0046)
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }
    */
    // Exporter definition and use

    /*


        if (verbose)
        {
            std::cout << " -- Defining the exporter ... " << std::flush;
        }
        ExporterVTK<mesh_Type> exporter(dataFile,"ADR3D");
        //ExporterEnsight<mesh_Type> exporter ( dataFile, meshPtr, "solution", Comm->MyPID() ) ;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Defining the exported quantities ... " << std::flush;
        }
        boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
        boost::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );
        boost::shared_ptr<vector_Type> betaPtr (new vector_Type (beta, Repeated) );

    if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Updating the exporter ... " << std::flush;
        }
        exporter.addVariable (  ExporterData<mesh_Type>::ScalarField, "solution", uFESpace,
                            solutionPtr, UInt(0),ExporterData<mesh_Type>::SteadyRegime, ExporterData<mesh_Type>::Node);


        exporter.addVariable (  ExporterData<mesh_Type>::VectorField, "trasport", betaFESpace,
                            betaPtr, UInt(0),ExporterData<mesh_Type>::SteadyRegime, ExporterData<mesh_Type>::Node);

     exporter.setMeshProcId ( uFESpace->mesh() , 0 );


    //    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "solution", uFESpace, solutionPtr, UInt (0) );
      //  exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "error", uFESpace, solutionErrPtr, UInt (0) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Exporting ... " << std::flush;
        }
        exporter.postProcess (0);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }*/

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


