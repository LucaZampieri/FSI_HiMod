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
    @brief Test convergence of FEM

    @author Matteo Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi89@gmail.com>
    @date 01-09-2013

	We use this case to test the convergence of FEM by increment
	DOF only on the slice of the domain. This is usefull to make some
	intresting matches with HiMod. For the moment the comparison is done
	again case test DDDD_Adr.
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

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"
#pragma GCC diagnostic warning "-Wdeprecated-declarations"

#include <lifev/core/array/MatrixEpetra.hpp>

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
        return 5;
    }
    if (i == 1)
    {
        return 1;
    }
    if (i == 2)
    {
        return 1;
    }

    return NAN;
}

Real g_adr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real x = 0.0;
    Real Lx = 0.2;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}



Real f_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
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

Real ues_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
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

    const UInt maxit = 6;


    //-----------------------------------------------------------------------------------------------------
    boost::shared_ptr< mesh_Type > fullMeshPtraus ( new RegionMesh<LinearTetra> ( Comm ) );
    regularMesh3D ( *fullMeshPtraus, 1, 30, 20, 20, false,
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
    uFESpaceaus->interpolate ( static_cast<feSpace_Type::function_Type> ( ues_adr ), vect, 0.0 );

    Real l2uex (uFESpaceaus->l2Norm (vect) );
    //-----------------------------------------------------------------------------------------------------

       std::fstream fileOut;
        fileOut.open ("ConvFEM.dat", std::fstream::out);

        if (fileOut.is_open() == false)
        {
            std::cerr << "File not opened" << std::endl;
            exit (1);
        }

		LifeChrono FEMchrono;

    for (UInt i (1); i < maxit; ++i)
    {

        FEMchrono.start();

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

        adrAssembler.addMass (systemMatrix, 3.0);

        systemMatrix->globalAssemble();

        // Definition and assembly of the RHS

        //vector_Type rhs(uFESpace->map(),Unique);
        vector_Type rhs (uFESpace->map(), Repeated);
        rhs *= 0.0;

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        fInterpolated *= 0.0;
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( f_adr ), fInterpolated, 0.0 );
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();

        // Definition and application of the BCs

        BCHandler bchandler;
        BCFunctionBase BCu ( g_adr );

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

		FEMchrono.stop();

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
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( ues_adr ), solutionErr, 0.0 );
        solutionErr -= solution;
        solutionErr.abs();

        Real l2error (uFESpace->l2Error (ues_adr, vector_Type (solution, Repeated), 0.0) );

        fileOut << NelementsY* NelementsZ << '\t' << l2error / l2uex <<'\t'<<FEMchrono.diff()<< std::endl;

        NelementsY += i;
        NelementsZ += i;
    }

fileOut.close();

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


