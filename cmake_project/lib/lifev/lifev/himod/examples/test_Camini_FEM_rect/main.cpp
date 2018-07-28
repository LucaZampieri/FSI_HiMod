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
    @brief Test camini FEM

    @author Matteo Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi89@gmail.com>
    @date 01-09-2013

	We use this case to compare 3D picture with HiMod Camini casetest
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




#include <lifev/core/array/MatrixEpetra.hpp>

//#include <lifev/core/filter/ExporterEnsight.hpp>
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"

#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"
#pragma GCC diagnostic warning "-Wdeprecated-declarations"

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


Real epsilon (0.05);

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
        return 0;
    }

    return NAN;
}

Real fRhs ( const Real& /* t */, const Real& x, const Real& y, const Real&  z  , const ID& /* i */ )
{

    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    Real frac1X = 1. / 8.;
    Real frac1Y = 1. / 8.;
    Real frac1Z = 1. / 8.;

    Real frac2X = 2. / 8.;
    Real frac2Y = 7. / 8.;
    Real frac2Z = 7. / 8.;

    Real thick = 900;

    using std::exp;
    using std::pow;

    return  1e2 * ( exp ( -thick * ( pow (x - Lx * frac1X, 2) + pow (y - Ly * frac1Y, 2) + pow (z - Lz * frac1Z, 2) ) )
                    + exp ( -thick * ( pow (x - Lx * frac2X, 2) + pow (y - Ly * frac2Y, 2) + pow (z - Lz * frac2Z, 2) ) ) );

}

Real g_inflow ( const Real& /* t */, const Real& /*x*/, const Real& /* y */, const Real& /* z */, const ID& /* i */ )
{
    return  0;
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
    const UInt Nelements (dataFile ("mesh/nelements", 10) );
    const UInt NelementsSlice (dataFile ("mesh/nelementsSLICE", 10) );
    // Build and partition the mesh
    boost::shared_ptr< mesh_Type > fullMeshPtr ( new RegionMesh<LinearTetra> ( Comm ) );
    regularMesh3D ( *fullMeshPtr, 1, Nelements, NelementsSlice, NelementsSlice, false,
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

    //Chrono definition
    LifeChrono ADDproblem;
    ADDproblem.start();

    // Perform the assembly of the matrix
    adrAssembler.addDiffusion (systemMatrix, epsilon);
    vector_Type beta (betaFESpace->map(), Repeated);
    betaFESpace->interpolate (betaFct, beta, 0.0);
    adrAssembler.addAdvection (systemMatrix, beta);
    adrAssembler.addMass (systemMatrix, 0.3);
    systemMatrix->globalAssemble();



    // Definition and assembly of the RHS
    vector_Type rhs (uFESpace->map(), Repeated);
    rhs *= 0.0;
    vector_Type fInterpolated (uFESpace->map(), Repeated);
    fInterpolated *= 0.0;
    uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
    adrAssembler.addMassRhs (rhs, fInterpolated);
    rhs.globalAssemble();

    // Definition and application of the BCs
    BCHandler bchandler;
    BCFunctionBase BCu ( g_inflow );
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

    ADDproblem.stop();
    std::cout << "************************************" << std::endl;
    std::cout << "GDL 2D = " << NelementsSlice* NelementsSlice << std::endl;
    std::cout << "Tempo di assemblaggio = " << ADDproblem.diff() << std::endl;
    std::cout << "************************************" << std::endl;

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

    //Exporter
    ExporterVTK<mesh_Type> exporter (dataFile, "ADR3D");
    boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
    boost::shared_ptr<vector_Type> betaPtr (new vector_Type (beta, Repeated) );
    exporter.addVariable (  ExporterData<mesh_Type>::ScalarField, "solution", uFESpace,
                            solutionPtr, UInt (0), ExporterData<mesh_Type>::SteadyRegime, ExporterData<mesh_Type>::Node);
    exporter.addVariable (  ExporterData<mesh_Type>::VectorField, "trasport", betaFESpace,
                            betaPtr, UInt (0), ExporterData<mesh_Type>::SteadyRegime, ExporterData<mesh_Type>::Node);
    exporter.setMeshProcId ( uFESpace->mesh() , 0 );
    exporter.postProcess (0);


    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


