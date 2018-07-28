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
    @brief Tutorial showing the possibilities to use non constant coefficients

    @author Matteo Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi89@gmail.com>
    @date 01-09-2013

    In this tutorial we show how to use non constant coefficients
 */

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

#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>
#include <lifev/himod/modalbasis/HiModAssembler.hpp>

using namespace LifeV;

typedef RegionMesh<LinearLine>              mesh_Type;
typedef MatrixEpetraStructured<Real>        matrix_Type;
typedef VectorEpetraStructured              vector_Type;
typedef VectorSmall<3>                      TreDvector_type;
typedef LifeV::Preconditioner                      basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
typedef PreconditionerIfpack                prec_Type;
typedef boost::shared_ptr<prec_Type>        precPtr_Type;

// This file is the basic example to start the testing phase of
// the AddADRProblem with non constant coefficient.

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

// Now we define the coefficients as functions, now they are constants

// Diffusion
Real mu (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1.;
}
// Advection
Real beta (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 0:
            return 0.0;
        case 1:
            return 0.0;
        case 2:
            return 0.0;
    }
	return 0.0;
}
// Reaction
Real sigma (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

// the following code is the same as the first tutorial
int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif
    {

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

        GetPot command_line ( argc, argv );
        GetPot dataFile ( "data" );
        UInt verbose = dataFile ("set/verbose", 0);
        UInt mtot = dataFile ("himod/m", 20);
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
        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
        ( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, "P1", 1, Comm) );
        DOF DataFESpace (uSpace->dof() );
        UInt numdof = DataFESpace.numTotalDof();
        if (verbose)
        {
            std::cout << "Degrees of freedom FEM = " << numdof << std::endl;
        }
        LifeChrono himodchrono;
        himodchrono.start();
        Real Ly = dataFile ("himod/ly", 2.);
        Real Lz = dataFile ("himod/lz", 2.);
        if (verbose)
        {
            std::cout << "ModalBasis definition" << std::endl;
        }
        boost::shared_ptr<ModalSpaceRectangular> MB ( new ModalSpaceRectangular (Ly, Lz, mtot) );
        if (verbose)
        {
            std::cout << "Adding BC on the slices" << std::endl;
        }
        MB->addSliceBCY ("dir", "dir");
        MB->addSliceBCZ ("dir", "dir");
        MB->evaluateBasis();
        if (verbose)
        {
            std::cout << "HiModAssembler definition" << std::endl;
        }
        HiModAssembler< mesh_Type , matrix_Type , vector_Type, Rectangular > HM (uSpace, MB, Comm);
        if (verbose)
        {
            std::cout << "Creation of mapepetra for vector structured" << std::endl;
        }
        MapEpetra Map (mtot * DataFESpace.numTotalDof(), Comm);
        if (verbose)
        {
            std::cout << "Assembling systemMatrix with an ADR problem" << std::endl;
        }
        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( Map ) );
        std::vector<UInt> block_row (MB->mtot(), DataFESpace.numTotalDof() );
        std::vector<UInt> block_col (MB->mtot(), DataFESpace.numTotalDof() );
        systemMatrix->setBlockStructure (block_row, block_col);
        *systemMatrix *= 0.0;

		//There are no changes at all, not even now.
		//We should point out that, here, the code will be very slow, because we didn't make 
		//Some basic optization, like evaluating the coefficients on the quad grid before the loop.
        HM.addADRProblem (systemMatrix, mu, beta, sigma);
        himodchrono.stop();
        if (verbose)
        {
            std::cout << "Time of assemble: " << himodchrono.diff() << std::endl;
        }
        if (verbose)
        {
            std::cout << "Assembling rhs" << std::endl;
        }
        LifeChrono rhschrono;
        rhschrono.start();
        boost::shared_ptr<vector_Type> rhs (new vector_Type ( Map, Repeated ) );
        *rhs *= 0.0;
        rhs->setBlockStructure (block_row);
        boost::shared_ptr<vector_Type> f_interpolated (new vector_Type ( Map, Repeated ) );
        *f_interpolated *= 0.0;
        f_interpolated->setBlockStructure (block_row);
        HM.interpolate (f, f_interpolated);
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
        HM.addDirichletBC_In (systemMatrix, rhs, g);
        if (verbose)
        {
            std::cout << "Closing rhs and the matrix" << std::endl;
        }
        rhs->globalAssemble();
        systemMatrix->globalAssemble();
        // Some export for post processing, might be interesting to explore the sparsity pattern of the matrix
        systemMatrix->spy ("SystemMatrix");
        rhs->spy ("rhs");

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

        UInt nquadY = MB->qrY().nbQuadPt();
        UInt nquadZ = MB->qrZ().nbQuadPt();
        MapEpetra Map_3D (nquadY * nquadZ * numdof, Comm);
        vector_Type solution_3D (Map_3D, Unique);
        solution_3D = HM.evaluateBase3DGrid (*solution);
        vector_Type ues_3D (Map_3D, Unique);
        ues_3D = HM.evaluateBase3DGrid (ues);
        vector_Type err_3D (Map_3D, Unique);
        err_3D = solution_3D;
        err_3D += ues_3D * (-1);
        Real normL2 = HM.normL2 (err_3D);
        Real norm = HM.normL2 (ues_3D);
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
