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

/* ========================================================

Simple ETA test to compute a block matrix associated to the Stokes problem

*/


/**
   @file main.hpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
*/


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>

#include "ETA_Blocks2DTest.hpp"


// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;

// ===================================================
//! Main
// ===================================================
int main (int argc, char* argv[])
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif
    bool verbose (Comm->MyPID() == 0);

    // Known errors
    const LifeV::Real knownMatrixError ( 4.55 );
    const LifeV::Real knownRhsError ( 0.0192705923625467 );


    // Tolerance
    const LifeV::Real tolerance ( 1e-10 );

    ETA_Blocks2DTest eta_blocks2DTest;

    // Error of the problem
    const std::vector<Real> errors = eta_blocks2DTest.run();
    const bool unsuccess = ( std::abs (errors[0] - knownMatrixError) > tolerance ||
                             std::abs (errors[1] - knownRhsError) > tolerance);

    if (unsuccess)
    {
        if (verbose)
        {
            std::cout << "End Result: TEST NOT PASSED" << std::endl;
        }
    }
    else
    {
        if (verbose)
        {
            std::cout << "End Result: TEST PASSED" << std::endl;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    if ( unsuccess )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }

}
