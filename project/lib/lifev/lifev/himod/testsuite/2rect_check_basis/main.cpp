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
    @brief Tutorial

    @author Matteo Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi89@gmail.com>
    @date 01-09-2013
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

#include <boost/shared_ptr.hpp>

#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>
#include <lifev/himod/util/CheckModalBasis1D.hpp>

using namespace LifeV;
typedef ModalSpaceRectangular              modalbasis_type;
typedef boost::shared_ptr<modalbasis_type> modalbasis_ptrType;

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

#ifdef HAVE_MPI
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
    ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif


    Int mtot = 15;
    Real Ly = 1;
    Real Lz = 1;

	//Change this loop to test different number of modes or different length of the domain
	//Look at the explanation that will be displayed
    for (UInt test (0); test < 1; ++test)
    {
        mtot *= 2;
        std::cout << "**********************************" << std::endl;
        modalbasis_ptrType MB (new modalbasis_type (Ly, Lz, mtot) );
        MB->addSliceBCY ("dir", "dir");
        MB->addSliceBCZ ("dir", "dir");
        MB->evaluateBasis();
        CheckModalBasis1D checker (MB);
        checker.VerifyBC (1.0, 0.0);
        checker.VerifyOrthonormality();
        std::cout << std::endl;
    }
	std::cout<<" Try to modify the main.cpp to perform different tests with different BC"<<std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
