#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
//#pragma GCC diagnostic warning "-Wunused-variable"
//#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings:
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

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/util/LifeChrono.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>

#include <lifev/himod/modalbasis/HiModAssembler.hpp>

#include <lifev/himod/util/GeneralTest.hpp>

#include <lifev/himod/util/CaseTest.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

using namespace LifeV;

typedef RegionMesh<LinearLine>              mesh_Type;
typedef MatrixEpetraStructured<Real>        matrix_Type;
typedef VectorEpetraStructured              vector_Type;
typedef VectorSmall<3>                      TreDvector_type;

typedef LifeV::Preconditioner                      basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
typedef PreconditionerIfpack                prec_Type;
typedef boost::shared_ptr<prec_Type>        precPtr_Type;

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
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "-                             A simple test                                      -" << std::endl;
    std::cout << "-               PACS project, by Matteo Aletti & Andrea Bortolossi               -" << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;


    // ********** GetPot **********
    GetPot command_line ( argc, argv );
    GetPot dataFile ( "data" );
    //****************************


    std::cout << "#----------------------------------------------#" << std::endl;
    std::cout << "           Starting show test  : "  << dataFile ("himod/casetest", "") << std::endl;
    std::cout << "#----------------------------------------------#" << std::endl;


	//The use of the GeneralTest class is very simple
	//You just write the following 3 lines.
	//What you really have to do is modify the data for getpot
	/*
	[set]
		verbose = 0
	[../]

	[himod]
		casetest=DD2D// there are several test you can try
		#	 ___3
		#  4	|   |
		#	|___| 2
		#	  1
		#
		# DDDD		
    	 # DDDD_Adv	
    	 # DDDD_ADR	
    	 # DRDR		
    	 # RRRR
		# RR2D
		# BDRR	
		# BDDD	
		# DD2D
		# Camini (DDDD with omogeneous inflow)
	
		m = 25  //number of modes
		Nel = 10 //number of elements
	
		[./export]
			nx	= 21  //Data for the grid
			ny	= 32
			nz	= 32
			solution = Soluzione  //part of the name of the output files
			forceterm = Forzante
		[../]
	[../]
	*/
    GeneralTest Test (&quadRuleLobSeg32pt, &quadRuleLobSeg32pt);

    Test.setDatafile (dataFile);

    Test.run();

	// If you want to see some results try to visualize the solution in paraview.


    std::cout << "#----------------------------------------------#" << std::endl;
    std::cout << "           End of show test   : " << dataFile ("himod/casetest", "") << std::endl;
    std::cout << "#----------------------------------------------#" << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
