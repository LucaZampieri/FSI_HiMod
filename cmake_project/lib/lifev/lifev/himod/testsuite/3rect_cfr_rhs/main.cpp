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

using namespace LifeV;

typedef RegionMesh<LinearLine>              mesh_Type;
typedef MatrixEpetraStructured<Real>        matrix_Type;
typedef VectorEpetraStructured              vector_Type;
typedef VectorSmall<3>                      TreDvector_type;

typedef LifeV::Preconditioner                      basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
typedef PreconditionerIfpack                prec_Type;
typedef boost::shared_ptr<prec_Type>        precPtr_Type;

Real f (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return     2 * y * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * (Ly - y)
               + 2 * z * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * (Lz - z)
               + 4 * z * std::pow (y, 2) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow ( (Lx - x), 2) * (Ly - y)
               + 4 * y * std::pow (z, 2) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow ( (Lx - x), 2) * (Lz - z)
               - 4 * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow (y, 2) * std::pow (z, 2) * (Ly - y) * (Lz - z)
               - 4 * std::pow (y, 2) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4 * std::pow (z, 2) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4 * y * std::pow (z, 3) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4 * z * std::pow (y, 3) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * std::pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4 * std::pow (y, 3) * std::pow (z, 3) * std::exp (2 * y * z * std::pow ( (Lx - x), 2) ) * (Ly - y) * (Lz - z) * std::pow (2 * Lx - 2 * x, 2);

}

//In this example we show the possibility to use a a different system to assemble the right hand side which is slower because it
// evaluates the function over the quadrature nodes, while the faster evaluate the function over the 1D mesh and it interpolate it linearly.
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
    std::cout << "-                        Testing Addrsh vs Addrhs_highprec                       -" << std::endl;
    std::cout << "-               PACS project, by Matteo Aletti & Andrea Bortolossi               -" << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;

    // ********** GetPot **********
    GetPot command_line ( argc, argv );
    // const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( "data" );
    //****************************


    // number of modes in y e z direction
    UInt mtot = dataFile ("himod/m", 20);

    //-----------------------------------------------------------------------
    //          SPACES AND PROBLEM
    //----------------------------------------------------------------------

    //Mesh definition

    const UInt Nelements (dataFile ("mesh/num_elements", 10) );

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh1D ( *fullMeshPtr, 0, Nelements, false, dataFile ("mesh/lx", 1.), 0.0);

    //FESpace and ETFESpace
    std::string Poly_type = dataFile ("mesh/Poly_type", "P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, Poly_type, 1, Comm) );

    //CHECK DOF
    DOF DataFESpace (uSpace->dof() );
    UInt numdof = DataFESpace.numTotalDof();
    std::cout << "FEM dof = " << numdof << std::endl;

    //MODAL BASIS CLASS
    Real Ly = dataFile ("himod/ly", 2.);
    Real Lz = dataFile ("himod/lz", 2.);

    boost::shared_ptr<ModalSpaceRectangular> MB ( new ModalSpaceRectangular (Ly, Lz, mtot) );
    std::cout << "Dimension of the modal basis = " << mtot << std::endl;

    MB->addSliceBCY ("dir", "dir");
    MB->addSliceBCZ ("dir", "dir");

    MB->evaluateBasis();

    //HIMOD CLASS
    HiModAssembler< mesh_Type , matrix_Type , vector_Type, Rectangular > HM (uSpace, MB, Comm);
    std::vector<UInt> block_row (MB->mtot(), DataFESpace.numTotalDof() );

    MapEpetra Map (mtot * DataFESpace.numTotalDof(), Comm);

    // Using the functor
    LifeChrono functorcrono;
    functorcrono.start();

    boost::shared_ptr<vector_Type> rhs_functor (new vector_Type ( Map, Repeated ) );
    *rhs_functor *= 0.0;
    rhs_functor->setBlockStructure (block_row);

    HM.addrhs_HiPrec (rhs_functor, f);

    functorcrono.stop();

    // Using the standard Addrhs
    LifeChrono standardcrono;
    standardcrono.start();

    boost::shared_ptr<vector_Type> f_interpolated (new vector_Type ( Map, Repeated ) );
    *f_interpolated *= 0.0;
    f_interpolated->setBlockStructure (block_row);

    boost::shared_ptr<vector_Type> rhs_standard (new vector_Type ( Map, Repeated ) );
    *rhs_standard *= 0.0;
    rhs_standard->setBlockStructure (block_row);

    HM.interpolate (f, f_interpolated);
    HM.addrhs (rhs_standard, f_interpolated);

    standardcrono.stop();

    // Compute the error
    *rhs_standard += *rhs_functor * (-1);
    std::cout << "Difference between rhs_standard and rhs_functor = " << rhs_standard->norm2() / rhs_functor->norm2() << std::endl;
    std::cout << "Time to assemble rhs_standard = " << standardcrono.diff() << std::endl;
    std::cout << "Time to assemble rhs_functor = " << functorcrono.diff() << std::endl;

    /*!
        COMMENT:
        Notes that the functor uses exactly double time in respect to the interpolate approach, this because the Addrhs_HiPrec evaluate the function over all quadrature nodes, which are more than the actual nodes of the 1D mesh.
    */
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
