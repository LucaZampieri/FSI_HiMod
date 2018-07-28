#include <boost/lexical_cast.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include "include/ReferenceMap.hpp"
#include "include/FSIData.hpp"
#include "include/NSModalSpaceCircular.hpp"
#include "include/NSHiModAssembler.hpp"
#include "include/FSISolver.hpp"
#include "include/QuadratureRule.hpp"
#include <boost/shared_ptr.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <vector>
#include <lifev/core/fem/QuadratureRule.hpp>

#include "include/poiseuille.hpp" // animalata ma rgossa

#include <iostream>
#include <stdexcept>

#define BASEDIR "./"
#define VERBOSE 1

using namespace LifeV;

typedef RegionMesh<LinearLine> mesh_Type;
typedef MatrixEpetraStructured<Real> matrix_Type;
typedef VectorEpetraStructured vector_Type;

typedef NSHiModAssembler<mesh_Type, matrix_Type, vector_Type, Circular> assembler_Type;



inline void file_exists (const std::string& name) {
  // Function that test if a file is findable, throws an exception otherwise
    ifstream f(name.c_str());
    if (!f.good()){
      throw "The file cannot be found! \nInterupting the program";
    }
}


int main( int argc, char* argv[] )
{
  // ------------- Set some variables
  bool verbose = VERBOSE;
  std::string baseDir = BASEDIR; // set the base directory


  // ------- Set comunicators  -----------------------
  if (verbose) { std::cout <<"Definition of the communicator: ";}
#ifdef HAVE_MPI
   MPI_Init (&argc, &argv);
#endif
#ifdef HAVE_MPI
   boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
   ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
   boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif
  if (verbose) {std::cout <<"DONE" <<std::endl;}

  GetPot commandLine( argc, argv );

  // ----------  Get Data  ----------
  std::cout << "    GETTING THE DATA \n\n\n";
  std::string dataFileName = baseDir+"data/data";
  std::cout <<"Searching for data in :"
            <<dataFileName<<std::endl;
  // Test if the dataFile exists
  try{
    file_exists(dataFileName);
  }  catch (const char* msg) {std::cerr << msg << endl;return 0;}
  // Get data
  GetPot dataFile( dataFileName );
  if (verbose) { std::cout <<"Definition of FSIData: "; }
  FSIData data( dataFile );
  data.printAll();

  // ---------- Define Fespaces/ModalSpaceCircular  ----------
  std::cout <<"Definition of NSModalSpaceCircular: ";
  boost::shared_ptr< mesh_Type > fullMeshPtr( new mesh_Type );
  regularMesh1D( *fullMeshPtr, 0, data.Nelements(), false, data.L(), 0.0 );
  // FE spaces
  boost::shared_ptr<FESpace<mesh_Type, MapEpetra>> uSpace( new FESpace<mesh_Type, MapEpetra>( fullMeshPtr, "P2", 1, Comm ) );
  boost::shared_ptr<FESpace<mesh_Type, MapEpetra>> pSpace( new FESpace<mesh_Type, MapEpetra>( fullMeshPtr, "P1", 1, Comm ) );

  std::vector<Real> uNodes( data.Nelements()*2+1 ); // 2n+1 nodes if it is P2
  for ( UInt i = 0; i < (data.Nelements()+1); ++i ) // n+1 nodes for P1
  {
    uNodes[2*i] = (uSpace->mesh()->pointList)(i).x();
    if (i != data.Nelements())
      uNodes[2*i+1] = ( (uSpace->mesh()->pointList)(i+1).x() + (uSpace->mesh()->pointList)(i).x() ) / 2;
  }
  // Define the map to the reference domain
  /*
  onst function_Type& jr,
                const function_Type& jtheta,
                const function_Type& dr,
                const function_Type& dtheta,
                const function_Type& drtheta,
                const function_Type& dthetar,
                const function_Type& jacobian,
                const function_Type& jacobianWall, // Luca
                const function_Type& inverseRhat,
                const QuadratureRule* quadrho = &quadRuleSeg32pt,
                const QuadratureRule* quadtheta = &quadRuleSeg32pt )
*/

  ReferenceMap    refMap( Jr, Jtheta, Dr, Dtheta, Drtheta, Dthetar, Jacobian, JacobianWall, inverseRhat, uSpace->map(),
                              &quadRuleSeg32pt, &quadRuleSeg64pt );

  //ReferenceMap refMap( uNodes, data.R(), &quadRuleSeg32pt, &quadRuleFQSeg32pt );
  // Define the Modal Basis
  boost::shared_ptr<NSModalSpaceCircular> MB( new NSModalSpaceCircular(
      data.Radius(), data.dRadius(), data.theta(), data.mx(), data.mr(), data.mtheta(), data.mp(),
      uNodes.size(), &refMap, &quadRuleBoundary, &(refMap.qrRho()), &(refMap.qrTheta()) ) );


  MB->evaluateBasisFSI( "Zernike", "ZernikeNatural", "Zernike", "ZernikeNatural", 0 );

  // Definition of NSHiModAssemblerCircular
  boost::shared_ptr<assembler_Type> HM( new assembler_Type( uSpace, pSpace, MB, Comm ) );


  // ---------- Create The solver and solves the problem ----------
  FSISolver solver( HM, dataFile, &refMap );
  solver.solve();

  std::cout <<"---------- THE END !!! ----------" <<std::endl;
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}


// End of file
