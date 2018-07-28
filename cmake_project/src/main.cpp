#include <boost/lexical_cast.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include "../include/ReferenceMap.hpp"
#include "../include/FSIData.hpp"
#include "../include/NSModalSpaceCircular.hpp"
#include "../include/NSHiModAssembler.hpp"
#include "../include/FSISolver.hpp"
#include "../include/QuadratureRule.hpp"
#include <boost/shared_ptr.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <vector>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <iostream>

using namespace LifeV;

typedef RegionMesh<LinearLine> mesh_Type;
typedef MatrixEpetraStructured<Real> matrix_Type;
typedef VectorEpetraStructured vector_Type;

typedef NSHiModAssembler<mesh_Type, matrix_Type, vector_Type, Circular> assembler_Type;

int main( int argc, char* argv[] )
{
  bool verbose = 1;

  if (verbose)
  {
    std::cout <<"Definition of the communicator: ";
  }
#ifdef HAVE_MPI
   MPI_Init (&argc, &argv);
#endif
#ifdef HAVE_MPI
   boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
   ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
   boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif
  if (verbose)
  {
    std::cout <<"DONE" <<std::endl;
  }

  GetPot commandLine( argc, argv );
  std::cout << "    GETTING THE DATA \n\n\n";
  std::string dataPath = "/home/zampieri/Desktop/FSI_HiMod/cmake_project/data/data";
  GetPot dataFile( dataPath );
  if (verbose)
  {
    std::cout <<"Definition of FSIData: ";
  }
  FSIData data( dataFile );
  if (verbose)
  {
    std::cout <<"DONE" <<std::endl;
  }
  data.printAll();

  if (verbose)
  {
    std::cout <<"Definition of NSModalSpaceCircular: ";
  }
  boost::shared_ptr< mesh_Type > fullMeshPtr( new mesh_Type );
  regularMesh1D( *fullMeshPtr, 0, data.Nelements(), false, data.L(), 0.0 );

  boost::shared_ptr<FESpace<mesh_Type, MapEpetra>> uSpace( new FESpace<mesh_Type, MapEpetra>( fullMeshPtr, "P2", 1, Comm ) );
  boost::shared_ptr<FESpace<mesh_Type, MapEpetra>> pSpace( new FESpace<mesh_Type, MapEpetra>( fullMeshPtr, "P1", 1, Comm ) );

  std::vector<Real> nodes( data.Nelements()*2+1 );
  for ( UInt i = 0; i < (data.Nelements()+1); ++i )
  {
    nodes[2*i] = (uSpace->mesh()->pointList)(i).x();
    if (i != data.Nelements())
      nodes[2*i+1] = ( (uSpace->mesh()->pointList)(i+1).x() + (uSpace->mesh()->pointList)(i).x() ) / 2;
  }
  ReferenceMap refMap( nodes, data.R(), &quadRuleSeg32pt, &quadRuleFQSeg32pt );

  boost::shared_ptr<NSModalSpaceCircular> MB( new NSModalSpaceCircular( data.R(), data.theta(), data.mx(), data.mr(), data.mtheta(), data.mp(), nodes.size(), &refMap, &quadRuleBoundary, &(refMap.qrRho()), &(refMap.qrTheta()) ) );
  MB->evaluateBasisFSI( "Zernike", "ZernikeNatural", "Zernike", "ZernikeNatural", 0 );
  if (verbose)
  {
    std::cout <<"DONE" <<std::endl;
  }

  if (verbose)
  {
    std::cout <<"Definition of NSHiModAssemblerCircular: ";
  }
  boost::shared_ptr<assembler_Type> HM( new assembler_Type( uSpace, pSpace, MB, Comm ) );
  {
    std::cout <<"DONE" <<std::endl;
  }

  if (verbose)
  {
    std::cout <<"Definition of FSISolver: ";
  }
  FSISolver solver( HM, dataFile, &refMap );
  {
    std::cout <<"DONE" <<std::endl;
  }

  solver.solve();

  std::cout <<"THE END !!!" <<std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
