#ifndef FSI_SOLVER_HPP
#define FSI_SOLVER_HPP

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

#include "QuadratureRule.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include "FSIData.hpp"

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include "NSModalSpaceCircular.hpp"
#include "NSHiModAssembler.hpp"
#include "ReferenceMap.hpp"


#include <lifev/himod/tools/HiModExporterVtk.hpp>

#include <fstream>
#include <sstream>

#include <vector>
#include <functional>

namespace LifeV
{

/* Class for the solution of a fluid structure interaction problem in a cylindrical domain */

class FSISolver
{
public:

  // Define some aliases
  typedef RegionMesh<LinearLine>              mesh_Type;
  typedef MatrixEpetraStructured<Real>        matrix_Type;
  typedef VectorEpetraStructured              vector_Type;
  typedef boost::shared_ptr<NSHiModAssembler<mesh_Type, matrix_Type, vector_Type, Circular>> assembler_ptrType;

  typedef PreconditionerIfpack                prec_Type;
  typedef boost::shared_ptr<prec_Type>        precPtr_Type;

  typedef std::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& )> function_Type;
  typedef std::function<Real ( const Real& )> timeFunction_Type;

  typedef std::vector<std::vector<std::vector<std::tuple<Real, Real, Real>>>> grid_type;

  typedef std::vector<std::vector<std::vector<Real>>> Vector3D_type;


  // Constructor
  FSISolver( const assembler_ptrType& HM_,
             const GetPot& dataFile,
                   ReferenceMap* refMap_ );

  // Function to solve the problem in general
  void solve();

private:

  // Function to compute the solution of the linear system
  void solveSystem(const Real& t);

  // Function to compute the velocity of the ALE map using the displacement of the grid
  void computeALEVelocity();

  // Function to compute the displacement of the border of the domain using the velocity at the border and assuming independence of the angle
  void computeDisplacement();

  // Function to update the reference map
  void updateMap( const Real& t );

  // Function to update the computational grid
  void updateGrid();

  // Functions to convert the indices corresponding to the three coordinates to a general index for the VectorEpetraStructured
  UInt xcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
  {
    return ( j + k*nQuadRho + i*nQuadRho*nQuadTheta );
  }
  UInt rcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
  {
    return ( udof*nQuadRho*nQuadTheta + j + k*nQuadRho + i*nQuadRho*nQuadTheta );
  }

  UInt thetacoord2index( const UInt& i, const UInt& j, const UInt& k ) const
  {
    return ( 2*udof*nQuadRho*nQuadTheta + j + k*nQuadRho + i*nQuadRho*nQuadTheta );
  }
  UInt pcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
  {
    return ( 3*udof*nQuadRho*nQuadTheta + j + k*nQuadRho + i*nQuadRho*nQuadTheta );
  }
  /*UInt coord2index( const UInt& i, const UInt& j, const UInt& k )
  {
    return ( j + k*nQuadRho + i*nQuadRho*nQuadTheta );
  }
  */
  UInt coord2indexWall( const UInt& i, const UInt& k )
  {
    return ( k + i*nQuadTheta );
  }

  // Function to set up the grid during the construction of the solver
  void setGrid();

  void setCartesianGrid();

  // Function to set up the 3D containers of the solution during the construction of the solver
  void setContainers3D();

  // Function to convert the VectorEpetraStructured solution in a 3D solution
  void expandSolution();

  //void save( const Real& t, const UInt& iter ) const;
  void save( const UInt& iter ) const;

  // Assembler for the fluid structure interaction problem in a circular domain
  assembler_ptrType HM;

  DOF DataVelFESpace;
  DOF DataPressFESpace;

  // Class which collect all the data useful to solve the problem
  boost::shared_ptr<FSIData> data;

  // Reference map
  ReferenceMap* refMap;

  // computational grid
  grid_type grid;
  //grid_type gridOld;
  grid_type cartesianGrid;

  // Number of quadrature points
  Real uMeshSize;
  UInt udof;
  UInt pdof;
  Real nQuadRho;
  Real nQuadTheta;

  // Number of time steps
  UInt numbStep;

  // Communicator
  boost::shared_ptr<Epetra_Comm> Comm;

  // All the MapEpetra used
  MapEpetra Map;
  MapEpetra Map_3D;
  //MapEpetra Map_ALE;
  MapEpetra Map_Wall;
  MapEpetra Map_Eta;

  std::vector<UInt> block_row;
  std::vector<UInt> block_col;

  boost::shared_ptr<prec_Type> precPtr;
  Teuchos::RCP< Teuchos::ParameterList > belosList;
  LinearSolver linearSolver;

  // Matrix and right hand side vector of the linear system
  boost::shared_ptr<matrix_Type> SystemMatrix;
  boost::shared_ptr<vector_Type> Rhs;

  // Solution of the linear system
  boost::shared_ptr<vector_Type> Solution;

  // Solution in the quadrature nodes expanding modal functions
  vector_Type solution_3D;
  vector_Type solution_3DOld;

  // Radial velocity on the lateral boundary
  vector_Type urWall;
  vector_Type urWallOld;

  // Radial velocity of the ALE map (assuming that the deformation is only on the radial direction wx = 0 and wtheta = 0)
  //vector_Type wr;

  // Radial displacement on the lateral boundary
  vector_Type etar;
  vector_Type etarOld;

  // Solution in the computational grid
  Vector3D_type ux;
  Vector3D_type ur;
  Vector3D_type utheta;
  Vector3D_type uy;
  Vector3D_type uz;
  Vector3D_type p;

  HiModExporterVtk exporterVtk;

};

}

#endif
