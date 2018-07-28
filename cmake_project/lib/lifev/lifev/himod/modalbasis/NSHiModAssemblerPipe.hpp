//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file NSHiModAssembler.hpp
     @brief This file contains a class for assemble an Navier-Sokes problem on a cylinder with circular section with the Hierarchical Model Reduction method

     @date 03/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
 */
#ifndef __NSHIMODASSEMBLERPIPE__
#define __NSHIMODASSEMBLERPIPE__


#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>
#include <omp.h>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/himod/modalbasis/NSModalSpacePipe.hpp>
#include <lifev/himod/modalbasis/Projector.hpp>
#include <lifev/himod/tools/ReferenceMap.hpp>
#include <lifev/himod/tools/BCstructure.hpp>

#include <lifev/core/array/VectorSmall.hpp>

#pragma GCC diagnostic ignored "-Wconversion"
#include <lifev/core/filter/ExporterVTK.hpp>
#pragma GCC diagnostic warning "-Wconversion"

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <Epetra_SerialComm.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

template< typename mesh_type, typename matrix_type, typename vector_type>
class NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>
{

public:

    //! @name Public Types
    //@{

    typedef    NSHiModAssembler< mesh_type , matrix_type , vector_type, 2 > HMA_type;
    //! Typedef for the map
    typedef MapEpetra                                       map_type;
    //! Typedef for finite element space
    typedef FESpace<mesh_type, map_type>                    fespace_type;
    //! Typedef for a pointer on the FEspace
    typedef boost::shared_ptr<fespace_type>                 fespace_ptrType;
    //! Typedef for a ETFEspace
    typedef ETFESpace<mesh_type, map_type, 1, 1>            etfespace_type;
    //! Typedef for a pointer on the FEspace
    typedef boost::shared_ptr<etfespace_type>               etfespace_ptrType;
    //! Typedef for the modal basis
    typedef NSModalSpacePipe                              modalbasis_type;
    //! Typedef for a pointer on the modal basis
    typedef boost::shared_ptr<modalbasis_type>              modalbasis_ptrType;
    //! Typedef for a pointer on the system matrix
    typedef boost::shared_ptr<matrix_type>                  matrix_ptrType;
    typedef std::vector<std::vector<Real> >              MBMatrix_type;
    //! Typedef for a pointer on the rhs vector
    typedef boost::shared_ptr<vector_type>                  vector_ptrType;
    //! Typedef for a three component vector
    typedef VectorSmall<3>                                  TreDvector_type;
    //! Typedef for a function type
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
    //! Typedef for the functor that compute the projection
    typedef Projector                                       projector_type;
    //! Typedef for a pointer to the functor that compute the projection
    typedef boost::shared_ptr<projector_type>               projector_ptrType;
    //! Typedef for the chrono
    typedef LifeChrono                                      chrono_type;
    //! Typedef for a pointer on the communicator
    typedef typename map_type::comm_ptrtype                 commPtr_Type;
    //! Typedef for the export
    typedef RegionMesh<LinearTetra>                         export_mesh_Type;
    typedef boost::shared_ptr<export_mesh_Type>             export_meshPtr_Type;
    typedef ExporterData<export_mesh_Type>                  exporterData_Type;
    typedef typename exporterData_Type::WhereEnum            WhereEnum;
    typedef typename exporterData_Type::FieldTypeEnum        FieldTypeEnum;
    typedef typename exporterData_Type::FieldRegimeEnum        FieldRegimeEnum;
    typedef VectorEpetra                                    export_vector_type;
    typedef boost::shared_ptr<export_vector_type>            export_vector_ptrType;
    typedef PreconditionerIfpack                            prec_Type;
    //@}


    //! @name Constructor & Destructor
    //@{
    /*!
    The construction of the HiMod Space based on a fespace1D and e ModalSpaceCircular2D. You need comm to construct ETFspace.
    \param fespace a pointer to the fespace1D
    \param modalbasis a pointer to the ModalSpaceCircular2D
    \param Comm a comunicator
    */
    NSHiModAssembler( const fespace_ptrType& ufespace, const fespace_ptrType& pfespace, const modalbasis_ptrType& modalbasis, commPtr_Type& Comm );

    NSHiModAssembler( const HMA_type& HMA, commPtr_Type& Comm  ):
                    M_modalbasis( new modalbasis_type(HMA.modalspaceObj()) ), 
                    M_velocityFespace( HMA.velocityFespace() ),
                    M_etufespace( new etfespace_type ( HMA.velocityFespace()->mesh(), 
                                                       &( HMA.velocityFespace()->refFE() ),
                                                       &( HMA.velocityFespace()->fe().geoMap() ), Comm ) ),
                    M_pressureFespace( HMA.pressureFespace() ),
                    M_etpfespace( new etfespace_type ( HMA.velocityFespace()->mesh(), 
                                                       &( HMA.velocityFespace()->refFE() ),
                                                       &( HMA.velocityFespace()->fe().geoMap() ), Comm ) ){};

    //! empty destructor
    virtual ~NSHiModAssembler() {};
    //@}

    //! @name Public Methods
    //@{
    /*!
    With this function you can set the coefficient of the variational form of the standard 3D ADR problem.
    \param systemmatrix a pointer to the systemmatrix
    \param mu   the diffusion coefficient
    \param beta the 3D-vector defining a, constant, advection field
    \param sigma the reaction coefficient
    */
    
    void addStokesProblem( const matrix_ptrType& systemMatrix,
                            const Real& nu, ReferenceMap& refMap,
                            const Real& t, const Real& alpha );
                            
    void addAdvection( const matrix_ptrType& systemMatrix,
                       const vector_type& adv );
                     
    // Advection field as input parameter ( Navier-Stokes equations ) 
    /*void addNavierStokesProblem( const matrix_ptrType& systemMatrix,
                        const Real& nu,
                        ReferenceMap& refMap, const Real& t, const Real& alpha,
                        const vector_type& advection
                        );*/
    void pressureMassMatrix( const matrix_ptrType& systemMatrix,
                        ReferenceMap& refMap
                        );
    //! @name Public Methods
    //@{
    /*!
    With this function you can set the coefficient of the variational form
    \param systemmatrix a pointer to the systemmatrix
    \param mu   the diffusion coefficient, a function
    \param beta the 3D-vector defining a, constant, advection field a function
    \param sigma the reaction coefficient a function
    */
    /*void addNavierStokesProblem( const matrix_ptrType& systemMatrix,
                        const function_Type& nu,
                        ReferenceMap& refMap, const Real& t, const Real& alpha
    //                    , const vector_type& u_kMinusOne
                                );*/
    /*!
    This method interpolates the Fourier coefficients of f in the nodes of 1D FESpace.
    @param f_interpolated a pointer to the vector which has a block structure, where each block refers to a different mode
    @param f the force term
    */
    void interpolate( const function_Type& fx, const function_Type& fr, const function_Type& ftheta,
                      const vector_ptrType& f_interpolated, const Real& t );
    
    void interpolate( const vector_type& fxrtheta,
                      const vector_ptrType& f_interpolated );

    /*!
    This methods add the rhs for a constante force term
    \param rhs a pointer to the rhs
    \param f the constant force term
    */
    void addrhs( const vector_ptrType& rhs,
                 const Real& fx, const Real& fr, const Real& ftheta );
    /*!
    This methods add the rhs for a generic force term \f$f=f(x,y,z)\f$ which has been already interpolated with the
    interpolate method
    @param rhs a pointer to the rhs
    @param fx_interpolated the x-force term which has been interpolated
    @param fr_interpolated the r-force term which has been interpolated
    @param ftheta_interpolated the theta-force term which has been interpolated
    */
    void addrhs( const vector_ptrType& rhs,
                 const vector_ptrType& f_interpolated, const vector_ptrType& uold_interpolated );
                 
    /*!
    This methods add the dirichlet boundary conditions at inflow.
    \param systemmatrix a pointer to the systemmatrix
    \param rhs a pointer to the rhs
    \param g = g(y,z) the Dirichlet data.

    It computes the fourier coefficients of g, and it imposes them
    essentially via penalization. ( 10^30 * u_kk= 10^30 * g_k )
    */
                            
    // ---------------  non constant radius ( not all the bcs are available )
    void addDirichletBC_xIn( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
                            
    void addDirichletBC_rIn( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
                            
    void addDirichletBC_thetaIn( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
                            
    void addNeumannBC_xIn( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
    
    void addNeumannBC_rIn( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
    
    void addNeumannBC_thetaIn( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );

    void addNeumannBC_xOut( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
    
    void addNeumannBC_rOut( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
    
    void addNeumannBC_thetaOut( const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
                            
    void addDirichletBC_xOut( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );

    void addDirichletBC_rOut( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );
                            
    void addDirichletBC_thetaOut( const matrix_ptrType& systemMatrix,
                            const vector_ptrType& rhs,
                            const function_Type& g, const Real& t );

    // Uneducated version for non constant radius: the last parameter enables the diagonalization of the sine matrix
    void addBC( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const BCdata& bc, const Real& t );
    
    void diagonalizeBlocks( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs );
    
    void diagonalizeBCmatrix( const matrix_ptrType& matrix, const vector_ptrType& rhs );

    //@}

    //!@name getMethods
    //@{
    modalbasis_ptrType modalspace() const
    {
        return M_modalbasis;
    }

    modalbasis_type modalspaceObj() const
    {
        return *M_modalbasis;
    }

    fespace_ptrType velocityFespace() const
    {
        return M_velocityFespace;
    }
    
    fespace_ptrType pressureFespace() const
    {
        return M_pressureFespace;
    }
    //@}

    /*!
        Compute the value of the output of the HiMod problem on the 3D grid made by quadrature nodes and mesh nodes FEM.
    */
    vector_type evaluateBase3DGrid( const vector_type& fun );

    /*!
        Compute the value of a function on the 3D grid made by quadrature nodes and mesh nodes FEM.
    */
    vector_type evaluateBase3DGrid( const function_Type& xVel, const function_Type& rVel, const function_Type& thetaVel,
                                    const function_Type& press, const Real& t );

    /*!
        Compute the L2 norm of the function
    */
    Real normL2( const vector_type& fun, const std::string& solution );

    /*!
        This method evaluates a vector HiMod function (the collect of the HiMod coefficients) in the point (x,y,z). It uses the         generator of basis stored in modalbasis.
    */
    Real evaluateHiModFunc( const vector_ptrType& HMcoef,
                            const Real& x, const Real& r, const Real& theta, const UInt& component );
    std::vector<Real> evaluateHiModFunc( const vector_ptrType& HMcoef,
                                         const std::vector<Real>& x, const std::vector<Real>& y, const std::vector<Real>& z );

    /*!
        Exports a vector on a grid with nx*ny*nz points with optimizations (the grid is structured)
    */
//    void exporterStructuredVTK( const UInt& nx, const UInt& ny, const UInt& nz, const vector_ptrType& fun, const GetPot& dfile, std::string prefix, std::string content = "solution");

    /*!
        Exports a function on a grid with nx*ny*nz points with optimizations (the grid is structured)
    */
//    void exporterFunctionVTK (const UInt& nx, const UInt& ny, const UInt& nz, const function_Type& f, const GetPot& dfile, std::string prefix, std::string content = "function");
private:

    //! A pointer to the modal basis
    modalbasis_ptrType  M_modalbasis;

    //! A pointer to the ETFespace
    etfespace_ptrType   M_etufespace;
    etfespace_ptrType   M_etpfespace;

    //! A pointer to the finite element space
    fespace_ptrType     M_velocityFespace;
    fespace_ptrType     M_pressureFespace;
    
};

//---------------------------------------------
//      IMPLEMENTATION
//---------------------------------------------

//Constructor
template< typename mesh_type, typename matrix_type, typename vector_type>
NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
NSHiModAssembler( const fespace_ptrType& ufespace, const fespace_ptrType& pfespace,
                  const modalbasis_ptrType& modalbasis, commPtr_Type& Comm ) :
                M_modalbasis( modalbasis ),
                M_etufespace( new etfespace_type ( ufespace->mesh(), &( ufespace->refFE() ), &( ufespace->fe().geoMap() ), Comm ) ),
                M_etpfespace( new etfespace_type ( ufespace->mesh(), &( pfespace->refFE() ), &( ufespace->fe().geoMap() ), Comm ) ),
                M_velocityFespace( ufespace ), M_pressureFespace( pfespace )
{}


// Implementation of addNavierStokesProblem for both the constant and the non constant version (TODO).
// Implementation of addrhs for both the constant and the non constant version.
//#include <lifev/himod/modalbasis/HMAaddNSAssembler.hpp>
#include <lifev/himod/modalbasis/HMAaddNSAssemblerPipe.hpp>

// Implementation of addBC and all the related functions (add_Dirichlet/Neumann_In/Out_x/r/theta)
#include <lifev/himod/modalbasis/HMAbcHandlerPipe.hpp>

// Implementation of interpolation functions
#include <lifev/himod/modalbasis/HMAinterpolationPipe.hpp>

// Utilities
#include <lifev/himod/modalbasis/HMAutilitiesPipe.hpp>

} // namespace LifeV

#endif // __NSHIMODASSEMBLERCIRCULAR__


