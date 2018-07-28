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
 *   @file HiModAssembler.hpp
     @brief This file contains a class for assemble an ADR problem with the Hierarchical Model Reduction method

     @date 06/2013
     @author A. Bortolossi <andrea.bortolossi@gmail.com>
     @author M. Aletti <teo.aletti@gmail.com>
 */
#ifndef __HIMODASSEMBLERRECTANGULAR__
#define __HIMODASSEMBLERRECTANGULAR__


#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>
#include <lifev/himod/modalbasis/Projector.hpp>

#include <lifev/core/array/VectorSmall.hpp>

#pragma GCC diagnostic ignored "-Wconversion"
#include <lifev/core/filter/ExporterVTK.hpp>
#pragma GCC diagnostic warning "-Wconversion"

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

template< typename mesh_type, typename matrix_type, typename vector_type>
class HiModAssembler<mesh_type, matrix_type, vector_type, 0>
{

public:

    //! @name Public Types
    //@{

    //! Typedef for the map
    typedef MapEpetra                                       map_type;

    //! Typedef for finite element space
    typedef FESpace<mesh_type, map_type>                     fespace_type;
    //! Typedef for a pointer on the FEspace
    typedef boost::shared_ptr<fespace_type>                 fespace_ptrType;

    //! Typedef for a ETFEspace
    typedef ETFESpace<mesh_type, map_type, 1, 1>              etfespace_type;
    //! Typedef for a pointer on the FEspace
    typedef boost::shared_ptr<etfespace_type>               etfespace_ptrType;

    //! Typedef for the modal basis
    typedef ModalSpaceRectangular                           modalbasis_type;
    //! Typedef for a pointer on the modal basis
    typedef boost::shared_ptr<modalbasis_type>              modalbasis_ptrType;

    //! Typedef for a pointer on the system matrix
    typedef boost::shared_ptr<matrix_type>                  matrix_ptrType;
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
    typedef typename map_type::comm_ptrtype commPtr_Type;

    //! Typedef for the export
    typedef RegionMesh<LinearTetra> export_mesh_Type;

    typedef boost::shared_ptr<export_mesh_Type>                export_meshPtr_Type;
    
    typedef ExporterData<export_mesh_Type>                     exporterData_Type;
    
    typedef typename exporterData_Type::WhereEnum       WhereEnum;
    
    typedef typename exporterData_Type::FieldTypeEnum   FieldTypeEnum;
    
    typedef typename exporterData_Type::FieldRegimeEnum FieldRegimeEnum;
    
    typedef VectorEpetra export_vector_type;
    
    typedef boost::shared_ptr<export_vector_type> export_vector_ptrType;
    //@}


    //! @name Constructor & Destructor
    //@{
    /*!
    The construction of the HiMod Space based on a fespace1D and e modalspace2D. You need comm to construct ETFspace.
    \param fespace a pointer to the fespace1D
    \param modalbasis a pointer to the modalspace2D
    \param Comm a comunicator
    */
    HiModAssembler (const fespace_ptrType& fespace, const modalbasis_ptrType& modalbasis, commPtr_Type& Comm);

    //! empty destructor
    virtual ~HiModAssembler() {};
    //@}

    //! @name Public Methods
    //@{
    /*!
    With this function you can set the coefficient of the variational form of the standard 3D ADR problem.
    \f[-\mu\Delta u+b\cdot\nabla u + \sigma u = f\f]a
    \param systemmatrix a pointer to the systemmatrix
    \param mu   the diffusion coefficient
    \param beta the 3D-vector defining a, constant, advection field
    \param sigma the reaction coefficient
    */
    void addADRProblem ( const matrix_ptrType& systemMatrix,
                         const Real& mu,
                         const TreDvector_type& beta,
                         const Real& sigma);
    //! @name Public Methods
    //@{
    /*!
    With this function you can set the coefficient of the variational form of the standard 3D ADR problem.
    \f[-\mu\Delta u+b\cdot\nabla u + \sigma u = f\f]a
    \param systemmatrix a pointer to the systemmatrix
    \param mu   the diffusion coefficient, a function
    \param beta the 3D-vector defining a, constant, advection field a function
    \param sigma the reaction coefficient a function
    */
    void addADRProblem ( const matrix_ptrType& systemMatrix,
                         const function_Type& mu,
                         const function_Type& beta,
                         const function_Type& sigma);
    /*!
    This method interpolates the Fourier coefficients of f in the nodes of 1D FESpace.
    @param f_interpolated a pointer to the vector which has a block structure, where each block refers to a different mode
    @param f the force term
    */
    void interpolate (   const function_Type& f,
                         const vector_ptrType& f_interpolated);
    /*!
    This methods add the rhs for a constante force term
    \param rhs a pointer to the rhs
    \param f the constant force term
    */
    void addrhs (    const vector_ptrType& rhs,
                     const Real& f);
    /*!
    This methods add the rhs for a generic force term \f$f=f(x,y,z)\f$ which has been already interpolated with the
    interpolate method
    @param rhs a pointer to the rhs
    @param f_interpolated the force term which has been interpolated
    */
    void addrhs (    const vector_ptrType& rhs,
                     const vector_ptrType& f_interpolated);
    /*!
    This method add the rhs for a generic force term \f$f=f(x,y,z)\f$, but it uses a functor. It's slower than the standard
    Addrhs method, but you don't have to call the interpolate method in the main. It's should also be preciser specially
    if you use high order quadrature rule in the FESpace.
    @param rhs a pointer to the rhs
    @param f the force term
    */
    void addrhs_HiPrec (    const vector_ptrType& rhs,
                            const function_Type& f);
    /*!
    This methods add the dirichlet boundary conditions at inflow.
    \param systemmatrix a pointer to the systemmatrix
    \param rhs a pointer to the rhs
    \param g = g(y,z) the Dirichlet data.

    It computes the fourier coefficients of g, and it imposes them
    essentially via penalization. ( 10^30 * u_kk= 10^30 * g_k )
    */
    void addDirichletBC_In ( const matrix_ptrType& systemMatrix,
                             const vector_ptrType& rhs,
                             const function_Type& g);

    void interp_Coefficients (const vector_ptrType& r11, const vector_ptrType& r10, const vector_ptrType& r00, const UInt& j, const UInt& k, const function_Type& mu, const function_Type& beta, const function_Type& sigma);
    //@}

    // get method

    modalbasis_ptrType modalspace()
    {
        return M_modalbasis;
    }

    fespace_ptrType fespace()
    {
        return M_fespace;
    }


    /*!
        Compute the value of the output of the HiMod problem on the 3D grid made by quadrature nodes and mesh nodes FEM.
    */
    vector_type evaluateBase3DGrid (const vector_type& fun);

    /*!
        Compute the value of a function on the 3D grid made by quadrature nodes and mesh nodes FEM.
    */
    vector_type evaluateBase3DGrid (const function_Type& fun);

    /*!
        Compute the L2 norm of the function
    */
    Real normL2 (const vector_type& fun);

    /*!
        This method evaluates a vector HiMod function (the collect of the HiMod coefficients) in the point (x,y,z). It uses the         generator of basis stored in modalbasis.
    */
    Real evaluateHiModFunc (const vector_ptrType& rhs, const Real& x, const Real& y, const Real& z);

    /*!
        Exports a vector on a grid with nx*ny*nz points with optimizations (the grid is structured)
    */
    void exporterStructuredVTK (const UInt& nx, const UInt& ny, const UInt& nz, const vector_ptrType& fun, const GetPot& dfile, std::string prefix, std::string content = "solution");

    /*!
        Exports a function on a grid with nx*ny*nz points with optimizations (the grid is structured)
    */
    void exporterFunctionVTK (const UInt& nx, const UInt& ny, const UInt& nz, const function_Type& f, const GetPot& dfile, std::string prefix, std::string content = "function");
private:

    //! A pointer to the modal basis
    modalbasis_ptrType  M_modalbasis;

    //! A pointer to the ETFespace
    etfespace_ptrType   M_etfespace;

    //! A pointer to the finite element space
    fespace_ptrType     M_fespace;
};

//---------------------------------------------
//      IMPLEMENTATION
//---------------------------------------------

//Constructor
template< typename mesh_type, typename matrix_type, typename vector_type>
HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
HiModAssembler (const fespace_ptrType& fespace, const modalbasis_ptrType& modalbasis, commPtr_Type& Comm) :
    M_modalbasis (modalbasis),
    M_etfespace ( new etfespace_type ( fespace->mesh(), & (fespace->refFE() ), & (fespace->fe().geoMap() ), Comm) ),
    M_fespace (fespace)
{}

//AddADRAssembler
/*
    Assemble the matrix of the problem, for the moment with costant coefficients (mu,beta,gamma)
*/
template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addADRProblem (const matrix_ptrType& systemMatrix, const Real& mu, const TreDvector_type& beta, const Real& sigma)
{

    //Cycling on sub-blocks for each pair of frequency
    for (UInt k = 0; k < M_modalbasis->mtot(); ++k)
    {
        for (UInt j = 0; j < M_modalbasis->mtot(); ++j)
        {

            //k refers to the test function, j refers to the solution

            //For each blocks compute the integral coefficient on the slice
            //std::vector<Real> Coeff(5,0.0);
            VectorSmall<5> Coeff;
            //j refers to the solution k refers to the test function
            Coeff[0] = M_modalbasis->compute_PhiPhi (j, k);     //mu sigma beta1
            Coeff[1] = M_modalbasis->compute_DyPhiPhi (j, k);   //beta2
            Coeff[2] = M_modalbasis->compute_DyPhiDyPhi (j, k); //mu
            Coeff[3] = M_modalbasis->compute_DzPhiPhi (j, k);   //beta3
            Coeff[4] = M_modalbasis->compute_DzPhiDzPhi (j, k); //mu
            //Assemble the (j,k)1D FEMproblem

            /*
                 Compute the coefficient which refers to the boundary contribute. Note that in the case of Dirichlet
                the coefficiente automatically become zero. In the case of Neumann condition chi_y or / and chi_z are equal
                to zero.
            */

            VectorSmall<2> Boundary;
            Boundary[0] = M_modalbasis->compute_y_PhiPhi (j, k);
            Boundary[1] = M_modalbasis->compute_z_PhiPhi (j, k);

            UInt p_j = (M_modalbasis->eigenvalues (j).p - 1);
            UInt p_k = (M_modalbasis->eigenvalues (k).p - 1);
            UInt q_j = (M_modalbasis->eigenvalues (j).q - 1);
            UInt q_k = (M_modalbasis->eigenvalues (k).q - 1);

            Real Ly = M_modalbasis->Ly();
            Real Lz = M_modalbasis->Lz();

            Real etaj_Ly = 1. / sqrt (Ly) * M_modalbasis->phiy (p_j, M_modalbasis->qrY().nbQuadPt() - 1);
            Real etaj_0  = 1. / sqrt (Ly) * M_modalbasis->phiy (p_j, 0);

            Real xij_Lz  = 1. / sqrt (Lz) * M_modalbasis->phiz (q_j, M_modalbasis->qrZ().nbQuadPt() - 1);
            Real xij_0   = 1. / sqrt (Lz) * M_modalbasis->phiz (q_j, 0);

            Real etak_Ly = 1. / sqrt (Ly) * M_modalbasis->phiy (p_k, M_modalbasis->qrY().nbQuadPt() - 1);
            Real etak_0  = 1. / sqrt (Ly) * M_modalbasis->phiy (p_k, 0);

            Real xik_Lz  = 1. / sqrt (Lz) * M_modalbasis->phiz (q_k, M_modalbasis->qrZ().nbQuadPt() - 1);
            Real xik_0   =  1. / sqrt (Lz) * M_modalbasis->phiz (q_k, 0);

            // Maybe in case of different basis this part have to change!!
            Real chi_y = M_modalbasis->gbY()->chi();
            Real chi_z = M_modalbasis->gbZ()->chi();
            /*

                 ___3
              4 |   |
                |___| 2
                  1
            z
            ^
            |
            *--> y

            std::cout<<"etaj_Ly "<<etaj_Ly<<std::endl;
            std::cout<<"etaj_0 "<<etaj_0<<std::endl;
            std::cout<<"etak_Ly "<<etak_Ly<<std::endl;
            std::cout<<"etak_0 "<<etak_0<<std::endl;
            std::cout<<"xij_Lz "<<xij_Lz<<std::endl;
            std::cout<<"xij_0 "<<xij_0<<std::endl;
            std::cout<<"xik_Lz "<<xik_Lz<<std::endl;
            std::cout<<"xik_0 "<<xik_0<<std::endl;
            std::cout<<"Chi Y "<<chi_y<<std::endl;
            std::cout<<"Chi Z "<<chi_z<<std::endl;
            exit(1);
            */
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate ( elements (M_etfespace->mesh() ),
                            M_fespace->qr(),
                            M_etfespace,
                            M_etfespace,

                            mu * Coeff[0] * dot ( grad (phi_i) , grad (phi_j) )
                            + beta[0]*Coeff[0] * dot ( grad (phi_j) , value (oneVector) ) *phi_i
                            + ( mu * (Coeff[2] + Coeff[4])
                                + beta[1]*Coeff[1] + beta[2]*Coeff[3]
                                + sigma * Coeff[0]
                                + chi_z * Boundary[0] * (xij_Lz * xik_Lz + xij_0 * xik_0)
                                + chi_y * Boundary[1] * (etaj_Ly * etak_Ly + etaj_0 * etak_0) ) * phi_i * phi_j


                          )
                        >> (systemMatrix->block (k, j) );
            }

        }
    }

    return;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addADRProblem ( const matrix_ptrType& systemMatrix,
                const function_Type& mu,
                const function_Type& beta,
                const function_Type& sigma)
{
    std::cout << "Non constant version, might be very slow.." << std::endl;
    //TODO std::cout<<"We should first evaluate the coefficients on the quadrature rule.."<<std::endl;

    //Cycling on sub-blocks for each pair of frequency
    for (UInt k = 0; k < M_modalbasis->mtot(); ++k)
    {
        for (UInt j = 0; j < M_modalbasis->mtot(); ++j)
        {

            //Assemble the (j,k)1D FEMproblem

            vector_ptrType R11 (new vector_type (M_fespace->map(), Repeated) );
            vector_ptrType R10 (new vector_type (M_fespace->map(), Repeated) );
            vector_ptrType R00 (new vector_type (M_fespace->map(), Repeated) );

            interp_Coefficients (R11, R10, R00, j, k, mu, beta, sigma);

            /*
                 Compute the coefficient which refers to the boundary contribute. Note that in the case of Dirichlet
                the coefficiente automatically become zero. In the case of Neumann condition chi_y or / and chi_z are equal
                to zero.
            */

            VectorSmall<2> Boundary;
            Boundary[0] = M_modalbasis->compute_y_PhiPhi (j, k);
            Boundary[1] = M_modalbasis->compute_z_PhiPhi (j, k);

            UInt p_j = (M_modalbasis->eigenvalues (j).p - 1);
            UInt p_k = (M_modalbasis->eigenvalues (k).p - 1);
            UInt q_j = (M_modalbasis->eigenvalues (j).q - 1);
            UInt q_k = (M_modalbasis->eigenvalues (k).q - 1);

            Real Ly = M_modalbasis->Ly();
            Real Lz = M_modalbasis->Lz();

            Real etaj_Ly = 1. / sqrt (Ly) * M_modalbasis->phiy (p_j, M_modalbasis->qrY().nbQuadPt() - 1);
            Real etaj_0  = 1. / sqrt (Ly) * M_modalbasis->phiy (p_j, 0);

            Real xij_Lz  = 1. / sqrt (Lz) * M_modalbasis->phiz (q_j, M_modalbasis->qrZ().nbQuadPt() - 1);
            Real xij_0   = 1. / sqrt (Lz) * M_modalbasis->phiz (q_j, 0);

            Real etak_Ly = 1. / sqrt (Ly) * M_modalbasis->phiy (p_k, M_modalbasis->qrY().nbQuadPt() - 1);
            Real etak_0  = 1. / sqrt (Ly) * M_modalbasis->phiy (p_k, 0);

            Real xik_Lz  = 1. / sqrt (Lz) * M_modalbasis->phiz (q_k, M_modalbasis->qrZ().nbQuadPt() - 1);
            Real xik_0   =  1. / sqrt (Lz) * M_modalbasis->phiz (q_k, 0);

            // Maybe in case of different basis this part have to change!!
            Real chi_y = M_modalbasis->gbY()->chi();
            Real chi_z = M_modalbasis->gbZ()->chi();
            /*

                 ___3
              4 |   |
                |___| 2
                  1
            z
            ^
            |
            *--> y

            std::cout<<"etaj_Ly "<<etaj_Ly<<std::endl;
            std::cout<<"etaj_0 "<<etaj_0<<std::endl;
            std::cout<<"etak_Ly "<<etak_Ly<<std::endl;
            std::cout<<"etak_0 "<<etak_0<<std::endl;
            std::cout<<"xij_Lz "<<xij_Lz<<std::endl;
            std::cout<<"xij_0 "<<xij_0<<std::endl;
            std::cout<<"xik_Lz "<<xik_Lz<<std::endl;
            std::cout<<"xik_0 "<<xik_0<<std::endl;
            std::cout<<"Chi Y "<<chi_y<<std::endl;
            std::cout<<"Chi Z "<<chi_z<<std::endl;
            exit(1);
            */
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate ( elements (M_etfespace->mesh() ),
                            M_fespace->qr(),
                            M_etfespace,
                            M_etfespace,
                            value (M_etfespace, *R11) * dot ( grad (phi_i) , grad (phi_j) )
                            + value (M_etfespace, *R10) * dot ( grad (phi_j) , value (oneVector) ) *phi_i
                            + value (M_etfespace, *R00) * phi_i * phi_j
							+ (chi_z * Boundary[0] * (xij_Lz * xik_Lz + xij_0 * xik_0)
                            +  chi_y * Boundary[1] * (etaj_Ly * etak_Ly + etaj_0 * etak_0) ) * phi_i * phi_j
                          )
                        >> (systemMatrix->block (k, j) );
            }

        }
    }
    return;
}




template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
interp_Coefficients (const vector_ptrType& r11, const vector_ptrType& r10, const vector_ptrType& r00, const UInt& j, const UInt& k, const function_Type& mu, const function_Type& beta, const function_Type& sigma)
{
    // First, we build a "quadrature" that consists in the nodes (0 weight)
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_fespace->refFEPtr()->shape() ), M_fespace->refFEPtr()->shape() );
    interpQuad.setPoints (M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_fespace->refFEPtr()->nbDof(), 0) );

    // Then, we define a currentFE with nodes on the reference nodes
    CurrentFE interpCFE (M_fespace->refFE(), getGeometricMap (*M_fespace->mesh() ), interpQuad);

    // Some constants
    UInt totalNumberElements (M_fespace->mesh()->numElements() );
    UInt numberLocalDof (M_fespace->dofPtr()->numLocalDof() );//Local to the element        Esempio P1  2
    //UInt numberTotalDof (M_fespace->dofPtr()->numTotalDof() );//Total           Esempio P1 N

    // Storage for the values
    std::vector<Real> nodalValues11 (numberLocalDof, 0);
    std::vector<Real> nodalValues10 (numberLocalDof, 0);
    std::vector<Real> nodalValues00 (numberLocalDof, 0);
    std::vector<Real> FEValues11 (numberLocalDof, 0);
    std::vector<Real> FEValues10 (numberLocalDof, 0);
    std::vector<Real> FEValues00 (numberLocalDof, 0);

    // Do the loop over the cells

    //Esempio N nodi totalNumberElements = N-1
    for (UInt iterElement (0); iterElement < totalNumberElements; ++iterElement)
    {
        // We update the CurrentFE so that we get the coordinates of the nodes
        interpCFE.update (M_fespace->mesh()->element (iterElement), UPDATE_QUAD_NODES);

        // Loop over the degrees of freedom (= quadrature nodes)
        for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
        {
            // Store the nodal value of the fourier coefficient const UInt& j, const UInt& k, const function_Type& mu, const Real& xinterpCFE.quadNode (iterDof, 0)
            nodalValues11[iterDof] =  M_modalbasis->compute_R11 (j, k, mu, interpCFE.quadNode (iterDof, 0) );
            nodalValues10[iterDof] =  M_modalbasis->compute_R10 (j, k, beta, interpCFE.quadNode (iterDof, 0) );
            nodalValues00[iterDof] =  M_modalbasis->compute_R00 (j, k, mu, beta, sigma, interpCFE.quadNode (iterDof, 0) );
        }

        // Transform the nodal values in FE values
        // Perchè ad esempio nel caso P1 tutti i nodi interni gli ho calcolati il doppio delle volte necessario

        FEValues11 = M_fespace->refFEPtr()->nodalToFEValues (nodalValues11); //In P1 case it's useless
        FEValues10 = M_fespace->refFEPtr()->nodalToFEValues (nodalValues10); //In P1 case it's useless
        FEValues00 = M_fespace->refFEPtr()->nodalToFEValues (nodalValues00); //In P1 case it's useless
        // Then on the dimension of the FESpace (scalar field vs vectorial field)
        for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
        {
            // Find the ID of the considered DOF
            ID globalDofID (M_fespace->dofPtr()->localToGlobalMap (iterElement, iterDof) );
            // Compute the value of the function and set it
            //#pragma GCC diagnostic ignored "-Wconversion"
            r11->setCoefficient (globalDofID, FEValues11[iterDof]);
            r10->setCoefficient (globalDofID, FEValues10[iterDof]);
            r00->setCoefficient (globalDofID, FEValues00[iterDof]);
            //#pragma GCC diagnostic warning "-Wconversion"
        }
    }
}
//Note:
//1) Added get methods in FESpace.hpp to get res FE pointer
//2) There's not a function to add a coefficient to a VectorEpetraStructured, it should be into VectorStructuredView.hpp like
//   in MatrixStructuredView.hpp
//3) Problem with vect.setCoefficient if you use parallel mode?
//4) MeshPartitioner.hpp contiene dei float!!!!
template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
interpolate (const function_Type& f, const vector_ptrType& f_interpolated)
{
    // First, we build a "quadrature" that consists in the nodes (0 weight)
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_fespace->refFEPtr()->shape() ), M_fespace->refFEPtr()->shape() );
    interpQuad.setPoints (M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_fespace->refFEPtr()->nbDof(), 0) );

    // Then, we define a currentFE with nodes on the reference nodes
    CurrentFE interpCFE (M_fespace->refFE(), getGeometricMap (*M_fespace->mesh() ), interpQuad);

    // Some constants
    UInt totalNumberElements (M_fespace->mesh()->numElements() );
    UInt numberLocalDof (M_fespace->dofPtr()->numLocalDof() );//Local to the element
    UInt numberTotalDof (M_fespace->dofPtr()->numTotalDof() );//Total

    // Storage for the values
    std::vector<Real> nodalValues (numberLocalDof, 0);
    std::vector<Real> FEValues (numberLocalDof, 0);

    //Do the loop over the frequencies
    for (UInt k (0); k < M_modalbasis->mtot() ; ++k)
    {
        // Do the loop over the cells
        for (UInt iterElement (0); iterElement < totalNumberElements; ++iterElement)
        {
            // We update the CurrentFE so that we get the coordinates of the nodes
            interpCFE.update (M_fespace->mesh()->element (iterElement), UPDATE_QUAD_NODES);

            // Loop over the degrees of freedom (= quadrature nodes)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                // Store the nodal value of the fourier coefficient
                nodalValues[iterDof] =  M_modalbasis->fourierCoeffPointWise (interpCFE.quadNode (iterDof, 0), f, k);
            }

            // Transform the nodal values in FE values
            FEValues = M_fespace->refFEPtr()->nodalToFEValues (nodalValues); //In P1 case it's uselesss

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                // Find the ID of the considered DOF
                ID globalDofID (M_fespace->dofPtr()->localToGlobalMap (iterElement, iterDof) );
                // Compute the value of the function and set it
                //#pragma GCC diagnostic ignored "-Wconversion"
                f_interpolated->setCoefficient (globalDofID + k * numberTotalDof, FEValues[iterDof]);
                //#pragma GCC diagnostic warning "-Wconversion"
            }
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addrhs (const vector_ptrType& rhs, const Real& f)
{
    //Cycling on sub-blocks for each pair of frequency
    for (UInt k = 0; k < M_modalbasis->mtot(); ++k)
    {
        Real Coeff = M_modalbasis->compute_Phi (k);

        {
            using namespace ExpressionAssembly;

            integrate ( elements (M_etfespace->mesh() ),
                        M_fespace->qr(),
                        M_etfespace,

                        f * Coeff * phi_i
                      )
                    >> (rhs->block (k) );
        }

    }
    return;
}

//
//1) Verify efficency of assignment on the little rhs_k
template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addrhs (const vector_ptrType& rhs, const vector_ptrType& f_interpolated)
{


    //Cycling on sub-blocks for each pair of frequency
    for (UInt k = 0; k < M_modalbasis->mtot(); ++k)
    {
        VectorEpetra rhs_k (M_fespace->map(), Repeated);
        for (UInt s (0); s < M_fespace->dofPtr()->numTotalDof(); ++s)
        {
            rhs_k[s] = (*f_interpolated) [s + f_interpolated->block (k)->firstIndex()];
        }


        {
            using namespace ExpressionAssembly;

            integrate  ( elements (    M_etfespace->mesh() ),
                         M_fespace->qr(),
                         M_etfespace,
                         value (M_etfespace, rhs_k) *phi_i
                       )
                    >> (rhs->block (k) );
        }

    }

    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addrhs_HiPrec ( const vector_ptrType& rhs,  const function_Type& f)
{
    projector_ptrType F (new projector_type (M_modalbasis) );
    F->setFunction (f);

    for (UInt k = 0; k < M_modalbasis->mtot(); ++k)
    {
        F->setMode (k);
        {
            using namespace ExpressionAssembly;

            integrate ( elements (M_etfespace->mesh() ),
                        M_fespace->qr(),
                        M_etfespace,
                        eval (F, X) * phi_i
                      )
                    >> (rhs->block (k) );
        }

    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
addDirichletBC_In (const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g)
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->fourierCoefficients (g);
    UInt dof = M_etfespace->dof().numTotalDof();
    for ( UInt j = 0; j < M_modalbasis->mtot(); ++j)
    {
        systemMatrix->setCoefficient ( j * dof, j * dof, 1e+30 );
        rhs->setCoefficient (j * dof, 1e+30 * FCoefficients_g[j]);
    }
}

// ----------------------  UTILITY FUNCTIONS -------------------------------------

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
evaluateBase3DGrid (const vector_type& fun)
{

    UInt nquadY = M_modalbasis->qrY().nbQuadPt();
    UInt nquadZ = M_modalbasis->qrZ().nbQuadPt();

    DOF DataFESpace (M_fespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);

    MapEpetra Map (nquadY * nquadZ * ndofFE, Comm);
    vector_type fcoeff ( Map, Unique );
    fcoeff *= 0.0;

    /*
         sum(j,s){ u_js * phi(y,z)_j * psi(x)_s } = sum(j,s){ u_js * alpha(y)_p(j) * beta(z)_q(j) * psi(x)_s }

        Il vettore soluzione (u_js) può essere rimappato con l'indice k = j*ndofFe + s

        Fisso s:        ovvero ciclo continuamente su una slice, su di essa insiste una sola funzione di base fem e vale 1
                    questo in generale anche nel caso P2 o più (vero?).
        Fisso j:        data una frequenza calcolo il contributo su ogni nodo della slice
        Itero su j: aggiungo ai nodi tutti i contributi di tutte le frequenze

        Itero su s: passo alla slice successiva

        7 8 9
        4 5 6
        1 2 3   ->  1 2 3 4 5 6 7 8 9
    */

    UInt p_j;
    UInt q_j;

    Real normy = 1.0 / std::sqrt (M_modalbasis->Ly() );
    Real normz = 1.0 / std::sqrt (M_modalbasis->Lz() );

    for (UInt s (0); s < ndofFE; ++s) //on FE nodes
    {
        for (UInt j (0); j < M_modalbasis->mtot(); ++j) //Ciclying over all modes contribute
        {
            p_j = (M_modalbasis->eigenvalues (j).p - 1);
            q_j = (M_modalbasis->eigenvalues (j).q - 1);
            for (UInt nz (0); nz < nquadZ; ++nz) //on quadrature node nquadZ
            {
                for (UInt ny (0); ny < nquadY; ++ny) //on quadrature node nqadY
                {

                    fcoeff[ny + nz * nquadY + s * nquadZ * nquadY] +=     fun[j * ndofFE + s] *
                                                                          M_modalbasis->phiy (p_j, ny) * normy *
                                                                          M_modalbasis->phiz (q_j, nz) * normz;
                }
            }
        }
    }

    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
evaluateBase3DGrid (const function_Type& fun)
{
    UInt nquadY = M_modalbasis->qrY().nbQuadPt();
    UInt nquadZ = M_modalbasis->qrZ().nbQuadPt();

    DOF DataFESpace (M_fespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);

    MapEpetra Map (nquadY * nquadZ * ndofFE, Comm);
    vector_type fcoeff ( Map, Unique );
    fcoeff *= 0.0;

    Real x, y, z;

    for (UInt s (0); s < (ndofFE - 1); ++s) //on FE nodes in realtà devo ciclare sugli elementi e non sui nodi
    {
        //TODO if you use P2 you have to modify this part of the code
        QuadratureRule interpQuad;
        interpQuad.setDimensionShape (shapeDimension (M_fespace->refFEPtr()->shape() ), M_fespace->refFEPtr()->shape() );
        interpQuad.setPoints (M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_fespace->refFEPtr()->nbDof(), 0) );
        CurrentFE interpCFE (* ( M_fespace->refFEPtr() ), getGeometricMap (* ( M_fespace->mesh() ) ), interpQuad);
        interpCFE.update (M_fespace->mesh()->element (s), UPDATE_QUAD_NODES);

        x = interpCFE.quadNode (0, 0);  //0 -> nodo di sinistra dell'elemento
        //0 -> legato alla dimensione

        for (UInt nz (0); nz < nquadZ; ++nz) //on quadrature node nquadZ
        {

            z = M_modalbasis->Lz() * M_modalbasis->qrZ().quadPointCoor (nz, 0);

            for (UInt ny (0); ny < nquadY; ++ny) //on quadrature node nqadY
            {
                y = M_modalbasis->Ly() * M_modalbasis->qrY().quadPointCoor (ny, 0);

                fcoeff[ny + nz * nquadY + s * nquadZ * nquadY] = fun (0, x, y, z, 0);

                if (s == (ndofFE - 2) ) //se sono sull'ultimo elemento
                {
                    x = interpCFE.quadNode (1, 0);  //1 -> nodo di destra dell'elemento
                    //0 -> legato alla dimensione
                    fcoeff[ny + nz * nquadY + (s + 1) *nquadZ * nquadY] +=     fun (0, x, y, z, 0);
                }
            }
        }
    }

    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
Real HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
normL2 (const vector_type& fun)
{
    Real norm = 0.0;
    UInt nquadY = M_modalbasis->qrY().nbQuadPt();
    UInt nquadZ = M_modalbasis->qrZ().nbQuadPt();

    DOF DataFESpace (M_fespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_fespace->refFEPtr()->shape() ), M_fespace->refFEPtr()->shape() );
    interpQuad.setPoints (M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_fespace->refFEPtr()->nbDof(), 0) );
    CurrentFE interpCFE (* ( M_fespace->refFEPtr() ), getGeometricMap (* ( M_fespace->mesh() ) ), interpQuad);
    interpCFE.update (M_fespace->mesh()->element (0), UPDATE_QUAD_NODES);

    Real h   = interpCFE.quadNode (1, 0) - interpCFE.quadNode (0, 0);
    Real w_x = 0.5 * h;
    Real w_y;
    Real w_z;

    //Ciclo sulla griglia proprio nello stesso modo in cui ho ciclato per assegnare i valori in funCoeff_3D
    for (UInt s (0); s < (ndofFE - 1); ++s) //ndofFE
    {

        for (UInt nz (0); nz < nquadZ; ++nz) //on quadrature node nquadZ
        {
            for (UInt ny (0); ny < nquadY; ++ny) //on quadrature node nqadY
            {

                w_y = M_modalbasis->qrY().weight (ny);
                w_z = M_modalbasis->qrZ().weight (nz);
                norm +=     (fun[ny + nz * nquadY + s * nquadZ * nquadY] * fun[ny + nz * nquadY + s * nquadZ * nquadY] +
                             fun[ny + nz * nquadY + (s + 1) * nquadZ * nquadY] * fun[ny + nz * nquadY + (s + 1) * nquadZ * nquadY]) *
                            w_y * w_z *
                            M_modalbasis->Ly() * M_modalbasis->Lz() *
                            w_x;
            }

        }
    }

    return sqrt (norm);

}


template< typename mesh_type, typename matrix_type, typename vector_type>
Real HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
evaluateHiModFunc (const vector_ptrType& fun, const Real& x, const Real& y, const Real& z)
{
    Real yh = y / M_modalbasis->Ly();
    Real zh = z / M_modalbasis->Lz();


    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_fespace->refFEPtr()->shape() ), M_fespace->refFEPtr()->shape() );
    interpQuad.setPoints (M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_fespace->refFEPtr()->nbDof(), 0) );
    CurrentFE interpCFE (* ( M_fespace->refFEPtr() ), getGeometricMap (* ( M_fespace->mesh() ) ), interpQuad);
    interpCFE.update (M_fespace->mesh()->element (0), UPDATE_QUAD_NODES);
    DOF DataFESpace (M_fespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    Real h   = interpCFE.quadNode (1, 0) - interpCFE.quadNode (0, 0);

    UInt i = std::floor ( x / h );

    Real fem_i   = ( h * ( i + 1 ) - x ) / h;
    Real fem_ip1 = ( x - i * h ) / h;

    for (UInt m (0); m < M_modalbasis->mtot(); ++m)
    {
        return  ( (*fun) [ i + m * ndofFE] * fem_i +  (*fun) [ i + 1 + m * ndofFE] * fem_ip1 ) *
                M_modalbasis->gbY()->evalSinglePoint (M_modalbasis->eigenvalues (m).wp, yh) / sqrt (M_modalbasis->M_Ly) *
                M_modalbasis->gbZ()->evalSinglePoint (M_modalbasis->eigenvalues (m).wq, zh) / sqrt (M_modalbasis->M_Lz);
    }
}

//-------------------------------- EXPORT -----------------------------------------------------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
exporterStructuredVTK (const UInt& nx, const UInt& ny, const UInt& nz, const vector_ptrType& fun, const GetPot& dfile, std::string prefix, std::string content)
{
    Real Lx = M_fespace->mesh()->lastPoint().x() - M_fespace->mesh()->firstPoint().x();
    Real Ly = M_modalbasis->Ly();
    Real Lz = M_modalbasis->Lz();

    UInt mtot = M_modalbasis->mtot();

    DOF DataFESpace (M_fespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    Real h = Lx / (ndofFE - 1); //Only P1

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
    boost::shared_ptr< export_mesh_Type > MeshPtr ( new RegionMesh<LinearTetra> ( Comm ) );
    regularMesh3D ( *MeshPtr, 1, nx - 1, ny - 1, nz - 1, false,
                    Lx, Ly , Lz ,
                    0.0,  0.0,  0.0);

    boost::shared_ptr<FESpace< export_mesh_Type, MapEpetra > > fespace ( new FESpace< export_mesh_Type, MapEpetra > (MeshPtr, "P1", 1, Comm) );
    export_vector_ptrType fun_out ( new export_vector_type (fespace->map(), Repeated) );
    (*fun_out) *= 0 ;

    std::vector<Real> phiy;
    phiy.resize (ny);

    std::vector<Real> phiz;
    phiz.resize (nz);

    std::vector<Real> psi;
    psi.resize (nx);

    // For every frequency compute all the contributes one time only
    for (UInt m (0); m < mtot; ++m)
    {
        //i <-> x      j <-> y     k <-> z
        for (UInt j (0); j < ny; ++j)
        {
            Real y = MeshPtr->peak (nx * j).y() ; //0 + j*nx + 0*ny*nx
            phiy[ j ] = M_modalbasis->gbY()->evalSinglePoint (M_modalbasis->eigenvalues (m).wp, y / Ly ) / sqrt (Ly);
        }

        for (UInt k (0); k < nz; ++k)
        {
            Real z = MeshPtr->peak (nx * ny * k).z() ; //0 + 0*nx + k*ny*nx
            phiz[ k ] = M_modalbasis->gbZ()->evalSinglePoint (M_modalbasis->eigenvalues (m).wq, z / Lz ) / sqrt (Lz);
        }

        for ( UInt i (0); i < nx; ++i )
        {
            Real x = MeshPtr->peak (i).x();

            UInt id = static_cast<UInt> ( std::floor ( x / h ) );

            Real fem_i   = ( h * ( id + 1 ) - x ) / h;//only P1
            Real fem_ip1 = ( x - id * h ) / h;

            psi[i] = ( (*fun) [ id + m * ndofFE] * fem_i  +  (*fun) [ id + 1 + m * ndofFE] * fem_ip1 );
        }

        for ( UInt k (0); k < nz; ++k )
        {
            for ( UInt j (0); j < ny; ++j )
            {
                for ( UInt i (0); i < nx; ++i )
                {

                    (*fun_out) [i + j * nx + k * ny * nx] +=  psi[i] * phiy[j] * phiz[k];
                }
            }
        }
    }
    //End of the evaluation phase
    //Now we define an exporter and save the function
    ExporterVTK<export_mesh_Type> exporter ( dfile, prefix);
    exporter.addVariable (exporterData_Type::ScalarField, content, fespace, fun_out, 0, exporterData_Type::SteadyRegime, exporterData_Type::Node );
    exporter.setMeshProcId ( MeshPtr , 0 );

    exporter.postProcess (0);
}


template< typename mesh_type, typename matrix_type, typename vector_type>
void HiModAssembler<mesh_type, matrix_type, vector_type, 0>::
exporterFunctionVTK (const UInt& nx, const UInt& ny, const UInt& nz, const function_Type& f, const GetPot& dfile, std::string prefix, std::string content)
{
    Real Lx = M_fespace->mesh()->lastPoint().x() - M_fespace->mesh()->firstPoint().x();
    Real Ly = M_modalbasis->Ly();
    Real Lz = M_modalbasis->Lz();
    Real dx = Lx / (nx-1);
    Real dy = Ly / (ny-1);
    Real dz = Lz / (nz-1);

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
    boost::shared_ptr< export_mesh_Type > MeshPtr ( new RegionMesh<LinearTetra> ( Comm ) );
    regularMesh3D ( *MeshPtr, 1, nx - 1, ny - 1, nz - 1, false,
                    Lx, Ly , Lz ,
                    0.0,  0.0,  0.0);

    boost::shared_ptr<FESpace< export_mesh_Type, MapEpetra > > fespace ( new FESpace< export_mesh_Type, MapEpetra > (MeshPtr, "P1", 1, Comm) );
    export_vector_ptrType fun_out ( new export_vector_type (fespace->map(), Repeated) );
    (*fun_out) *= 0 ;

    for ( UInt k (0); k < nz; ++k )
    {
        for ( UInt j (0); j < ny; ++j )
        {
            for ( UInt i (0); i < nx; ++i )
            {

                (*fun_out) [i + j * nx + k * ny * nx] =  f (0, i * dx, j * dy, k * dz, 0);
            }
        }
    }

    ExporterVTK<export_mesh_Type> exporter ( dfile, prefix);
    exporter.addVariable (   exporterData_Type::ScalarField, content, fespace, fun_out,
                             0, exporterData_Type::SteadyRegime, exporterData_Type::Node );
    exporter.setMeshProcId ( MeshPtr , 0 );

    exporter.postProcess (0);
}

} // namespace LifeV

#endif // __MODALBASIS__
