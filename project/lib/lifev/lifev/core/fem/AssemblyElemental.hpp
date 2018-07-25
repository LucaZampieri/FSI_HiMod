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
    @file
    @brief File containing the procedures for the local assembly of the differential operators

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    To delete:
    L. 195 : stiff_curl
    L. 271 : grad_div
    L. 326 : stab_stokes
    L. 434 : grad_ss
    L. 493 : source_fhn
    L. 617 : choldc
    L. 618 : cholsl
 */


#ifndef _ELEMOPER_H_INCLUDED
#define _ELEMOPER_H_INCLUDED


#include <boost/numeric/ublas/matrix.hpp>


#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/core/fem/CurrentFEManifold.hpp>
#include <lifev/core/fem/ReferenceFEHybrid.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/DOF.hpp>

namespace LifeV
{
//! @name Public typedefs
//@{
typedef boost::numeric::ublas::matrix<Real> Matrix;
typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrix;
//@}

/*! /namespace AssemblyElemental

  This namespace is specially designed to contain the elementary
  operations (corresponding to differential operators) that build
  the local contributions to be used in the assembly procedures.

 */
namespace AssemblyElemental
{
//! @name Public typedefs
//@{
//! Use the portable syntax of the boost function
typedef boost::function5 < const Real,  const Real&, const Real&,
        const Real&, const Real&, const ID&  > function_Type;
//@}

//! Elementary mass for constant mass coefficient
/*!
  This function assembles the local mass matrix when the mass coefficient is constant.

  @param localMass The local matrix to be filled (not cleaned by this function)
  @param massCFE The currentFE structure already updated for the assembly. It requires
  phi and wDetJacobian to be accessible.
  @param coefficient The mass coefficient
  @param fieldDim The dimension of the FE space (scalar/vectorial)
 */
void mass (MatrixElemental& localMass,
           const CurrentFE& massCFE,
           const Real& coefficient,
           const UInt& fieldDim);

//! Elementary stiffness for constant coefficient
/*!
  This function assembles the local stiffness matrix when the coefficient is constant.

  @param localStiff The local matrix to be filled (not cleaned by this function)
  @param stiffCFE The currentFE structure already updated for the assembly. It requires
  dphi and wDetJacobian to be accessible.
  @param coefficient The coefficient
  @param fieldDim The dimension of the FE space (scalar/vectorial)
 */
void stiffness (MatrixElemental& localStiff,
                const CurrentFE& stiffCFE,
                const Real& coefficient,
                const UInt& fieldDim);


//! Interpolation procedure
template<typename localVector, typename globalVector>
void interpolate (localVector& localValues,
                  const CurrentFE& interpCFE,
                  const UInt& spaceDim,
                  const DOF& betaDof,
                  const UInt& elementID,
                  const globalVector& beta)
{
    const UInt nbQuadPt (interpCFE.nbQuadPt() );
    const UInt nbFEDof (interpCFE.nbFEDof() );
    const UInt totalDof (betaDof.numTotalDof() );

    for (UInt iterDim (0); iterDim < spaceDim; ++iterDim)
    {
        // Loop on the quadrature nodes
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            localValues[iQuadPt][iterDim] = 0.0;

            // Loop over the basis functions
            for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
            {
                localValues[iQuadPt][iterDim] +=
                    beta[ betaDof.localToGlobalMap (elementID, iDof) + iterDim * totalDof]
                    * interpCFE.phi (iDof, iQuadPt);
            }
        }
    }
}

//! Interpolation of the gradient
template<typename localVector, typename globalVector>
void interpolateGradient (localVector& localGradient,
                          const CurrentFE& interpCFE,
                          const UInt& spaceDim,
                          const DOF& betaDof,
                          const UInt& elementID,
                          const globalVector& beta)
{
    const UInt nbQuadPt (interpCFE.nbQuadPt() );
    const UInt nbFEDof (interpCFE.nbFEDof() );
    const UInt totalDof (betaDof.numTotalDof() );

    for (UInt iterDim (0); iterDim < spaceDim; ++iterDim)
    {
        // Loop on the quadrature nodes
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {

            for (UInt jDim (0); jDim < nDimensions; ++jDim)
                //for ( jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
            {
                localGradient[ iQuadPt ][ iterDim ][ jDim ] = 0.0;
                for ( UInt i = 0; i < nbFEDof; ++i )
                    localGradient[ iQuadPt ][ iterDim ][ jDim ] += interpCFE.phiDer ( i, jDim, iQuadPt ) *
                                                                   beta[ betaDof.localToGlobalMap (elementID, i) + iterDim * totalDof];
            }
        }
    }

}


//! Interpolation of the divergence
template<typename localVector, typename globalVector>
void interpolateDivergence (localVector& localDivergence,
                            const CurrentFE& interpCFE,
                            const DOF& betaDof,
                            const UInt& elementID,
                            const globalVector& beta)
{
    const UInt nbQuadPt (interpCFE.nbQuadPt() );
    const UInt nbFEDof (interpCFE.nbFEDof() );
    const UInt totalDof (betaDof.numTotalDof() );

    // Loop on the quadrature nodes
    for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
    {
        localDivergence[ iQuadPt ] = 0.0;

        for ( UInt i = 0; i < nbFEDof; ++i )
        {
            for (UInt jDim (0); jDim < nDimensions; ++jDim)
                localDivergence[ iQuadPt ] += interpCFE.phiDer ( i, jDim, iQuadPt ) *
                                              // here we are assuming that beta belongs to the FE space of interpCFE
                                              beta[ betaDof.localToGlobalMap (elementID, i) + jDim * totalDof ];
        }
    }

}


template<typename localVector>
void massDivW (MatrixElemental& localMass,
               const CurrentFE& massCFE,
               const Real& coefficient,
               const localVector& localValues,
               const UInt& fieldDim)
{
    const UInt nbFEDof (massCFE.nbFEDof() );
    const UInt nbQuadPt (massCFE.nbQuadPt() );
    Real localValue (0);

    // Assemble the local mass
    for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        MatrixElemental::matrix_view localView = localMass.block (iterFDim, iterFDim);

        // Loop over the basis functions
        for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
        {
            // Build the local matrix only where needed:
            // Lower triangular + diagonal parts
            for (UInt jDof (0); jDof <= iDof; ++jDof)
            {
                localValue = 0.0;

                // Loop on the quadrature nodes
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += localValues[iQuadPt]
                                  * massCFE.phi (iDof, iQuadPt)
                                  * massCFE.phi (jDof, iQuadPt)
                                  * massCFE.wDetJacobian (iQuadPt);
                }

                localValue *= coefficient;

                // Add on the local matrix
                localView (iDof, jDof) += localValue;

                if (iDof != jDof)
                {
                    localView (jDof, iDof) += localValue;
                }
            }
        }
    }
}

//! Elementary advection \beta \grad u v
template<typename localVector>
void advection (MatrixElemental& localAdv,
                const CurrentFE& advCFE,
                const Real& coefficient,
                const localVector& localValues,
                const UInt& fieldDim)
{
    const UInt nbFEDof (advCFE.nbFEDof() );
    const UInt nbQuadPt (advCFE.nbQuadPt() );
    Real localValue (0.0), advGrad (0.0);
    MatrixElemental matTmp ( nbFEDof, 1, 1 );
    matTmp.zero();
    MatrixElemental::matrix_view matTmpView = matTmp.block (0, 0);

    // Loop over the basis functions
    for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
    {
        // Build the local matrix
        for (UInt jDof (0); jDof < nbFEDof; ++jDof)
        {
            localValue = 0.0;

            // Loop on the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                advGrad = 0.;
                for (UInt iDim (0); iDim < advCFE.nbLocalCoor(); ++iDim)
                {
                    advGrad += localValues[iQuadPt][iDim]
                               * advCFE.dphi (jDof, iDim, iQuadPt);
                }

                localValue += advGrad * advCFE.phi (iDof, iQuadPt)
                              * advCFE.wDetJacobian (iQuadPt);
            }
            // Add on the local matrix
            matTmpView (iDof, jDof) = coefficient * localValue;
        }
    }
    for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        MatrixElemental::matrix_view localView = localAdv.block (iterFDim, iterFDim);

        // Copy on the components
        localView = matTmpView;
    }
}

/*
 * Added by Gwenol Grandperrin, August 2011
 */
//! Assemble the term \f$ \int_\Omega \phi_j\cdot\mathbf{u}\phi_i\f$
void advectionNewton ( Real coef, VectorElemental& vel,
                       MatrixElemental& elmat, const CurrentFE& fe,
                       int iblock, int jblock );

//! Elementary advection, term u\grad \beta v
template<typename localTensor>
void symmetrizedAdvection (MatrixElemental& localAdv,
                           const CurrentFE& advCFE,
                           const Real& coefficient,
                           const localTensor& localGradient,
                           const UInt& fieldDim)
{
    const UInt nbFEDof (advCFE.nbFEDof() );
    const UInt nbQuadPt (advCFE.nbQuadPt() );
    Real localValue (0.0);

    ASSERT (fieldDim == nDimensions, "Symmetrized operator works only with vectors of the same dimension as the space" );


    for (UInt iCoor (0); iCoor < fieldDim; ++iCoor)
    {
        for (UInt jCoor (0); jCoor < fieldDim; ++jCoor)
        {
            // Extract the view of the matrix
            MatrixElemental::matrix_view localView = localAdv.block (iCoor, jCoor);

            // Loop over the basis functions
            for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
            {
                // Build the local matrix
                for (UInt jDof (0); jDof < nbFEDof; ++jDof)
                {
                    localValue = 0.0;

                    // Loop on the quadrature nodes
                    for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                    {
                        localValue +=
                            localGradient[iQuadPt][iCoor][jCoor]
                            * advCFE.phi (jDof, iQuadPt)
                            * advCFE.phi (iDof, iQuadPt)
                            * advCFE.wDetJacobian (iQuadPt);

                    }

                    // Add on the local matrix
                    localView (iDof, jDof) = coefficient * localValue;
                }
            }
        }
    }
}

void grad (MatrixElemental& localGrad,
           const CurrentFE& uCFE,
           const CurrentFE& pCFE,
           const UInt& fieldDim);

void divergence (MatrixElemental& localDivergence,
                 const CurrentFE& uCFE,
                 const CurrentFE& pCFE,
                 const UInt& fieldDim,
                 const Real& coefficient = 1.0);

void stiffStrain (MatrixElemental& localStiff,
                  const CurrentFE& stiffCFE,
                  const Real& coefficient,
                  const UInt& fieldDim);

void bodyForces (VectorElemental& localForce,
                 const CurrentFE& massRhsCFE,
                 const function_Type& fun,
                 const Real& t,
                 const UInt& fieldDim);
}

//----------------------------------------------------------------------
//
//!@name               Operators for classical finite elements
//@{
//----------------------------------------------------------------------
//!coef(t,x,y,z,u)

void mass ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock = 0, int jblock = 0 );
void mass ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock, int jblock, UInt nb );
//! Mass term with coefficients given for each quadrature point
void mass ( const std::vector<Real>& qpt_coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock, int jblock, UInt nb );

void stiff ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock = 0, int jblock = 0 );
void stiff ( Real coef, Real ( *fct ) ( Real, Real, Real ), MatrixElemental& elmat,
             const CurrentFE& fe, int iblock, int jblock );
void stiff ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock, int jblock, int nb );
//! Stiff term with coefficient given for each quadrature node
void stiff ( const std::vector<Real>& qpt_coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock, int jblock, int nb );
void stiff_curl ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
                  int iblock, int jblock, int nb );





//! \f$ coef \cdot ( e(u) , e(v) )\f$
void stiff_strain ( Real coef, MatrixElemental& elmat, const CurrentFE& fe );

//! \f$ coef \cdot ( div u , div v )\f$
void stiff_div ( Real coef, MatrixElemental& elmat, const CurrentFE& fe );

//! \f$ coef \cdot ( [\nabla u^k]^T \nabla u : \nabla v  )\f$
void stiff_dergradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

//! \f$ coef \cdot ( [\nabla u]^T \nabla u^k + [\nabla u^k]^T \nabla u : \nabla v  )\f$ for Newton on St-Venant
void stiff_dergrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

//! \f$ coef \cdot ( trace { [\nabla u^k]^T \nabla u }, \nabla\cdot  v  ) \f$ for Newton on St-Venant
void stiff_derdiv ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// -----------added Rita 2008   for non linear St-Venant----------------------------------------------------------

// coef * ( (\div u_k) \grad u : \grad v  )--------------------------------------------------------------------controllato!!!
void stiff_divgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( (\div u) \grad u_k : \grad v  )
// part of the jacobian of stiff_divgrad
void stiff_divgrad_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad u_k : \grad u_k) * ( \grad u : \grad v  )---------------------------------------------controllato!!!
void stiff_gradgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad u_k : \grad u) *( \grad u_k : \grad v  )
// part of the jacobian stiff_gradgrad
void stiff_gradgrad_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad u^k \grad u : \grad v  )------------------------------------------------------------------controllato!!!
void stiff_dergrad_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad \delta u \grad u^k : \grad v  )
// part of the jacobian of stiff_dergrad_gradbis
void stiff_dergrad_gradbis_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad u^k [\grad u]^T : \grad v  )------------------------------------------------------------controllato!!!
void stiff_dergrad_gradbis_Tr ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * ( \grad \delta u [\grad u^k]^T : \grad v  )------------------------------------------------------------controllato!!!
// part of the jacobian of stiff_dergrad_gradbis_Tr
void stiff_dergrad_gradbis_Tr_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * (  \grad u^k [\grad u^k]^T \grad u : \grad v  )------------------------------------------------------------controllato!!!
void stiff_gradgradTr_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

// coef * (  \grad u^k [\grad u]^T \grad u^k : \grad v  )------------------------------------------------------------controllato!!!
// part of the jacobian of  stiff_gradgradTr_gradbis
void stiff_gradgradTr_gradbis_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

//  coef * (  \grad u [\grad u^k]^T \grad u^k : \grad v  )------------------------------------------------------------controllato!!!
// secondo part of the jacobian of stiff_gradgradTr_gradbis
void stiff_gradgradTr_gradbis_3 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

//------------------------------------------------------------------------------------------------------------------------------------------


void grad ( const int icoor, Real coef, MatrixElemental& elmat,
            const CurrentFE& fe_u, const CurrentFE& fe_p,
            int iblock = 0, int jblock = 0 );
void div ( const int icoor, Real coef, MatrixElemental& elmat,
           const CurrentFE& fe_u, const CurrentFE& fe_p,
           int iblock = 0, int jblock = 0 );
void grad_div ( Real coef_grad, Real coef_div, MatrixElemental& elmat,
                const CurrentFE& fe_u, const CurrentFE& fe_p,
                int block_pres );

template<typename UsrFct>
void advection ( Real coef, const UsrFct& beta,
                 MatrixElemental& elmat, const CurrentFE& fe, int iblock, int jblock, int nb, Real t = 0. )
{
    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    Real v, s;
    Real x, y, z;
    Matrix v_grad (ZeroMatrix (fe.nbFEDof(), fe.nbQuadPt() ) );

    //Evaluate the velocity field at the quadrature nodes
    for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
    {
        fe.coorQuadPt (x, y, z, iq);
        for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
        {
            v = beta (t, x, y, z, icoor);
            for ( UInt j = 0; j < fe.nbFEDof(); ++j)
            {
                v_grad (j, iq) += v * fe.phiDer (j, icoor, iq );
            }
        }
    }


    //Assemble the local matrix
    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
        {
            s = 0.;
            for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
            {
                s += v_grad (j, iq) * fe.phi ( i, iq ) * fe.weightDet ( iq );
            }
            mat_tmp ( i, j ) = s * coef;
        }
    }

    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( UInt i = 0; i < fe.nbDiag(); i++ )
        {
            for ( UInt j = 0; j < fe.nbDiag(); j++ )
            {
                mat_icomp ( i, j ) += mat_tmp ( i, j );
            }
        }
    }
}

void stab_stokes ( Real visc, Real coef_stab, MatrixElemental& elmat,
                   const CurrentFE& fe, int block_pres );
void advection ( Real coef, VectorElemental& vel, MatrixElemental& elmat,
                 const CurrentFE& fe, int iblock, int jblock, int nb );

void source ( Real constant, VectorElemental& elvec, const CurrentFE& fe, int iblock );

void source ( Real constant, VectorElemental& elvec, const CurrentFE& fe, Real t, int iblock );


// right-hand sides for Chorin-Teman projection scheme
void source_divuq (Real alpha, VectorElemental& uLoc,  VectorElemental& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock = 0 );
void source_gradpv (Real alpha, VectorElemental& pLoc,  VectorElemental& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock );

//! Assembly for the source term \f$ \int c v \f$ where \f$c\f$ is a given by the values in the quadrature nodes.
/*!
  This function add in the elementary vector the term \f$ \int c v \f$.
  The function \f$c\f$ is given by its values in the quadrature nodes.

  @param constant Values of the function in the quadrature nodes
  @param elvec The local vector where to add the values
  @param currentFe The currentFE associated to the cell where to assemble
  @param iblock The component of v that is concerned
*/
void source_mass (const std::vector<Real>& constant, VectorElemental& elvec, const CurrentFE& currentFe, const int& iblock);

//! Assembly for the source term \f$ \int \nabla c \cdot \nabla v \f$ where \f$c\f$ is a given by the values in the quadrature nodes.
/*!
  The function \f$\nabla c\f$ is given by its values in the quadrature nodes, coordinate after coordinate (first,
  the values for the first componant of the gradient in all the quadrature nodes, then second component,...).

  @param constant Values of the gradient in the quadrature nodes
  @param elvec The local vector where to add the values
  @param currentFe The currentFE associated to the cell where to assemble
  @param iblock The component of v that is concerned
*/
void source_stiff (const std::vector<Real>& constant, VectorElemental& elvec, const CurrentFE& currentFe, const int& iblock);

//!@}
//!@name Elementary operations for the interior penalty stabilization
//!@{
//

//! \f$ coef < \nabla p1, \nabla q2 >\f$
void ipstab_grad ( const Real coef, MatrixElemental& elmat,
                   const CurrentFE& fe1, const CurrentFE& fe2,
                   const CurrentFEManifold& bdfe, int iblock = 0, int jblock = 0 );

//! \f$ coef < \nabla u1, \nabla v2 >\f$
void ipstab_grad ( const Real coef, MatrixElemental& elmat,
                   const CurrentFE& fe1, const CurrentFE& fe2,
                   const CurrentFEManifold& bdfe, int iblock, int jblock, int nb );

//! \f$ coef < \nabla\cdot  u1, \nabla\cdot  v2 >\f$
void ipstab_div ( const Real coef, MatrixElemental& elmat,
                  const CurrentFE& fe1, const CurrentFE& fe2,
                  const CurrentFEManifold& bdfe, int iblock = 0, int jblock = 0 );
//! \f$ coef < \beta1 . \nabla u1, \beta2 . \nabla v2 >\f$
void ipstab_bgrad ( const Real coef, MatrixElemental& elmat,
                    const CurrentFE& fe1, const CurrentFE& fe2,
                    const VectorElemental& beta, const CurrentFEManifold& bdfe,
                    int iblock, int jblock, int nb );
//! \f$ coef < |\beta . n|^2 / |\beta| \nabla p1, \nabla q2 >\f$
void ipstab_bagrad ( const Real coef, MatrixElemental& elmat,
                     const CurrentFE& fe1, const CurrentFE& fe2,
                     const VectorElemental& beta, const CurrentFEManifold& bdfe,
                     int iblock = 0, int jblock = 0 );

//!\f$ coef < |\beta\cdot n| \nabla p1, \nabla q2 >\f$
//! p1 lives in fe1
//! q2 lives in fe2
//! beta lives in fe3

void ipstab_bagrad ( const Real           coef,
                     MatrixElemental&             elmat,
                     const CurrentFE&     fe1,
                     const CurrentFE&     fe2,
                     const CurrentFE&     fe3,
                     const VectorElemental&       beta,
                     const CurrentFEManifold&   bdfe,
                     int iblock = 0, int jblock = 0 );

//!@}
///////////////////////////////////////


////////////////////////////////////////
//! Convective term with a local vector coefficient (useful for Navier-Stokes problem)
void grad ( const int icoor, const VectorElemental& vec_loc, MatrixElemental& elmat,
            const CurrentFE& fe1, const CurrentFE& fe2,
            int iblock, int jblock );

//! Convective term with a local vector coefficient (useful for Navier-Stokes problem+adv-diff)
void grad ( const int icoor, const VectorElemental& vec_loc, MatrixElemental& elmat,
            const CurrentFE& fe1, const CurrentFE& fe2, const CurrentFE& fe3,
            int iblock = 0, int jblock = 0 );

//! Conective term with a local vector given by quadrature node
/*!
  To use this function, we must ensure that the velocity is stored in a good way in the std::vector:
  If there are nQ quadrature nodes, the i-th component (starting from 0) of the velocity in the iq-th quadrature node
  (also starting from 0) has to be stored in the ( i*nQ + iq)-th element of the std::vector.
*/
void grad ( const int& icoor, const std::vector<Real>& localVector, MatrixElemental& elmat,
            const CurrentFE& currentFE1, const CurrentFE& currentFE2,
            const int& iblock = 0, const int& jblock = 0);

//! Convective term with a local vector coefficient for Navier-Stokes problem in Skew-Symmetric form
void grad_ss ( const int icoor, const VectorElemental& vec_loc, MatrixElemental& elmat,
               const CurrentFE& fe1, const CurrentFE& fe2,
               int iblock = 0, int jblock = 0 );

//! StreamLine Diffusion
void stiff_sd ( Real coef, const VectorElemental& vec_loc, MatrixElemental& elmat, const CurrentFE& fe,
                const CurrentFE& fe2, int iblock = 0, int jblock = 0, int nb = 1 );

/////////////////////////////////////////////
//
//! source  \f$ \int fct \phi_i \f$
//
template <typename UsrFct>
void source ( const UsrFct& fct, VectorElemental& elvec, const CurrentFE& fe, int iblock = 0 )
{
    UInt i, ig;
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real s;
    for ( i = 0; i < fe.nbFEDof(); i++ )
    {
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( i, ig ) * fct ( fe.quadPt ( ig, 0 ), fe.quadPt ( ig, 1 ), fe.quadPt ( ig, 2 ),
                                          iblock) * fe.weightDet ( ig );
        }
        vec ( i ) += s;
    }
}

/////////////////////////////////////////////
//
//! source  \f$ \int fct(t) \phi_i\f$
//
template <typename UsrFct>
void source ( const UsrFct& fct, VectorElemental& elvec, const CurrentFE& fe, Real t, int iblock = 0 )
{
    UInt i, ig;
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real s;
    for ( i = 0; i < fe.nbFEDof(); i++ )
    {
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( i, ig ) * fct (t, fe.quadPt ( ig, 0 ), fe.quadPt ( ig, 1 ), fe.quadPt ( ig, 2 ),
                                         iblock ) * fe.weightDet ( ig );
        }
        vec ( i ) += s;
    }
}


void source ( Real coef, VectorElemental& f, VectorElemental& elvec, const CurrentFE& fe,
              int fblock = 0, int eblock = 0 );

void source_fhn ( Real coef_f, Real coef_a, VectorElemental& u, VectorElemental& elvec, const CurrentFE& fe,
                  int fblock = 0, int eblock = 0 );

//! \f$( beta\cdot\nabla u^k, v  )\f$
void source_advection ( const Real& coefficient, const VectorElemental& beta_loc, const VectorElemental& uk_loc,
                        VectorElemental& elvec, const CurrentFE& fe );

//!@name Shape derivative terms for Newton FSI
//!@{
//
//! \f$ coef \cdot ( \nabla (-w^k):[I\nabla\cdot  d - (\nabla d)^T] u^k + convect^T[I\nabla\cdot  d - (\nabla d)^T] (\nabla u^k)^T , v  ) \f$ for Newton FSI
//
//! Remark: convect = \f$u^n-w^k\f$
//
void source_mass1 ( Real coef,
                    const VectorElemental& uk_loc,
                    const VectorElemental& wk_loc,
                    const VectorElemental& convect_loc,
                    const VectorElemental& d_loc,
                    VectorElemental& elvec,
                    const CurrentFE& fe );

//
//! \f$coef \cdot ( \nabla u^k dw, v  )\f$ for Newton FSI
//
//
void source_mass2 ( Real coef, const VectorElemental& uk_loc, const VectorElemental& dw_loc,
                    VectorElemental& elvec, const CurrentFE& fe );



void  source_mass3 ( Real coef, const VectorElemental& un_loc, const VectorElemental& uk_loc, const VectorElemental& d_loc,
                     VectorElemental& elvec, const CurrentFE& fe );





//
//! \f$coef \cdot ( [-p^k I + 2\mu e(u^k)] [I\nabla\cdot d - (\nabla d)^T] , \nabla v  )\f$ for Newton FSI
//
void source_stress ( Real coef, Real mu, const VectorElemental& uk_loc, const VectorElemental& pk_loc,
                     const VectorElemental& d_loc, VectorElemental& elvec, const CurrentFE& fe_u,
                     const CurrentFE& fe_p );

//
//! \f$+ \mu ( \nabla u^k \nabla d + [\nabla d]^T[\nabla u^k]^T : \nabla v )\f$
//
void source_stress2 ( Real coef, const VectorElemental& uk_loc, const VectorElemental& d_loc, VectorElemental& elvec, const CurrentFE& fe_u );



/**
 *Shape terms for the CE system in FSI Newton. It is the sum of the terms \c source_mass1\c, \c source_stress\c and \c source_stress2\c.
 *the difference between this method and the previous ones is that here the shape terms are assembled in a matrix, instead
 *of a vector. This implies some extra loop and the explicit construction of the tensor \f$[I\nabla\cdot  d - (\nabla d)^T]\f$
 *instead of it's multiplication times a vector. However the computation of these terms need to be done once per Newton
 *iterations, instead of once per Jacobian-vector multiplication.

 *Note that the term \c source_mass2\c is not considered if the fluid domain velocity \f w\f is trated explicitly.
 *This method is currently tested only for the P1-P1 stabilized space discretization.
 */
void shape_terms (
    //const VectorElemental& d_loc,
    Real coef,
    Real mu,
    const VectorElemental& un_loc,
    const VectorElemental& uk_loc,
    const VectorElemental& wk_loc,
    const VectorElemental& convect_loc,
    const VectorElemental& pk_loc,
    MatrixElemental& elmat,
    const CurrentFE& fe,
    const CurrentFE& fe_p,
    ID mmDim,
    MatrixElemental& /*elmatP*/,
    int iblock = 0,
    bool wImplicit = false,
    Real alpha = 0.,
    boost::shared_ptr<MatrixElemental> elmat_convect = boost::shared_ptr<MatrixElemental>()
);

//
//! \f$coef * (  (\nabla u^k):[I\nabla\cdot  d - (\nabla d)^T] , q  )\f$ for Newton FSI
//
void source_press ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat,
                    const CurrentFE& fe_u, const CurrentFE& fe_p, ID mmDim, int iblock = 0 );



//
//! \f$coef * (  (\nabla u^k):[I\nabla\cdot  d - (\nabla d)^T] , q  )\f$ for Newton FSI
//
void source_press ( Real coef, const VectorElemental& uk_loc, const VectorElemental& d_loc, VectorElemental& elvec,
                    const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock = 0 );


void source_press2 ( Real coef, const VectorElemental& p_loc, const VectorElemental& d_loc, VectorElemental& elvec,
                     const CurrentFE& fe, int iblock = 0 );
//@}
//!@name Mass matrix
//!@{
//-------------Mass matrix---------------------------------------
/*!
  Weighted Mass matrix with a permeability tensor which is a constant scalar matrix
  (i.e. \f$K^{-1} = coef \cdot Id\f$, coef being the inverse of the permeability).
*/

void mass_divw ( Real coef, const VectorElemental& w_loc, MatrixElemental& elmat, const CurrentFE& fe,
                 int iblock, int jblock, UInt nb );

//! Idem \c mass_divw \c, but with coefficient given by quadrature node
void mass_divw (const std::vector<Real>& coef, const VectorElemental& w_loc, MatrixElemental& elmat, const CurrentFE& fe,
                int iblock, int jblock, UInt nb );


void mass_gradu ( Real coef, const VectorElemental& u0_loc, MatrixElemental& elmat, const CurrentFE& fe );



//!@}

//!@name Cholesky
//!@{
//-------------Cholesky---------------------------------------
/*!
  Cholesky decomposition and solution for a KNM matrix.
*/
void choldc ( KNM<Real>& a, KN<Real>& p );
void cholsl ( KNM<Real>& a, KN<Real>& p, KN<Real>& b, KN<Real>& x );
//!@}


//!@name Operators for H(div) finite elements
//!@{
//-------------Operators for H(div) finite elements------------

// Gradient matrix.
/*!
Compute the gradient of an element in \f$ H(div, K ) \f$ space, i.e. the opposite of the transpose divergence
matrix, with \f$ K \f$ the current element. In formula
\f[
\mathrm{coef}  < \nabla q, w > \equiv - \mathrm{coef} < q, \nabla \cdot w > \,,
\f]
for \f$ q \in L^2(K) \f$, \f$ w \in H(div, K) \f$ and \f$ \mathrm{coef} \f$ a real scalar.
@param coef Constant real coefficient.
@param elmat Mixed element matrix.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param primalFE Current primal finite element in \f$ L^2(K) \f$.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
*/
void grad_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE,
                 int iblock = 0, int jblock = 0 );

// Divergence matrix.
/*!
Compute the divergence of an element in \f$ H(div,K ) \f$, with \f$ K \f$ the current element.
In formula
\f[
\mathrm{coef} < q, \nabla \cdot w > \,,
\f]
for \f$ q \in L^2(K) \f$, \f$ w \in H(div, K) \f$ and \f$ \mathrm{coef} \f$ a real scalar.
@param coef Constant real coefficient.
@param elmat Mixed element matrix.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param primalFE Current primal finite element in \f$ L^2(K) \f$.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
*/
void div_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE,
                int iblock = 0, int jblock = 0 );

// Hybrid variable times dual dot product outward unit normal.
/*!
Compute the product between an hybrid variable and a dual variable dot product outward unit
normal, in the current element \f$ K \f$. In formula
\f[
\mathrm{coef} < \lambda, v \cdot n >\,,
\f]
for \f$ \lambda \in H^{1/2}(\partial K) \f$, \f$ v \in H(div, K) \f$, \f$ coef \f$ a real scalar and \f$ n \f$
the normal unit vector oriented outward of the current element \f$ K \f$.
<BR>
The \f$ \lambda \f$ are the Lagrange multiplier basis functions that enforce continuity
of the normal component of the vectorial functions across two neighbouring elements.
They can be interprated as trace of primal variable.
<BR>
See Hybridization for Mixed Hybrid Finite Element Method.
<BR>
Thanks to the Piola transform, the computation is performed on the boundary of the reference element.
But in general, the boundary of a 3D reference element is not a 2D reference element.
<BR>
Example:
REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
REFERENCE HEXA  -> 6 REFERENCE QUAD.

@param coef Constant real coefficient.
@param elmat Mixed element matrix.
@param hybridFE Reference hybrid finite element.
@param dualDotNFE Reference dual dot product outward unit normal finite element.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
@note The previous way of construction, worked only for RTO hexa, calls
\code
TP_TP_Hdiv(coef, elmat, hybridFE, iblock, jblock);
\endcode
*/
void TP_VdotN_Hdiv ( Real coef, MatrixElemental& elmat, const ReferenceFEHybrid& hybridFE,
                     const ReferenceFEHybrid& dualDotNFE, int iblock = 0, int jblock = 0 );

// Mass matrix for the hybrid variable.
/*!
Compute the mass matrix for the hybrid variable, in the current element \f$ K \f$.
In formula
\f[
\mathrm{coef} < \lambda, \mu > \,,
\f]
for \f$ \lambda, \mu \in H^{1/2}(\partial K) \f$ and \f$ coef \f$ a real scalar.
<BR>
The \f$ \lambda \f$ and \f$ \mu \f$ are the Lagrange multiplier basis functions that enforce continuity
of the normal component of the vectorial functions across two neighbouring elements.
They can be interprated as trace of primal variable.
<BR>
See Hybridization for Mixed Hybrid Finite Element Method.
<BR>
Thanks to the Piola transform, the computation is performed on the boundary of the reference element.
But in general, the boundary of a 3D Reference element is not a 2D reference element.
<BR>
Example:
REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
REFERENCE HEXA  -> 6 REFERENCE QUAD.
@param coef Constant real coefficient.
@param elmat Mixed element matrix.
@param hybridFE Reference hybrid finite element.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
@note This is an obsolete function, call TP_VdotN_Hdiv instead.
*/
void TP_TP_Hdiv ( Real coef, MatrixElemental& elmat, const ReferenceFEHybrid& hybridFE, int iblock = 0, int jblock = 0);

// Mass matrix for dual variable with constant real permeability.
/*!
Compute the mass matrix in \f$ H(div, K ) \f$ with constant real permeability, with \f$ K \f$ the
current element. In formula
\f[
\mathrm{coef} < u, w > \,,
\f]
for \f$ u, w \in H(div, K) \f$ and \f$ coef \f$ a real scalar.
<BR>
Weighted Mass matrix with a permeability tensor which is a constant scalar matrix, i.e.
\f$ K^{-1} = \mathrm{coef} \, I \f$, \f$ \mathrm{coef} \f$ being the inverse of the permeability.
\attention It is \f$ \mathrm{coef} \f$ that is used, and not its inverse.
@param coef Constant real coefficient.
@param elmat Mixed element matrix.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
@note \f$ \mathrm{coeff} \f$ is the inverse of the permeability coefficient.
*/
void mass_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, int iblock = 0, int jblock = 0 );

// Mass matrix for dual variable with constant matrix permeability.
/*!
Compute the mass matrix in \f$ H(div, K ) \f$ with constant matrix permeability, with \f$ K \f$ the
current element. In formula
\f[
< \mathrm{Invperm}\, u, w > \,,
\f]
for \f$ u, w \in H(div, K) \f$ and \f$ \mathrm{Invperm} \f$ a real constant matrix.
<BR>
Weighted mass matrix in \f$ H(div, K) \f$ with permeability matrix which is a constant
per element symmetric positive definite matrix, non diagonal a priori,
and already inverted, with Lapack LU or Choleski for instance.
@param Invperm Constant coefficient tensor, constant means constant over the current element.
@param elmat Mixed element matrix.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
@note \f$ \mathrm{Invperm} \f$ is the inverse of the permeability matrix.
*/
void mass_Hdiv ( Matrix const& Invperm, MatrixElemental& elmat, const CurrentFE& dualFE,
                 int iblock = 0, int jblock = 0 );

// Mass matrix for dual variable with real function permeability.
/*!
Compute the mass matrix in \f$ H(div, K ) \f$ with real function permeability, with \f$ K \f$ the
current element. In formula
\f[
< \mathrm{InvpermFun}\, u, w > \,,
\f]
for \f$ u, w \in H(div, K) \f$ and \f$ \mathrm{InvpermFun} \f$ a real function.
<BR>
Weighted mass matrix with a permeability that is a scalar function. The inverse function
of the permeability should be provided. We note again that it is the inverse of the
permeability that is provided directly \f$ \mathrm{InvpermFun} = K^{-1} \f$.
@param InvpermFun Scalar function inverse of the permeability.
@param elmat Mixed element matrix.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param iblock Subarray index where to store the integral just computed.
@param jblock Subarray index where to store the integral just computed.
@note \f$ \mathrm{InvpermFun} \f$ is the inverse of the permeability function.
*/
void mass_Hdiv ( Real ( *InvpermFun ) ( const Real&, const Real&, const Real& ),
                 MatrixElemental& elmat, const CurrentFE& dualFE, int iblock = 0, int jblock = 0 );

// Source vector for dual variable with constant vector source.
/*!
Compute the source vector in \f$ H(div, K ) \f$ with constant vector permeability, with \f$ K \f$ the
current element. In formula
\f[
< g, w > \,,
\f]
for \f$ w \in H(div, K) \f$ and \f$ g \f$ a constant vector.
<BR>
@param source constant vector as the source.
@param elvect element vector.
@param dualFE Current dual finite element in \f$ H(div, K) \f$.
@param iblock Subarray index where to store the integral just computed.
*/
void source_Hdiv ( const Vector& source , VectorElemental& elvec, const CurrentFE& dualFE, int iblock = 0 );

//!@}

} // namespace LifeV
#endif
