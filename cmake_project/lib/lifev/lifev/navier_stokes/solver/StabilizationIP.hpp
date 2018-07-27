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
    @brief Interior Penality Stabilization

    @author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @contributor Umberto Villa <uvilla@emory.edu>
    @maintainer Umberto Villa <uvilla@emory.edu>

    @date 08-10-2004

    Implementation of I.P. stabilization for inf-sup incompatible finite elements for the Navier-Stokes equations.
 */

//TODO Change name to StabilizationIP

#ifndef _NSIPTERMS_HPP
#define _NSIPTERMS_HPP

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <boost/shared_ptr.hpp>

#define USE_OLD_PARAMETERS 0
#define WITH_DIVERGENCE 1

namespace LifeV
{

namespace details
{
//! StabilizationIP Class
/*!
 * @brief Interior Penality Stabilization
 * @author C. Winkelmann <christoph.winkelmann@epfl.ch>
 *
 * Implementation of I.P. stabilization for inf-sup incompatible finite elements
 * for the Navier-Stokes equations. <br>
 *
 *  This function adds the following stabilization terms to the Navier-Stokes monolithic matrix:
 *  <ol>
 *  <li> Block(1,1): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_\beta [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}] + \Sigma_{f\in\mathcal{F}}\int_{f} c_d [div \mathbf{u}] [div \mathbf{v}]@f$
 *  <li> Block(1,2): 0
 *  <li> Block(2,1): 0
 *  <li> Block(2,2): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_p [\nabla p] \cdot [\nabla q]@f$
 *  </ol>
 *
 *  where
 *  <ol>
 *  <li> @f$\mathcal{F}@f$ is the set of all the internal facets in the mesh;
 *  <li> @f$[\cdot]@f$ is the jump operator: @f$[x] = x^+ - x^-@f$.
 *  </ol>
 *  and the parameter @f$ c_\beta@f$, @f$c_d@f$, @f$c_p@f$, are defined as follows:
 *  <ol>
 *  <li> @f$\displaystyle c_\beta = \frac{\gamma_\beta h}  {\|\beta\|_\infty} @f$, being @f$h@f$ the facet measure;
 *  <li> @f$\displaystyle c_d = \gamma_d h \|\beta\|_\infty @f$;
 *  <li> @f$\displaystyle c_p = \frac{\gamma_p h}{\max(\|\beta\|_\infty, \nu/h^2)}@f$.
 *  </ol>
 *  Both high Pechlet numbers and inf-sup incompatible FEM are stabilized.
 *
 */

template<typename MeshType, typename DofType>
class StabilizationIP
{
public:

    //! @name Public Types
    //@{
    typedef boost::shared_ptr<MeshType> mesh_type; //deprecated
    typedef MeshType  mesh_Type;
    typedef DofType   dof_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef boost::shared_ptr<dof_Type>  dofPtr_Type;
    //@}

    //! @name Constructor and Destructor
    //@{

    //! Default Constructor
    StabilizationIP();

    virtual ~StabilizationIP() {};
    //@}

    //! @name Methods
    //@{
    //! compute IP stabilization terms and add them into matrix
    /*!
     *  This function adds the following stabilization terms to the Navier-Stokes monolithic matrix:
     *  <ol>
     *  <li> Block(1,1): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_\beta [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}] + \Sigma_{f\in\mathcal{F}}\int_{f} c_d [div \mathbf{u}] [div \mathbf{v}]@f$
     *  <li> Block(1,2): 0
     *  <li> Block(2,1): 0
     *  <li> Block(2,2): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_p [\nabla p] \cdot [\nabla q]@f$
     *  </ol>
     *  where
     *  <ol>
     *  <li> @f$\displaystyle c_\beta = \frac{\gamma_\beta h}  {\|\beta\|_\infty} @f$, being @f$h@f$ the facet measure;
     *  <li> @f$\displaystyle c_d = \gamma_\beta h \|\beta\|_\infty @f$;
     *  <li> @f$\displaystyle c_p = \frac{\gamma_p h}{\max(\|\beta\|_\infty, \nu/h^2)}@f$.
     *  </ol>
     *  Both high Pechlet numbers and inf-sup incompatible FEM are stabilized.
     *
     *  PREREQUISITE: The velocity and the pressure field should belong to the same finite element space
     *
     *  Parameters are the followings:
     *  @param dt      Real   timestep (INPUT)
     *  @param matrix  MatrixType where the stabilization terms are added into. (OUTPUT)
     *  @param state   VectorType velocity field for the linearization of the stabilization (INPUT)
     *  @param verbose bool   whenever of not to print on screen
     */
    template<typename MatrixType, typename VectorType>
    void apply ( MatrixType& matrix, const VectorType& state, bool verbose = true );

    //! Display class informations
    /*!
     * Write information relative to the class on output
     * @param output ostream ostream were to write (Default cout)
     */
    void showMe (std::ostream& output = std::cout) const;
    //@}

    //! @name Set Methods
    //@{
    //! Set the stabilization parameter @f$\gamma_\beta@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}]@f$
    void setGammaBeta (const Real& gammaBeta)
    {
        M_gammaBeta  = gammaBeta;
    }
    //! Set the stabilization parameter @f$\gamma_d@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [div \mathbf{u}] [div \mathbf{v}]@f$
    void setGammaDiv  (const Real& gammaDiv)
    {
        M_gammaDiv   = gammaDiv;
    }
    //! Set the stabilization parameter @f$\gamma_p@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\nabla p] \cdot [\nabla q]@f$
    void setGammaPress (const Real& gammaPress)
    {
        M_gammaPress = gammaPress;
    }
    //! Set the fluid viscosity @f$\nu@f$
    void setViscosity (const Real& viscosity)
    {
        M_viscosity = viscosity;
    }
    //! Set the mesh file
    void setMesh (const meshPtr_Type mesh)
    {
        M_mesh = mesh;
    }
    //! Set Discretization
    void setDiscretization (const dofPtr_Type& dof, const ReferenceFE& refFE, CurrentFEManifold& feBd, const QuadratureRule& quadRule);
    //! Set the fespace
    template<typename MapType>
    void setFeSpaceVelocity (FESpace<mesh_Type, MapType>& feSpaceVelocity);
    //@}
private:

    //! @name Private Types
    //@{
    //! facetToPoint(i,j) = localId of jth point on ith local facet
    typedef ID ( *FTOP ) ( ID const& localFacet, ID const& point );
    //@}

    //! @name Private Constructor
    //@{
    //! Copy Constructor
    StabilizationIP (const StabilizationIP<mesh_Type, dof_Type>& original);
    //@}

    //! @name Private Attributes
    //@{
    //! Pointer to the mesh object
    meshPtr_Type  M_mesh;
    //! reference to the DofType data structure
    dofPtr_Type   M_dof;
    //! current Fe on side 1 of the current facet
    boost::shared_ptr<CurrentFE>    M_feOnSide1;
    //! current Fe on side 2 of the current facet
    boost::shared_ptr<CurrentFE>    M_feOnSide2;
    //! current boundary FE
    CurrentFEManifold*  M_feBd;
    //! Stabilization parameter @f$\gamma_\beta@f$ for @f$\int_{facet} [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}]@f$
    Real         M_gammaBeta;
    //! Stabilization parameter @f$\gamma_d@f$ for @f$\int_{facet} [div \mathbf{u}] [div \mathbf{v}]@f$
    Real         M_gammaDiv;
    //! Stabilization parameter @f$\gamma_p@f$ for @f$\int_{facet} [\nabla p] \cdot [\nabla q]@f$
    Real         M_gammaPress;
    //! Fluid viscosity @f$\nu@f$
    Real         M_viscosity;
    //! facetToPoint(i,j) = localId of jth point on ith local facet
    FTOP         M_facetToPoint;
    //@}
}; // class StabilizationIP


//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename DofType>
StabilizationIP<MeshType, DofType>::StabilizationIP() :
    M_gammaBeta ( 0.0 ),
    M_gammaDiv  ( 0.0 ),
    M_gammaPress ( 0.0 ),
    M_viscosity ( 0.0 )
{}

//=============================================================================
// Method
//=============================================================================

template<typename MeshType, typename DofType>
template<typename MatrixType, typename VectorType>
void StabilizationIP<MeshType, DofType>::apply ( MatrixType& matrix,  const VectorType& state, const bool verbose )
{
    if ( M_gammaBeta == 0. && M_gammaDiv == 0. && M_gammaPress == 0. )
    {
        return;
    }

    LifeChronoFake chronoUpdate;
    LifeChronoFake chronoBeta;
    LifeChronoFake chronoElemComp;
    LifeChronoFake chronoAssembly1;
    LifeChronoFake chronoAssembly2;
    LifeChronoFake chronoAssembly3;
    LifeChronoFake chronoAssembly4;
    LifeChronoFake chronoAssembly5;
    LifeChronoFake chronoAssembly6;
    LifeChronoFake chronoAssembly7;
    LifeChronoFake chronoAssembly8;
    LifeChronoFake chronoAssembly9;
    LifeChrono chronoAssembly;

    ID geoDimensions = MeshType::S_geoDimensions;
    MatrixElemental elMatU ( M_feOnSide1->nbFEDof(), geoDimensions    , geoDimensions   );
    MatrixElemental elMatP ( M_feOnSide1->nbFEDof(), geoDimensions + 1, geoDimensions + 1 );


    const UInt nDof = M_dof->numTotalDof();

    Real normInf;
    state.normInf (&normInf);

    // local trace of the velocity
    VectorElemental beta ( M_feBd->nbFEDof(), geoDimensions );

    UInt myFacets (0);

    chronoAssembly.start();
    // loop on interior facets
    for ( UInt iFacet ( M_mesh->numBoundaryFacets() ); iFacet < M_mesh->numFacets();
            ++iFacet )
    {
        const UInt iElAd1 ( M_mesh->facet ( iFacet ).firstAdjacentElementIdentity()  );
        const UInt iElAd2 ( M_mesh->facet ( iFacet ).secondAdjacentElementIdentity() );

        if ( Flag::testOneSet ( M_mesh->facet ( iFacet ).flag(),
                                EntityFlags::SUBDOMAIN_INTERFACE | EntityFlags::PHYSICAL_BOUNDARY ) )
        {
            //std::cout << "iElAd1 = " << iElAd1 << "; iElAd2 = " << iElAd2 << std::endl;
            continue;
        }
        ++myFacets;

        chronoUpdate.start();
        // update current finite elements
#if WITH_DIVERGENCE
        M_feBd->update ( M_mesh->facet ( iFacet ), UPDATE_W_ROOT_DET_METRIC );
#else
        M_feBd->updateMeasNormal ( M_mesh->facet ( iFacet ) );
        KNM<Real>& normal = M_feBd->normal;
#endif
        const Real hK2 = M_feBd->measure();


        M_feOnSide1->updateFirstDeriv ( M_mesh->element ( iElAd1 ) );
        M_feOnSide2->updateFirstDeriv ( M_mesh->element ( iElAd2 ) );
        chronoUpdate.stop();

        Real bmax (0);
        if (normInf != 0.)
        {
            chronoBeta.start();
            // determine bmax = ||\beta||_{0,\infty,K}
            // first, get the local trace of the velocity into beta

            // local id of the facet in its adjacent element
            UInt iFaEl ( M_mesh->facet ( iFacet ).firstAdjacentElementPosition() );
            for ( UInt iNode ( 0 ); iNode < M_feBd->nbFEDof(); ++iNode )
            {
                UInt iloc ( M_facetToPoint ( iFaEl, iNode ) );
                for ( UInt iCoor ( 0 ); iCoor < M_feOnSide1->nbLocalCoor(); ++iCoor )
                {
                    UInt ig ( M_dof->localToGlobalMap ( iElAd1, iloc ) + iCoor * nDof );

                    if (state.blockMap().LID (ig) >= 0)
                    {
                        beta.vec() [ iCoor * M_feBd->nbFEDof() + iNode ] = state ( ig);
                    }
                }
            }

            // second, calculate its max norm
            for ( UInt l ( 0 ); l < static_cast<UInt> ( M_feOnSide1->nbLocalCoor() *M_feBd->nbFEDof() ); ++l )
            {
                if ( bmax < std::fabs ( beta.vec() [ l ] ) )
                {
                    bmax = std::fabs ( beta.vec() [ l ] );
                }
            }

            chronoBeta.stop();
        }


        // pressure stabilization
        if ( M_gammaPress != 0.0 )
        {
#if USE_OLD_PARAMETERS
            Real coeffPress ( M_gammaPress * hK2 ); // P1, P2 (code)
            //Real coeffPress = M_gammaPress * sqrt( hK2 ); // P1 p nonsmooth (code)
#else
            Real coeffPress = M_gammaPress * hK2 / // Pk (paper)
                              std::max<Real> ( bmax, M_viscosity / sqrt ( hK2 ) );
#endif

            elMatP.zero();
            chronoElemComp.start();
            // coef*\int_{facet} grad u1 . grad v1
            ipstab_grad ( coeffPress, elMatP, *M_feOnSide1, *M_feOnSide1, *M_feBd,
                          geoDimensions, geoDimensions);
            chronoElemComp.stop();
            chronoAssembly1.start();

            assembleMatrix (matrix, elMatP, *M_feOnSide1, *M_dof,
                            geoDimensions, geoDimensions, geoDimensions * nDof, geoDimensions * nDof);
            chronoAssembly1.stop();

            elMatP.zero();
            chronoElemComp.start();
            // coef*\int_{face} grad u2 . grad v2
            ipstab_grad ( coeffPress, elMatP, *M_feOnSide2, *M_feOnSide2, *M_feBd,
                          geoDimensions, geoDimensions);
            chronoElemComp.stop();
            chronoAssembly2.start();

            assembleMatrix (matrix, elMatP, *M_feOnSide2, *M_dof,
                            geoDimensions, geoDimensions, geoDimensions * nDof, geoDimensions * nDof);
            chronoAssembly2.stop();

            elMatP.zero();
            chronoElemComp.start();
            // - coef*\int_{facet} grad u1 . grad v2
            ipstab_grad (- coeffPress, elMatP, *M_feOnSide1, *M_feOnSide2, *M_feBd,
                         geoDimensions, geoDimensions);
            chronoElemComp.stop();
            chronoAssembly3.start();

            assembleMatrix (matrix, elMatP, *M_feOnSide1, *M_feOnSide2, *M_dof, *M_dof,
                            geoDimensions, geoDimensions, geoDimensions * nDof, geoDimensions * nDof);
            chronoAssembly3.stop();

            elMatP.zero();
            chronoElemComp.start();
            // - coef*\int_{facet} grad u2 . grad v1
            ipstab_grad (- coeffPress, elMatP, *M_feOnSide2, *M_feOnSide1, *M_feBd,
                         geoDimensions, geoDimensions);
            chronoElemComp.stop();
            chronoAssembly4.start();

            assembleMatrix (matrix, elMatP, *M_feOnSide2, *M_feOnSide1, *M_dof, *M_dof,
                            geoDimensions, geoDimensions, geoDimensions * nDof, geoDimensions * nDof);
            chronoAssembly4.stop();
        }

        // velocity stabilization
        if ( ( M_gammaDiv != 0 || M_gammaBeta != 0 ) && bmax > 0 )
        {
#if WITH_DIVERGENCE
#if USE_OLD_PARAMETERS
            Real coeffBeta ( M_gammaBeta * hK2 / std::max<Real> (bmax, hK2) ); // code
#else
            Real coeffBeta ( M_gammaBeta * hK2 / bmax ); // paper
#endif

            Real coeffDiv ( M_gammaDiv * hK2 * bmax ); // (code and paper)
            //Real coeffDiv ( M_gammaDiv * sqrt( hK2 ) * bmax ); // ? (code)
#else
            // determine bnmax = ||\beta \cdot n||_{0,\infty,K}
            // and       bcmax = ||\beta \cross n||_{0,\infty,K}

            chronoBeta.start();
            Real bnmax ( 0. );
            Real bcmax ( 0. );
            for ( UInt iNode (0); iNode < M_feBd->nbNode; ++iNode )
            {
                Real bn ( 0 );
                for ( UInt iCoor (0); iCoor < M_feOnSide1->nbLocalCoor(); ++iCoor )
                {
                    bn += normal (iNode, iCoor) *
                          beta.vec() [ iCoor * M_feBd->nbNode + iNode ];
                    bcmax = std::max<Real>
                            (bcmax, normal (iNode, (iCoor) % 3) *
                             beta.vec() [ (iCoor + 1) % 3 * M_feBd->nbNode + iNode ] -
                             normal (iNode, (iCoor + 1) % 3) *
                             beta.vec() [ (iCoor) % 3 * M_feBd->nbNode + iNode ]);
                }
                bnmax = std::max<Real> (bnmax, bn);
            }
            chronoBeta.stop();

            Real coeffGrad = hK2 * (M_gammaBeta * bnmax + M_gammaDiv * bcmax);
#endif
            elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // coef*\int_{facet} (\beta1 . grad u1) (\beta2 . grad v2)
            ipstab_bgrad ( coeffBeta, elMatU, *M_feOnSide1, *M_feOnSide1, beta,
                           *M_feBd, 0, 0, geoDimensions );
            // coef*\int_{facet} div u1 . div v1
            ipstab_div ( coeffDiv, elMatU, *M_feOnSide1, *M_feOnSide1, *M_feBd );
#else
            // coef*\int_{facet} grad u1 . grad v1
            ipstab_grad ( coeffGrad, elMatU, *M_feOnSide1, *M_feOnSide1, *M_feBd, 0, 0,
                          geoDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly5.start();
            for ( UInt iComp ( 0 ); iComp < geoDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp < geoDimensions; ++jComp )
                {
                    assembleMatrix ( matrix, elMatU, *M_feOnSide1, *M_dof,
                                     iComp, jComp, iComp * nDof, jComp * nDof );
                }
            chronoAssembly5.stop();

            elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // coef*\int_{facet} (\beta2 . grad u2) (\beta2 . grad v2)
            ipstab_bgrad ( coeffBeta, elMatU, *M_feOnSide2, *M_feOnSide2, beta,
                           *M_feBd, 0, 0, geoDimensions );
            // coef*\int_{facet} div u2 . div v2
            ipstab_div ( coeffDiv, elMatU, *M_feOnSide2, *M_feOnSide2, *M_feBd );
#else
            // coef*\int_{facet} grad u2 . grad v2
            ipstab_grad ( coeffGrad, elMatU, *M_feOnSide2, *M_feOnSide2, *M_feBd, 0, 0,
                          geoDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly6.start();
            for ( UInt iComp ( 0 ); iComp < geoDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp < geoDimensions; ++jComp )
                {
                    assembleMatrix ( matrix, elMatU, *M_feOnSide2, *M_dof,
                                     iComp, jComp, iComp * nDof, jComp * nDof );
                }
            chronoAssembly6.stop();

            elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // - coef*\int_{facet} (\beta1 . grad u1) (\beta2 . grad v2)
            ipstab_bgrad ( -coeffBeta, elMatU, *M_feOnSide1, *M_feOnSide2, beta,
                           *M_feBd, 0, 0, geoDimensions );
            // - coef*\int_{facet} div u1 . div v2
            ipstab_div ( -coeffDiv, elMatU, *M_feOnSide1, *M_feOnSide2, *M_feBd );
#else
            // - coef*\int_{facet} grad u1 . grad v2
            ipstab_grad ( -coeffGrad, elMatU, *M_feOnSide1, *M_feOnSide2, *M_feBd, 0, 0,
                          geoDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly7.start();
            for ( UInt iComp = 0; iComp < geoDimensions; ++iComp )
                for ( UInt jComp = 0; jComp < geoDimensions; ++jComp )
                {
                    assembleMatrix ( matrix, elMatU, *M_feOnSide1, *M_feOnSide2, *M_dof, *M_dof,
                                     iComp, jComp, iComp * nDof, jComp * nDof );
                }
            chronoAssembly7.stop();

            elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // - coef*\int_{facet} (\beta2 . grad u2) (\beta1 . grad v1)
            ipstab_bgrad ( -coeffBeta, elMatU, *M_feOnSide2, *M_feOnSide1, beta,
                           *M_feBd, 0, 0, geoDimensions );
            // - coef*\int_{facet} div u2 . div v1
            ipstab_div ( -coeffDiv, elMatU, *M_feOnSide2, *M_feOnSide1, *M_feBd );
#else
            // - coef*\int_{facet} grad u2 . grad v1
            ipstab_grad ( -coeffGrad, elMatU, *M_feOnSide2, *M_feOnSide1, *M_feBd, 0, 0,
                          geoDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly8.start();
            for ( UInt iComp ( 0 ); iComp < geoDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp < geoDimensions; ++jComp )
                {
                    assembleMatrix ( matrix, elMatU, *M_feOnSide2, *M_feOnSide1, *M_dof, *M_dof,
                                     iComp, jComp, iComp * nDof, jComp * nDof );
                }
            chronoAssembly8.stop();
        }

    } // loop on interior facets
    chronoAssembly.stop();
    if (verbose)
    {
        debugStream (7101) << "\n";
        debugStream (7101) << static_cast<UInt> (state.blockMap().Comm().MyPID() )
                           <<  "  .   Updating of element   done in "
                           << chronoUpdate.diffCumul()   << " s." << "\n";
        debugStream (7101) << "   .   Determination of beta done in "
                           << chronoBeta.diffCumul()     << " s." << "\n";
        debugStream (7101) << "   .   Element computations  done in "
                           << chronoElemComp.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 1              done in "
                           << chronoAssembly1.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 2              done in "
                           << chronoAssembly2.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 3              done in "
                           << chronoAssembly3.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 4              done in "
                           << chronoAssembly4.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 5              done in "
                           << chronoAssembly5.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 6              done in "
                           << chronoAssembly6.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 7              done in "
                           << chronoAssembly7.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   chrono 8              done in "
                           << chronoAssembly8.diffCumul() << " s." << "\n";
        debugStream (7101) << "   .   total                                   "
                           << chronoAssembly.diffCumul() << " s."
                           << " myFacets = " << myFacets << "\n";
    }

} // apply(...)

template<typename MeshType, typename DofType>
void StabilizationIP<MeshType, DofType>::showMe (std::ostream& output) const
{
    output << "StabilizationIP::showMe() " << std::endl;
    output << "Fluid Viscosity: " << M_viscosity << std::endl;
    output << "Stabilization coefficient velocity SD jumps:         " << M_gammaBeta  << std::endl;
    output << "Stabilization coefficient velocity divergence jumps: " << M_gammaDiv   << std::endl;
    output << "Stabilization coefficient pressure gradient jumps:   " << M_gammaPress << std::endl;
    M_mesh->showMe (output);
    M_dof->showMe (output);
}

//=============================================================================
// Setters method
//=============================================================================
template<typename MeshType, typename DofType>
void StabilizationIP<MeshType, DofType>::setDiscretization (const dofPtr_Type& dof, const ReferenceFE& refFE, CurrentFEManifold& feBd, const QuadratureRule& quadRule)
{
    M_dof = dof;
    M_feOnSide1.reset ( new CurrentFE (refFE, getGeometricMap (*M_mesh), quadRule) );
    M_feOnSide2.reset ( new CurrentFE (refFE, getGeometricMap (*M_mesh), quadRule) );
    M_feBd = &feBd;

    M_facetToPoint = MeshType::elementShape_Type::facetToPoint;
}

template<typename MeshType, typename DofType>
template<typename MapType>
void StabilizationIP<MeshType, DofType>::setFeSpaceVelocity (FESpace<mesh_Type, MapType>& feSpaceVelocity)
{
    setMesh (feSpaceVelocity.mesh() );
    setDiscretization (feSpaceVelocity.dofPtr(), feSpaceVelocity.refFE(),
                       feSpaceVelocity.feBd(), feSpaceVelocity.qr() );
}

} // namespace details

} // namespace LifeV

#endif /* _NSIPTERMS_HPP */
