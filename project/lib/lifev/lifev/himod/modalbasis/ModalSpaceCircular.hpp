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
    @file ModalSpace.hpp
    @brief This file contains the description of the class that
    handle all the utilities related to the modal expansion
    on the transversal circular fiber

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>

 */

#ifndef __MODALSPACECIRCULAR_HPP__
#define __MODALSPACECIRCULAR_HPP__

#include <lifev/himod/basis/Basis2DAbstract.hpp>
#include <lifev/himod/modalbasis/ModalSpaceAbstract.hpp>

namespace LifeV
{

//! ModalSpace - Class for Hi-Mod modal basis
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
 *  Class representing the modalbasis in a circle.
 *
 */
class ModalSpaceCircular: public ModalSpaceAbstract
{

public:

    //! typedef for the educatedbasis
    typedef Basis2DAbstract                         		basis2d_type;

    //! typedef for the pointer to the educatedbasis
    typedef basis2d_type*                                   basis2d_ptrType;

    //! @name Constructor & Destructor
    //@{
    /*!
      The only constructor available
      @param Rho ray of the transversal section of the domain
      @param Theta angle on the transversal section of the domain
      @param M dimension of the modal basis
      @param quadrho,quadtheta the quadrature rules
    */
    ModalSpaceCircular
    (const Real& Rho, const Real& Theta, const UInt& M, const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt) :
    ModalSpaceAbstract(M), M_quadruleRho (quadrho), M_quadruleTheta (quadtheta), M_Rho (Rho), M_Theta (Theta) {};

    //! Destroys the bidimensional basis generator
    virtual ~ModalSpaceCircular()
    {
	delete M_genbasisRhoTheta;
    };
    //@}

    //! @name Methods
    //@{

    //! Add the BC on the trasversal section. Set the basis generator on the section.
    void addSliceBC ( const std::string& BC, const Real& mu = 1., const Real& Chi = 1. );

    void evaluateBasis();

    virtual    std::vector<Real> fourierCoefficients (const function_Type& g) const;

    virtual    Real fourierCoeffPointWise (const Real& x, const function_Type& f, const UInt& k) const;

    virtual    Real compute1_PhiPhi (const UInt& j, const UInt& k) const;
    virtual    Real compute2_PhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Rho]x[0,Theta]} Drho(\psi_j) \psi_k \f$
    virtual    Real compute1_DrhoPhiPhi( const UInt& j, const UInt& k ) const;
    virtual    Real compute2_DrhoPhiPhi( const UInt& j, const UInt& k ) const;

    //! This methods compute the \f$\int_{[0,Rho]x[0,Theta]} Drho(\psi_j) Drho(\psi_k)\f$
    virtual    Real compute1_DrhoPhiDrhoPhi (const UInt& j, const UInt& k) const;
    virtual    Real compute2_DrhoPhiDrhoPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Rho]x[0,Theta]} Dtheta(\psi_j) \psi_k\f$
    virtual    Real compute1_DthetaPhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Rho]x[0,Theta]} Dtheta(\psi_j) Dtheta(\psi_k)\f$
    virtual    Real compute0_DthetaPhiDthetaPhi (const UInt& j, const UInt& k) const;

    virtual    Real compute_Phi (const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Rho]} \eta_j \eta_k\f$
    virtual Real compute_rho_PhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Theta]} \xi_j \xi_k\f$
    virtual Real compute_theta_PhiPhi (const UInt& j, const UInt& k) const;

    virtual Real compute_R11 (const UInt& j, const UInt& k,
                              const function_Type& mu, const Real& x) const;

    virtual Real compute_R10 (const UInt& j, const UInt& k,
                              const function_Type& beta, const Real& x) const;

    virtual Real compute_R00 (const UInt& j, const UInt& k,
                              const function_Type& mu, const function_Type& beta, const function_Type& sigma, const Real& x) const;

    //! ShowMe method
    virtual    void showMe() const;

    //@}

    //! @name Get Methods
    //@{
    //! return the radius of the section
    Real Rho() const
    {
        return M_Rho;
    }

    //! return the angle (2*pi)
    Real Theta() const
    {
        return M_Theta;
    }


    //! return the value of the sub-base in the node n at the mode j
    virtual Real phirho (const UInt& j, const UInt& n) const
    {
        return M_phirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real dphirho (const UInt& j, const UInt& n) const
    {
        return M_dphirho[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real phitheta (const UInt& j, const UInt& n) const
    {
        return M_phitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real dphitheta (const UInt& j, const UInt& n) const
    {
        return M_dphitheta[j][n];
    }

    EigenMap2D eigenvalues ( UInt j ) const
    {
		return M_eigenvalues[j];
    }

    //! return the Quadrature rule
    const QuadratureRule&  qrRho() const
    {
        return *M_quadruleRho;
    }

    //! return the Quadrature rule
    const QuadratureRule&  qrTheta() const
    {
        return *M_quadruleTheta;
    }

    const basis2d_ptrType  gbRhoTheta() const
    {
        return M_genbasisRhoTheta;
    }

protected:

    virtual void eigensProvider();

    EducatedBasisAbstract::EigenContainer			M_eigenvalues;

    const QuadratureRule*	M_quadruleRho;

    const QuadratureRule*	M_quadruleTheta;

    basis2d_ptrType			M_genbasisRhoTheta;

    Real					M_Rho;
    
    Real					M_Theta;

    MBMatrix_type			M_phirho;

    MBMatrix_type			M_phitheta;

    MBMatrix_type			M_dphirho;

    MBMatrix_type			M_dphitheta;    

};
}// namespace LifeV
#endif // __MODALSPACECIRCULAR_HPP_
