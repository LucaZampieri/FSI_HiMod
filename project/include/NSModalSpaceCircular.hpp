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
    on the transversal circular fiber for the Navier-Stokes Equations

     @date 03/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>

 */

#ifndef __NSMODALSPACECIRCULAR_HPP__
#define __NSMODALSPACECIRCULAR_HPP__

#include "QuadratureRule.hpp"
#include <lifev/himod/basis/Basis2DAbstract.hpp>
#include <lifev/himod/modalbasis/NSModalSpaceAbstract.hpp>
#include "ReferenceMap.hpp"

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

namespace LifeV
{

//! NSModalSpace - Class for Hi-Mod modal basis for Navier-Stokes Equations
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  Class representing the modalbasis in a circle.
 *
 */
class NSModalSpaceCircular: public NSModalSpaceAbstract
{

public:

    //! typedef for the educatedbasis
    typedef Basis2DAbstract                         		basis2d_type;

    //! typedef for the pointer to the educatedbasis
    typedef basis2d_type*                                   basis2d_ptrType;

    typedef ReferenceMap									referenceMap_type;

    typedef boost::shared_ptr<MBVector_type>				MBVector_ptrType;

    typedef VectorEpetraStructured              			vector_Type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //! @name Constructor & Destructor
    //@{
    /*!
      The only constructor available
      @param Rho radius of the transversal section of the domain
      @param Theta angle on the transversal section of the domain
      @param M dimension of the modal basis
      @param quadrho,quadtheta the quadrature rules
    */
    NSModalSpaceCircular
    ( const Real& Rho, const Real& dRho, const Real& Theta, const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp,
    referenceMap_type* map, const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
    const function_Type& R = fake, const function_Type& dR = fake ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ), M_Rho( Rho ), M_dRho( dRho ), M_Theta( Theta ),
    M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ), M_map( map ), M_fRho( R ), M_fdRho( dR ) {};

    NSModalSpaceCircular
    ( const function_Type& Rho, const function_Type& dRho, const Real& Theta, const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp,
    referenceMap_type* map, const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
    const Real& R = 1, const Real& dR = 1 ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ), M_fRho( Rho ), M_fdRho( dRho ), M_Theta( Theta ),
    M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ), M_map( map ), M_Rho( R ), M_dRho( dR ) {};

    //! Destroys the bidimensional basis generator
    virtual ~NSModalSpaceCircular()
    {
	delete M_xGenbasisRhoTheta;
	delete M_rGenbasisRhoTheta;
	delete M_thetaGenbasisRhoTheta;
	delete M_pGenbasisRhoTheta;
    };
    //@}

    //! @name Methods
    //@{
	static Real fake( const Real& t, const Real& x, const Real& r, const Real& th, const ID& i ) { return 0; };
    //! Add the BC on the trasversal section. Set the basis generator on the section.
    void addSliceBC ( const std::string& BCx, const Real& mux, const Real& Chix,
					    const std::string& BCr, const Real& mur, const Real& Chir,
					    const std::string& BCtheta, const Real& mutheta, const Real& Chitheta );

    void evaluateBasis();
    void evaluateBasis( const std::string& xPoly,
    					const std::string& rPoly,
    					const std::string& thetaPoly,
    					const std::string& pPoly,
    					const int trigonometricBasis );

	// ---------------         Fourier Coefficients         -----------------------
    virtual    MBVector_type xFourierCoefficients( const function_Type& g, const Real& t ) const;
    virtual    MBVector_type rFourierCoefficients( const function_Type& g, const Real& t ) const;
    virtual    MBVector_type thetaFourierCoefficients( const function_Type& g, const Real& t ) const;
    virtual    MBVector_type pFourierCoefficients( const function_Type& g, const Real& t ) const;

    virtual    MBVector_type xFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const;
    virtual    MBVector_type rFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const;
    virtual    MBVector_type thetaFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const;
    virtual    MBVector_type pFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const;

    virtual    Real xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const;
    virtual    Real rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const;
    virtual    Real thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const;
    virtual    Real pFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const;

    virtual    Real xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const;
    virtual    Real rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const;
    virtual    Real thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const;
    virtual    Real pFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const;

  	// ---------------         Section lumping         -----------------------
  	// nonLinear term for NS
    virtual		Real compute_r000xxx( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000xxr( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000xxt( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r100xxx( const UInt& k, const UInt& j, const UInt& s ) const;

    virtual		Real compute_r000rrx( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000rrr( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000rrt( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r100rrx( const UInt& k, const UInt& j, const UInt& s ) const;

    virtual		Real compute_r000ttx( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000ttr( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000ttt( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r100ttx( const UInt& k, const UInt& j, const UInt& s ) const;

    virtual		Real compute_r000rtt( const UInt& k, const UInt& j, const UInt& s ) const;
    virtual		Real compute_r000trt( const UInt& k, const UInt& j, const UInt& s ) const;

	virtual		Real compute_r000xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const;
    virtual		Real compute_r001xx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r010xx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r100xx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r101xx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r110xx( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const;
    virtual		Real compute_r001rr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r010rr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r100rr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r101rr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r110rr( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const;
    virtual		Real compute_r001tt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r010tt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r100tt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r101tt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r110tt( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000xr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r001xr( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r010xr( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000xt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r001xt( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r010xt( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000rx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r100rx( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000rt( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000tx( const UInt& k, const UInt& j, const Real& nu ) const;
    virtual		Real compute_r100tx( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000tr( const UInt& k, const UInt& j, const Real& nu ) const;

    virtual		Real compute_r000xp( const UInt& k, const UInt& j ) const;
    virtual		Real compute_r100xp( const UInt& k, const UInt& j ) const;

    virtual		Real compute_r000rp( const UInt& k, const UInt& j ) const;

    virtual		Real compute_r000tp( const UInt& k, const UInt& j ) const;

    virtual		Real compute_r000px( const UInt& k, const UInt& j ) const;
    virtual		Real compute_r001px( const UInt& k, const UInt& j ) const;
    virtual		Real compute_r010px( const UInt& k, const UInt& j ) const;

    virtual		Real compute_r000pr( const UInt& k, const UInt& j ) const;

    virtual		Real compute_r000pt( const UInt& k, const UInt& j ) const;
    virtual		Real compute_r000pp( const UInt& k, const UInt& j ) const;

    // ---------------         Section lumping: overload for x-dependence         -----------------------
	virtual		void compute_r000xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000xx ) const;
    virtual		void compute_r001xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xx ) const;
    virtual		void compute_r010xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xx ) const;
    virtual		void compute_r100xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100xx ) const;
    virtual		void compute_r101xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101xx ) const;
    virtual		void compute_r110xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110xx ) const;

    virtual		void compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000rr ) const;
    virtual		void compute_r001rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001rr ) const;
    virtual		void compute_r010rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010rr ) const;
    virtual		void compute_r100rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100rr ) const;
    virtual		void compute_r101rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101rr ) const;
    virtual		void compute_r110rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110rr ) const;

    virtual		void compute_r000tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000tt ) const;
    virtual		void compute_r001tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001tt ) const;
    virtual		void compute_r010tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010tt ) const;
    virtual		void compute_r100tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100tt ) const;
    virtual		void compute_r101tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101tt ) const;
    virtual		void compute_r110tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110tt ) const;

    virtual		void compute_r000xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000xr ) const;
    virtual		void compute_r001xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xr ) const;
    virtual		void compute_r010xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xr ) const;

    virtual		void compute_r000xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000xt ) const;
    virtual		void compute_r001xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xt ) const;
    virtual		void compute_r010xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xt ) const;

    virtual		void compute_r000rx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000rx ) const;
    virtual		void compute_r100rx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100rx ) const;

    virtual		void compute_r000rt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000rt ) const;

    virtual		void compute_r000tx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000tx ) const;
    virtual		void compute_r100tx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100tx ) const;

    virtual		void compute_r000tr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000tr ) const;

    virtual		void compute_r000xp( const UInt& k, const UInt& j,
										vector_Type& R000xp ) const;
    virtual		void compute_r100xp( const UInt& k, const UInt& j,
										vector_Type& R100xp ) const;

    virtual		void compute_r000rp( const UInt& k, const UInt& j,
										vector_Type& R000rp ) const;

    virtual		void compute_r000tp( const UInt& k, const UInt& j,
										vector_Type& R000tp ) const;

    virtual		void compute_r000px( const UInt& k, const UInt& j,
										vector_Type& R000px ) const;
    virtual		void compute_r001px( const UInt& k, const UInt& j,
										vector_Type& R001px ) const;
    virtual		void compute_r010px( const UInt& k, const UInt& j,
										vector_Type& R010px ) const;

    virtual		void compute_r000pr( const UInt& k, const UInt& j,
										vector_Type& R000pr ) const;

    virtual		void compute_r000pt( const UInt& k, const UInt& j,
										vector_Type& R000pt ) const;


	virtual    Real compute_Phix( const UInt& k ) const;
	virtual    void compute_Phix( const UInt& k,
									vector_Type& Phix ) const;

	virtual    Real compute_Phir( const UInt& k ) const;
	virtual    void compute_Phir( const UInt& k,
									vector_Type& Phir ) const;

	virtual    Real compute_Phitheta( const UInt& k ) const;
	virtual    void compute_Phitheta( const UInt& k,
    									vector_Type& Phitheta ) const;

    //! ShowMe method
    virtual    void showMe() const;

    //@}

    //! @name Get Methods
    //@{
    //! return the radius of the section
    const Real Rho() const
    {
        return M_Rho;
    }

    const function_Type fRho() const
    {
        return M_fRho;
    }

    const function_Type fdRho() const
    {
        return M_fdRho;
    }

    const Real dRho() const
    {
        return M_dRho;
    }

    //! return the angle (2*pi)
    const Real Theta() const
    {
        return M_Theta;
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real xphirho (const UInt& j, const UInt& n) const
    {
        return M_xphirho[j][n];
    }

	//! return the value of the sub-base in the node n at the mode j
    virtual Real rphirho (const UInt& j, const UInt& n) const
    {
        return M_rphirho[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real thetaphirho (const UInt& j, const UInt& n) const
    {
        return M_thetaphirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real xdphirho (const UInt& j, const UInt& n) const
    {
        return M_xdphirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real rdphirho (const UInt& j, const UInt& n) const
    {
        return M_rdphirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real thetadphirho (const UInt& j, const UInt& n) const
    {
        return M_thetadphirho[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real xphitheta (const UInt& j, const UInt& n) const
    {
        return M_xphitheta[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real rphitheta (const UInt& j, const UInt& n) const
    {
        return M_rphitheta[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real thetaphitheta (const UInt& j, const UInt& n) const
    {
        return M_thetaphitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real xdphitheta (const UInt& j, const UInt& n) const
    {
        return M_xdphitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real rdphitheta (const UInt& j, const UInt& n) const
    {
        return M_rdphitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real thetadphitheta (const UInt& j, const UInt& n) const
    {
        return M_thetadphitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pphirho (const UInt& j, const UInt& n) const
    {
        return M_pphirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pdphirho( const UInt& j, const UInt& n) const
    {
        return M_pdphirho[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pphitheta (const UInt& j, const UInt& n) const
    {
        return M_pphitheta[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pdphitheta (const UInt& j, const UInt& n) const
    {
        return M_pdphitheta[j][n];
    }

    EigenMap2D xEigenvalues( UInt j ) const
    {
	return M_xEigenvalues[j];
    }

    EigenMap2D rEigenvalues( UInt j ) const
    {
	return M_rEigenvalues[j];
    }

    EigenMap2D thetaEigenvalues( UInt j ) const
    {
	return M_thetaEigenvalues[j];
    }

    EigenMap2D pEigenvalues( UInt j ) const
    {
	return M_pEigenvalues[j];
    }

    //! return the Quadrature rule
    const QuadratureRule*  qrRho() const
    {
        return M_quadruleRho;
    }

    //! return the Quadrature rule
    const QuadratureRule*  qrTheta() const
    {
        return M_quadruleTheta;
    }

    //! return the reference map
    referenceMap_type*  map() const
    {
        return M_map;
    }

    const basis2d_ptrType  xGbRhoTheta() const
    {
        return M_xGenbasisRhoTheta;
    }

    const basis2d_ptrType  rGbRhoTheta() const
    {
        return M_rGenbasisRhoTheta;
    }

    const basis2d_ptrType  thetaGbRhoTheta() const
    {
        return M_thetaGenbasisRhoTheta;
    }

    const basis2d_ptrType  pGbRhoTheta() const
    {
        return M_pGenbasisRhoTheta;
    }

protected:

    virtual void eigensProvider();

    Basis2DAbstract::EigenContainer			M_xEigenvalues;
    Basis2DAbstract::EigenContainer			M_rEigenvalues;
    Basis2DAbstract::EigenContainer			M_thetaEigenvalues;
    Basis2DAbstract::EigenContainer			M_pEigenvalues;

    const QuadratureRule*	M_quadruleRho;
    const QuadratureRule*	M_quadruleTheta;

    basis2d_ptrType			M_xGenbasisRhoTheta;
    basis2d_ptrType			M_rGenbasisRhoTheta;
    basis2d_ptrType			M_thetaGenbasisRhoTheta;
    basis2d_ptrType			M_pGenbasisRhoTheta;

    const Real					M_Rho;
    const Real					M_dRho;
    const Real					M_Theta;

    const function_Type		M_fRho;
    const function_Type		M_fdRho;

    MBMatrix_type			M_xphirho;
    MBMatrix_type			M_xphitheta;
    MBMatrix_type			M_xdphirho;
    MBMatrix_type			M_xdphitheta;

    MBMatrix_type			M_rphirho;
    MBMatrix_type			M_rphitheta;
    MBMatrix_type			M_rdphirho;
    MBMatrix_type			M_rdphitheta;

    MBMatrix_type			M_thetaphirho;
    MBMatrix_type			M_thetaphitheta;
    MBMatrix_type			M_thetadphirho;
    MBMatrix_type			M_thetadphitheta;

    MBMatrix_type			M_pphirho;
    MBMatrix_type			M_pphitheta;
    MBMatrix_type			M_pdphirho;
    MBMatrix_type			M_pdphitheta;

    referenceMap_type*		M_map;

// -------------------------------- New FSI ------------------------------------

public:

    NSModalSpaceCircular
    ( const Real& Rho, const Real& dRho, const Real& Theta, const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp, const UInt& udof,
    referenceMap_type* map, const QuadratureRule* quadrhowall = &quadRuleBoundary, const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
    const function_Type& R = fake, const function_Type& dR = fake ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ), M_Rho( Rho ), M_dRho( dRho ), M_Theta( Theta ), M_udof( udof ),
    M_quadruleRhoWall(quadrhowall), M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ), M_map( map ), M_fRho( R ), M_fdRho( dR ),
    M_nQuadRho( M_quadruleRho->nbQuadPt() ), M_nQuadTheta( M_quadruleTheta->nbQuadPt() ) {};

    NSModalSpaceCircular
    ( const function_Type& Rho, const function_Type& dRho, const Real& Theta, const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp, const UInt& udof,
    referenceMap_type* map, const QuadratureRule* quadrhowall = &quadRuleBoundary, const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
    const Real& R = 1, const Real& dR = 1 ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ), M_fRho( Rho ), M_fdRho( dRho ), M_Theta( Theta ), M_udof( udof ),
    M_quadruleRhoWall(quadrhowall), M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ), M_map( map ), M_Rho( R ), M_dRho( dR ),
    M_nQuadRho( M_quadruleRho->nbQuadPt() ), M_nQuadTheta( M_quadruleTheta->nbQuadPt() ) {};

    //! return the Quadrature rule
    const QuadratureRule*  qrRhoWall() const
    {
        return M_quadruleRhoWall;
    }

    //! return the value of the sub-base at the mode j
    virtual Real rphirhoWall (const UInt& j) const
    {
        return M_rphirhoWall[j][5];
    }

    //! return the value of the sub-base derivate at the mode j
    virtual Real rdphirhoWall (const UInt& j) const
    {
        return M_rdphirhoWall[j][5];
    }

    void evaluateBasisFSI()
    void evaluateBasisFSI( const std::string& xPoly, const std::string& rPoly, const std::string& thetaPoly, const std::string& pPoly, const int trigonometricBasis );

    // Only x-dependence

    virtual void compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha, const Real& rho_s, const Real& h_s, const Real& e, vector_Type& R000rr ) const;

    virtual void compute_r00x( const UInt& j, const vector_Type& f, const vector_Type& u_old, const Real& alpha, vector_Type& R00x ) const;
    virtual void compute_r00r( const UInt& j, const vector_Type& f, const vector_Type& u_old, const vector_Type& urWall_old, const vector_Type& etar_old, const Real& alpha, const Real& rho_s, const Real& h_s, const Real& e, vector_Type& R00r ) const;
    virtual void compute_r00t( const UInt& j, const vector_Type& f, const vector_Type& u_old, const Real& alpha, vector_Type& R00t ) const;

    virtual void compute_b0( const UInt& j, const Real& p1, Real& B0 ) const;
    virtual void compute_bL( const UInt& j, const Real& p2, Real& BL ) const;

  protected:

      const QuadratureRule* M_quadruleRhoWall;

      MBMatrix_type     M_rphirhoWall;
      MBMatrix_type     M_rdphirhoWall;

      UInt M_udof;
      UInt M_nQuadRho;
      UInt M_nQuadTheta;

      // Helper functions
      UInt xcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
      {
        return ( j + k*M_nQuadRho + i*M_nQuadRho*M_nQuadTheta );
      }
      UInt rcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
      {
        return ( M_udof*M_nQuadRho*M_nQuadTheta + j + k*M_nQuadRho + i*M_nQuadRho*M_nQuadTheta );
      }
      UInt thetacoord2index( const UInt& i, const UInt& j, const UInt& k ) const
      {
        return ( 2*M_udof*M_nQuadRho*M_nQuadTheta + j + k*M_nQuadRho + i*M_nQuadRho*M_nQuadTheta );
      }
      UInt pcoord2index( const UInt& i, const UInt& j, const UInt& k ) const
      {
        return ( 3*M_udof*M_nQuadRho*M_nQuadTheta + j + k*M_nQuadRho + i*M_nQuadRho*M_nQuadTheta );
      }
      UInt coord2indexWall( const UInt& i, const UInt& k ) const
      {
        return ( k + i*M_nQuadTheta );
      }

// -----------------------------------------------------------------------------

};
}// namespace LifeV
#endif // __MODALSPACECIRCULAR_HPP_
