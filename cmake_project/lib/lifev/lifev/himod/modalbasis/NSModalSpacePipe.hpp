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

#ifndef __NSMODALSPACEPIPE_HPP__
#define __NSMODALSPACEPIPE_HPP__

#include <lifev/himod/basis/Basis2DAbstract.hpp>
#include <lifev/himod/modalbasis/NSModalSpaceAbstract.hpp>
#include <lifev/himod/tools/ReferenceMap.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

namespace LifeV
{

//! NSModalSpacePipe - Class for Hi-Mod modal basis for Navier-Stokes Equations
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  Class representing the modalbasis in a circle.
 *
 */
class NSModalSpacePipe: public NSModalSpaceAbstract
{

public:

    //! typedef for the educatedbasis
    typedef Basis2DAbstract                                 basis2d_type;
    //! typedef for the pointer to the educatedbasis
    typedef basis2d_type*                                   basis2d_ptrType;
    typedef ReferenceMap                                    referenceMap_type;
    typedef boost::shared_ptr<MBVector_type>                MBVector_ptrType;
    typedef VectorEpetraStructured                          vector_Type;
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
    
    // Constructor for R = R(x,t)
    NSModalSpacePipe
    ( const function_Type& Rho, const function_Type& dRho, const Real& Theta,
      const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp,
      referenceMap_type* map,
      const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
      const Real& R = 1, const Real& dR = 1 ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ),
    M_fRho( Rho ), M_fdRho( dRho ), M_Theta( Theta ),
    M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ),
    M_map( map ), M_Rho( R ), M_dRho( dR ) {};
    
    // Copy constructor
    NSModalSpacePipe( const NSModalSpacePipe& ms  ):
    NSModalSpaceAbstract( ms.mx(), ms.mr(), ms.mtheta(), ms.mp() ), 
    M_fRho( ms.fRho() ), M_fdRho( ms.fdRho() ), M_Theta( ms.Theta() ),
    M_quadruleRho( new QuadratureRule( *(ms.qrRho()) ) ), 
    M_quadruleTheta( new QuadratureRule( *(ms.qrTheta()) ) ),
    M_map( ms.map() ), M_Rho( ms.Rho() ), M_dRho( ms.dRho() ),
    M_xphirho( ms.xphirho() ), M_xphitheta( ms.xphitheta() ), 
    M_xdphirho( ms.xdphirho() ), M_xdphitheta( ms.xdphitheta() ), 
    M_rphirho( ms.rphirho() ), M_rphitheta( ms.rphitheta() ), 
    M_rdphirho( ms.rdphirho() ), M_rdphitheta( ms.rdphitheta() ), 
    M_thetaphirho( ms.thetaphirho() ), M_thetaphitheta( ms.thetaphitheta() ), 
    M_thetadphirho( ms.thetadphirho() ), M_thetadphitheta( ms.thetadphitheta() ), 
    M_pphirho( ms.pphirho() ), M_pphitheta( ms.pphitheta() ), 
    M_pdphirho( ms.pdphirho() ), M_pdphitheta( ms.pdphitheta() ),
    M_xGenbasisRhoTheta( ms.xGbRhoTheta() ),
    M_rGenbasisRhoTheta( ms.rGbRhoTheta() ),
    M_thetaGenbasisRhoTheta( ms.thetaGbRhoTheta() ),
    M_pGenbasisRhoTheta( ms.pGbRhoTheta() )
    {};

    // Constructor for patient-specific geometries
    NSModalSpacePipe
    ( const UInt& Mx, const UInt& Mr, const UInt& Mtheta, const UInt& Mp, 
      referenceMap_type* map,
      const QuadratureRule* quadrho = &quadRuleLobSeg64pt, const QuadratureRule* quadtheta = &quadRuleLobSeg64pt,
      const Real& R = 1, const Real& dR = 1 ) :
    NSModalSpaceAbstract( Mx, Mr, Mtheta, Mp ),
    M_Theta( 2*M_PI ),
    M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta ),
    M_map( map ), M_Rho( R ), M_dRho( dR ) 
    {
        M_fRho = M_map->fR();
        M_fdRho = M_map->fDxR();
    };


    //! Destroys the bidimensional basis generator
    virtual ~NSModalSpacePipe()
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
                        const std::string& pPoly );

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
    
    // ---------------         Section lumping         -----------------------
    virtual void compute_r00xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha, vector_Type& R00xx ) const;
    virtual void compute_r01xx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01xx ) const;
    virtual void compute_r10xx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10xx ) const;
    virtual void compute_r11xx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R11xx ) const;
    
    virtual void compute_r00rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha, vector_Type& R00rr ) const;
    virtual void compute_r01rr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01rr ) const;
    virtual void compute_r10rr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10rr ) const;
    virtual void compute_r11rr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R11rr ) const;
    
    virtual void compute_r00tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha, vector_Type& R00tt ) const;
    virtual void compute_r01tt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01tt ) const;
    virtual void compute_r10tt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10tt ) const;
    virtual void compute_r11tt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R11tt ) const;
    
    virtual void compute_r00xr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00xr ) const;
    virtual void compute_r01xr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01xr ) const;
        
    virtual void compute_r00xt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00xt ) const;
    virtual void compute_r01xt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01xt ) const;
    
    virtual void compute_r00rx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00rx ) const;
    virtual void compute_r10rx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10rx ) const;
        
    virtual void compute_r00rt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00rt ) const;
    virtual void compute_r01rt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01rt ) const;
    virtual void compute_r10rt( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10rt ) const;
    
    virtual void compute_r00tx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00tx ) const;
    virtual void compute_r10tx( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10tx ) const;
    
    virtual void compute_r00tr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R00tr ) const;
    virtual void compute_r01tr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R01tr ) const;
    virtual void compute_r10tr( const UInt& k, const UInt& j, const Real& nu, vector_Type& R10tr ) const;
    
    virtual void compute_r00xp( const UInt& k, const UInt& j, vector_Type& R00xp ) const;
    virtual void compute_r10xp( const UInt& k, const UInt& j, vector_Type& R10xp ) const;
    
    virtual void compute_r00rp( const UInt& k, const UInt& j, vector_Type& R00rp ) const;

    virtual void compute_r00tp( const UInt& k, const UInt& j, vector_Type& R00tp ) const;

    virtual void compute_r00px( const UInt& k, const UInt& j, vector_Type& R00px ) const;
    virtual void compute_r01px( const UInt& k, const UInt& j, vector_Type& R01px ) const;
    
    virtual void compute_r00pr( const UInt& k, const UInt& j, vector_Type& R00pr ) const;

    virtual void compute_r00pt( const UInt& k, const UInt& j, vector_Type& R00pt ) const;

    // Non linear term for NS
    virtual void compute_b00rr( const UInt& k, const UInt& j, vector_Type& b00rr, const vector_Type& adv ) const;
    virtual void compute_b10rr( const UInt& k, const UInt& j, vector_Type& b10rr, const vector_Type& adv ) const;
    virtual void compute_b00tr( const UInt& k, const UInt& j, vector_Type& b00tr, const vector_Type& adv ) const;
    virtual void compute_b10tt( const UInt& k, const UInt& j, vector_Type& b10tt, const vector_Type& adv ) const;
    virtual void compute_b00rt( const UInt& k, const UInt& j, vector_Type& b00rt, const vector_Type& adv ) const;
    virtual void compute_b00tt( const UInt& k, const UInt& j, vector_Type& b00tt, const vector_Type& adv ) const;
    virtual void compute_b10xx( const UInt& k, const UInt& j, vector_Type& b10xx, const vector_Type& adv ) const;
    virtual void compute_b00xx( const UInt& k, const UInt& j, vector_Type& b00xx, const vector_Type& adv ) const;
    
    // Int2d( phi_k r dr dtheta )
    virtual    void compute_Phix( const UInt& k, vector_Type& Phix ) const;
    virtual    void compute_Phir( const UInt& k, vector_Type& Phir ) const;
    virtual    void compute_Phitheta( const UInt& k, vector_Type& Phitheta ) const;
        
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
    virtual MBMatrix_type xphirho() const
    {
        return M_xphirho;
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real rphirho (const UInt& j, const UInt& n) const
    {
        return M_rphirho[j][n];
    }
    virtual MBMatrix_type rphirho() const
    {
        return M_rphirho;
    }
    
    //! return the value of the sub-base in the node n at the mode j
    virtual Real thetaphirho (const UInt& j, const UInt& n) const
    {
        return M_thetaphirho[j][n];
    }
    virtual MBMatrix_type thetaphirho() const
    {
        return M_thetaphirho;
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real xdphirho (const UInt& j, const UInt& n) const
    {
        return M_xdphirho[j][n];
    }
    virtual MBMatrix_type xdphirho() const
    {
        return M_xdphirho;
    }
    
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real rdphirho (const UInt& j, const UInt& n) const
    {
        return M_rdphirho[j][n];
    }
    virtual MBMatrix_type rdphirho() const
    {
        return M_rdphirho;
    }
    
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real thetadphirho (const UInt& j, const UInt& n) const
    {
        return M_thetadphirho[j][n];
    }
    virtual MBMatrix_type thetadphirho() const
    {
        return M_thetadphirho;
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real xphitheta (const UInt& j, const UInt& n) const
    {
        return M_xphitheta[j][n];
    }
    virtual MBMatrix_type xphitheta() const
    {
        return M_xphitheta;
    }
    
    //! return the value of the sub-base in the node n at the mode j
    virtual Real rphitheta (const UInt& j, const UInt& n) const
    {
        return M_rphitheta[j][n];
    }
    virtual MBMatrix_type rphitheta() const
    {
        return M_rphitheta;
    }
    
    //! return the value of the sub-base in the node n at the mode j
    virtual Real thetaphitheta (const UInt& j, const UInt& n) const
    {
        return M_thetaphitheta[j][n];
    }
    virtual MBMatrix_type thetaphitheta() const
    {
        return M_thetaphitheta;
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real xdphitheta (const UInt& j, const UInt& n) const
    {
        return M_xdphitheta[j][n];
    }
    virtual MBMatrix_type xdphitheta() const
    {
        return M_xdphitheta;
    }
    
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real rdphitheta (const UInt& j, const UInt& n) const
    {
        return M_rdphitheta[j][n];
    }
    virtual MBMatrix_type rdphitheta() const
    {
        return M_rdphitheta;
    }
    
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real thetadphitheta (const UInt& j, const UInt& n) const
    {
        return M_thetadphitheta[j][n];
    }
    virtual MBMatrix_type thetadphitheta() const
    {
        return M_thetadphitheta;
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pphirho (const UInt& j, const UInt& n) const
    {
        return M_pphirho[j][n];
    }
    virtual MBMatrix_type pphirho() const
    {
        return M_pphirho;
    }
 
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pdphirho( const UInt& j, const UInt& n) const
    {
        return M_pdphirho[j][n];
    }
    virtual MBMatrix_type pdphirho() const
    {
        return M_pdphirho;
    }
 
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pphitheta (const UInt& j, const UInt& n) const
    {
        return M_pphitheta[j][n];
    }
    virtual MBMatrix_type pphitheta() const
    {
        return M_pphitheta;
    }
 
    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real pdphitheta (const UInt& j, const UInt& n) const
    {
        return M_pdphitheta[j][n];
    }
    virtual MBMatrix_type pdphitheta() const
    {
        return M_pdphitheta;
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

    Basis2DAbstract::EigenContainer            M_xEigenvalues;
    Basis2DAbstract::EigenContainer            M_rEigenvalues;
    Basis2DAbstract::EigenContainer            M_thetaEigenvalues;
    Basis2DAbstract::EigenContainer            M_pEigenvalues;

    const QuadratureRule*    M_quadruleRho;
    const QuadratureRule*    M_quadruleTheta;

    basis2d_ptrType            M_xGenbasisRhoTheta;
    basis2d_ptrType            M_rGenbasisRhoTheta;
    basis2d_ptrType            M_thetaGenbasisRhoTheta;
    basis2d_ptrType            M_pGenbasisRhoTheta;

    const Real                    M_Rho;
    const Real                    M_dRho;
    const Real                    M_Theta;

    function_Type        M_fRho;
    function_Type        M_fdRho;

    MBMatrix_type            M_xphirho;
    MBMatrix_type            M_xphitheta;
    MBMatrix_type            M_xdphirho;
    MBMatrix_type            M_xdphitheta;
    
    MBMatrix_type            M_rphirho;
    MBMatrix_type            M_rphitheta;
    MBMatrix_type            M_rdphirho;
    MBMatrix_type            M_rdphitheta;
    
    MBMatrix_type            M_thetaphirho;
    MBMatrix_type            M_thetaphitheta;
    MBMatrix_type            M_thetadphirho;
    MBMatrix_type            M_thetadphitheta;
    
    MBMatrix_type            M_pphirho;
    MBMatrix_type            M_pphitheta;
    MBMatrix_type            M_pdphirho;
    MBMatrix_type            M_pdphitheta;
    
    referenceMap_type*        M_map;

};
}// namespace LifeV
#endif // __MODALSPACECIRCULAR_HPP_
