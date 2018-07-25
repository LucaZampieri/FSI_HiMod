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
 *   @file ChebyshevQuadraticShift.hpp
     @brief This file contains the interface of a two-dimensional uneducated basis
             with quadratically shifted Chebyshev polynomials T_j( 2*r^2 - 1 ) * sin( k * theta ), cos( k * theta ),
             multiplied by r if k is odd.

     @date 04/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
 */
 
#ifndef __CHEBYSHEVQUADRATICSHIFT_HPP__
#define __CHEBYSHEVQUADRATICSHIFT_HPP__

#include <lifev/himod/basis/UneducatedBasisAbstract.hpp>

namespace LifeV
{

class ChebyshevQuadraticShift : public UneducatedBasisAbstract
{
public:
    
    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
    The only available constructor: this class should be used only through the factory so it needs only to be
    default-constructable.
     */
    ChebyshevQuadraticShift() : UneducatedBasisAbstract() {};
    //@}
    
    //! @name Methods
    //@{
    
    virtual void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta );
    //@}
    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const;
    
/*    METHODS INHERITED FROM UNEDUCATEDBASISABSTRACT

    SETMETHODS
    setRho
    setTheta
    setNumberModes
    evaluateBasis
    evaluateSinglePoint
    
    GETMETHODS
    eigenvalues
*/

protected:
/*     PROTECTED MEMBERS INHERITED FROM UNEDUCATEDBASISABSTRACT
    Real M_Rho;
    Real M_Theta;
    UInt M_icurrent;
    UInt M_mtot;
    MBVector_type M_eigenvalues;
    Real M_normTheta;
    MBVector_type M_normRJ;
    void addRow( spy_type& matrix );
    void addColumn( spy_type& matrix );
    void deleteBesselRoots( MBMatrix_type& edges, MBMatrix_type& dedges );
*/
    // Chebyshev polynomial of order j
    Real T( const UInt& j, const Real& x ) const
    {
        return cos( j * acos( x ) );
    };

    Real dT( const UInt& j, const Real& x ) const
    {
        return j * sin( j * acos( x ) ) / sqrt( 1 - x * x );
    };

    Real radialFunction( const UInt& pr, const UInt& pt, const Real& rStar, const QuadratureRule* quadrulerho ) const;
    Real angularFunction( const UInt& pt, const Real& tStar ) const;

    Real xi( const UInt& kr, const UInt& kt, const Real& rStar, 
             const QuadratureRule* quadrulerho, const MBMatrix_type& phirho ) const;

};

// ===================================================
// Macros
// ===================================================

//! define the EducatedBasisD
inline
UneducatedBasisAbstract* createChebyshevQuadraticShift()
{
    return new ChebyshevQuadraticShift();
}

namespace
{
static bool registerChebyshevQuadraticShift = UneducatedBasisFactory::instance().registerProduct( "ChebyshevQuadraticShift",  &createChebyshevQuadraticShift );
}

} //End LifeV namespace
#endif
