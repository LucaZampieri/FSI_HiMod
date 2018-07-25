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
 *   @file EducatedBasisR.hpp
     @brief This file contains the interface of a two-dimensional educated basis
             with Robin boundary conditions.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
 */
 
#ifndef __EDUCATEDBASISR_HPP__
#define __EDUCATEDBASISR_HPP__

#include <lifev/himod/basis/EducatedBasisAbstract.hpp>

namespace LifeV
{

class EducatedBasisR : public EducatedBasisAbstract
{
public:

    //!@name Public Types
    //@{
    typedef EducatedBasisFunctorR            functor_type;

    typedef boost::shared_ptr<functor_type> functor_ptrType;
    //@}

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
    The only available constructor: this class should be used only through the factory so it needs only to be
    default-constructable.
    */
    EducatedBasisR() : EducatedBasisAbstract() {};
    //@}

    //! @name Methods
    //@{
    
       //! computeEigenvalues
    //! This method computes the eigenvalues of the bidimensional problem for Robin boundary conditions.
    virtual void computeEigenvalues();

    //@}

    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const {};
    //!@name SetMethods
    //@{
    //! setMu
    // override of the omonimous function in Basis2DAbstract
    void setMu( const Real& mu )
    {
        M_mu = mu;
    };
    
    //! setChi
    // override of the omonimous function in Basis2DAbstract    
    void setChi( const Real& chi )
    {
        M_chi = chi;
    };
    //@}

    // getMethod inherited from Basis2DAbstract - override
    Real chi() const
    {
        return M_chi;
    }

/*    METHODS INHERITED FROM EDUCATEDBASISABSTRACT

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

    Real            M_mu;
    Real            M_chi;
    functor_ptrType M_ptrFunctor;
    
    Real Bisection( functor_ptrType f, Real a, Real b, const Real tol, const int maxIt );
    
/*     PROTECTED MEMBERS INHERITED FROM EDUCATEDBASISABSTRACT
    Real             M_Rho;
    Real             M_Theta;
    UInt             M_mtot;
    EigenContainer     M_eigenvalues;
    Real             M_normTheta;
    MBVector_type     M_normRJ;
    void addRow( spy_type& matrix );
    void addColumn( spy_type& matrix );
    void deleteBesselRoots( MBMatrix_type& edges, MBMatrix_type& dedges );
*/

};

// ===================================================
// Macros
// ===================================================

//! define the EducatedBasisR
inline
EducatedBasisAbstract* createR()
{
    return new EducatedBasisR();
}

namespace
{
static bool registerR = EducatedBasisFactory::instance().registerProduct( "rob",  &createR );
}

} //End LifeV namespace
#endif //__EDUCATEDBASISR_HPP__
