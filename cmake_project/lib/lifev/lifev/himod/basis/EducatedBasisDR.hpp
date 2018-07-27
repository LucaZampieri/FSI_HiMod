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
 *   @file EducatedBasisRR.hpp
     @brief This file contains the interface of an educated basis (monodimensional) with Robin conditions on both side.
            \f[ \mu \partial_{\nu} u + \sigma u = 0\f]

     @date 06/2013
     @author M. Aletti <teo.aletti@gmail.com>
     @author A. Bortolossi


 */
#ifndef __EDUCATEDBASISDR_HPP__
#define __EDUCATEDBASISDR_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/himod/basis/Basis1DAbstract.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

class EducatedBasisDR : public Basis1DAbstract
{
public:
    typedef EducatedBasisFunctorAbstract    functor_type;
    typedef boost::shared_ptr<functor_type> functor_ptrType;
    //! @name Constructors & Destructor
    //@{

    /*! The only available constructor
     */
    EducatedBasisDR () : M_icurrent (0) {};
    //@}

    //! @name Methods
    //@{

    /*! This method provide the next EigenValue of the problem
        Compute the zero of the Functor RR in the correct domain, moreover Next() updates the value of M_icurrent, M_istar
        and M_star for the next applying of Next(). So if you want the first 20 zeros of the functor RR you have to apply next
        20 times and save the results in your variable.
    */
    Real Next();
    void evaluateBasis (         Basis1DAbstract::MBMatrix_type& phi,
                                 Basis1DAbstract::MBMatrix_type& dphi,
                                 const Basis1DAbstract::MBVector_type& eigenvalues,
                                 const QuadratureRule* quadrule) const;
    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const {};

    Real evalSinglePoint ( const Real& eigen, const Real& yh) const;

    void setMu (const Real& mu)
    {
        M_mu = mu;
    }
    void setChi (const Real& chi)
    {
        M_chi = chi;
    }
    //@}

    Real chi() const
    {
        return M_chi;
    }

private:

    UInt M_icurrent;
    Real M_mu;
    Real M_chi;
    functor_ptrType M_ptrFunctor;

};

// ===================================================
// Macros
// ===================================================

//! define the ModalSpaceDDDD
inline
Basis1DAbstract* createDR()
{
    return new EducatedBasisDR();
}

namespace
{
static bool registerDR = Basis1DFactory::instance().registerProduct ( "dirrob",  &createDR);
}

} //End LifeV namespace
#endif //__EDUCATEDBASISDR_HPP__
