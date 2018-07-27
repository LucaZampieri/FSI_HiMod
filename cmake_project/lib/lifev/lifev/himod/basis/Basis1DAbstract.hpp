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
 *   @file Basis1DAbstract.hpp
     @brief This is an abstract class that provides the generic interface
     that a monodimensional basis should have.

     @date 06/2013
     @author M. Aletti <teo.aletti@gmail.com>
     @author A. Bortolossi


 */
#ifndef __BASIS1DABSTRACT_HPP__
#define __BASIS1DABSTRACT_HPP__

#include <lifev/himod/basis/BasisAbstract.hpp>
#include <lifev/himod/basis/EducatedBasisFunctorImplemented.hpp>

namespace LifeV
{

class Basis1DAbstract:public BasisAbstract
{
public:

    Basis1DAbstract() {};

    virtual ~Basis1DAbstract() {};
    
    //! @name setMethods
    //@{
    //! With this method you can set the length of the domain, altough this class is meant to
    //! work on the reference domain in some cases the length of physical domain is still needed to adjust
    //! the boundary conditions coefficients.
    void setL( const Real& L )
    {
        M_L = L;
    };
    //@}



    //! @name Methods
    //@{
    //*! Any time you call this method it should return the next eigenvalue of the problem
    //! pay attention that even this class is meant to work on the reference domain (0,1)
    //! this method should return the eigenvalue of the physical problem, so you have to properly normalize it.
    virtual Real Next() = 0;

    //! This is the method that actually evaluates the modalbasis defined on the reference domain along with its
    //! derivative.
    virtual void evaluateBasis(     MBMatrix_type& phi,
                                    MBMatrix_type& dphi,
                                    const MBVector_type& eigenvalues,
                                    const QuadratureRule* quadrule      ) const = 0;

    //! This method takes only two parameters: a point in (0,1), the reference domain, and the
    //! frequency of sinusoidal function
    virtual Real evalSinglePoint ( const Real& eigen, const Real& yh) const = 0;
    //@}

protected:

    // The physical length of the 1D-domain.
    Real M_L;

};


// ===================================================
// Macros
// ===================================================

//! create factory for ModalSpace
typedef FactorySingleton< Factory < Basis1DAbstract,  std::string> > Basis1DFactory;

} // namespace LifeV

// Here the includes of object in the factory, we put them down here because the other options
// would have obliged the user to include all this header, that is too much complicated.
#include <lifev/himod/basis/EducatedBasisDD.hpp>
#include <lifev/himod/basis/EducatedBasisRR.hpp>
#include <lifev/himod/basis/EducatedBasisDR.hpp>
#include <lifev/himod/basis/EducatedBasisNN.hpp>
#include <lifev/himod/basis/FakeBasis.hpp>

#endif
