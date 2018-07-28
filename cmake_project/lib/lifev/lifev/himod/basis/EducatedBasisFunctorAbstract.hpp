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
    @file EducatedBasisFunctorAbstract.hpp
    @brief This is an abstract class that gives only the interface of the generic functor which hold
    the non-linear function whose zeros are the frequency of the problem.
    A functor was necessary because the function depends on some coefficients.
    remember that untill we have to use non-linear brent as a root finding algorithm you have to implement in
    the functor not f but f square.

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi

 */

#ifndef __EDUCATEDBASISFUNCTORABSTRACT__
#define __EDUCATEDBASISFUNCTORABSTRACT__

#include <lifev/core/LifeV.hpp>
//#include <boost/shared_ptr.hpp>

namespace LifeV
{
class EducatedBasisFunctorAbstract
{

public:
    //! @name Constructor and Distructor
    //@{
    //! The constructor
    EducatedBasisFunctorAbstract (const Real& mu, const Real& chi, const Real& L) : M_mu (mu), M_chi (chi), M_L (L) {};
    //! The distructor doesn't do anything
    virtual ~EducatedBasisFunctorAbstract() {};
    //@}

    //! @name Operators
    //@{
    //! The actual implementation of the functor
    virtual Real operator() (const Real& x) const = 0;
    //@}
protected:

    Real M_mu;
    Real M_chi;
    Real M_L;
};
} // End lifev namespace
#endif
