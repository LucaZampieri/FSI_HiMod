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
    @file EducatedBasisFunctorImplemented.hpp
    @brief This file contains all the different kind of functors needed to provide all the different combinations
    of BC. We put it all here, beacuse they are very short and we preferred to not produce many different, very short, files.

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi

 */

#ifndef __EDUCATEDBASISFUNCTORIMPLEMENTED_HPP__
#define __EDUCATEDBASISFUNCTORIMPLEMENTED_HPP__

#include <lifev/core/LifeV.hpp>

#include <lifev/himod/basis/EducatedBasisFunctorAbstract.hpp>
namespace LifeV
{
class EducatedBasisFunctorRR : public EducatedBasisFunctorAbstract
{
public:
    EducatedBasisFunctorRR (const Real& mu, const Real& chi, const Real& L) : EducatedBasisFunctorAbstract (mu, chi, L) {};
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() (const Real& x) const
    {
        return std::pow (2.0 * M_mu / M_L * x + std::tan (x) * (M_chi - M_mu / M_L * M_mu / M_L * x * x / M_chi), 2);
    }
    //@}
};

class EducatedBasisFunctorDR : public EducatedBasisFunctorAbstract
{
public:
    EducatedBasisFunctorDR (const Real& mu, const Real& chi, const Real& L) : EducatedBasisFunctorAbstract (mu, chi, L) {};
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() (const Real& x) const
    {
        return std::pow (M_mu / M_chi * x / M_L + std::tan (x), 2);
    }
    //@}
};

class EducatedBasisFunctorRD : public EducatedBasisFunctorAbstract
{
public:
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() (const Real& /*x*/ ) const
    {
        ERROR_MSG ("this type of BC ha not been implemented yet, please add it!");
        return -1;
    }
    //@}
};
} // End lifev namespace
#endif
