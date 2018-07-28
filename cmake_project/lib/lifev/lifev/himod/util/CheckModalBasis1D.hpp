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
    @file CheckModalBasis1D.hpp
    @brief A class that is very simple tool to check if the given basis is orthonormal and if it satisfies the BC.

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi@gmail.com>
    @contributor M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
    @contributor S. Guzzetti <sofia.guzzetti@gmail.com>

 */
#ifndef __CHECKMODALBASIS1D_HPP_
#define __CHECKMODALBASIS1D_HPP_

#include <lifev/core/LifeV.hpp>
#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>

namespace LifeV
{

class CheckModalBasis1D
{

public:

    typedef ModalSpaceRectangular                   modalbasis_type;
    typedef boost::shared_ptr<modalbasis_type>      modalbasis_ptrType;

    CheckModalBasis1D (const modalbasis_ptrType& modalbasis);
    void VerifyOrthonormality() const;
    void VerifyBC (const Real& mu = 1, const Real& Chi = 1) const;

private:

    modalbasis_ptrType M_modalbasis;
};

} // namespace LifeV

#endif // __CheckModalBasis1D__
