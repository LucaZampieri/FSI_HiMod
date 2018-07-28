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
 *   @file BasisAbstract.hpp
     @brief This is an abstract class that provides the generic interface
     that a 1D and 2D basis should have.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>


 */

#ifndef __BASISABSTRACT_HPP__
#define __BASISABSTRACT_HPP__

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <Epetra_SerialComm.h>


namespace LifeV
{

class BasisAbstract
{
public:

	typedef MatrixEpetraStructured<Real>        		 matrix_Type;
    typedef std::vector<std::vector<Real> >              MBMatrix_type;
    typedef std::vector<Real>                            MBVector_type;

	//!@name Constructor and Destructor
	//@{
	//! Constructor
    BasisAbstract() {};

	//! Virtual destructor
    virtual ~BasisAbstract() {};
    //@}

    //! @name Set Methods
    //@{
    //! With this method you can set the coefficient mu, typically needed when you apply a Robin condition
    // so remember to overload this methods anytime you need it. Viceversa if you do not need to set mu
    // this simply do nothing.
    virtual void setMu (Real const& /*mu*/) {};
 
    //! With this method you can set the coefficient chi, typically needed when you apply a Robin condition
    // so remember to overload this methods anytime you need it. Viceversa if you do not need to set mu
    // this simply do nothing.
    virtual void setChi (Real const& /*chi*/) {};
    //@}

    //! @name Methods
    //@{
    virtual Real chi() const = 0;
    //@}

protected:

};
}

#endif

