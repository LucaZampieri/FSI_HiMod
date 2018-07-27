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
 *   @file
     @brief This file contains the definition of the EvaluationScalar class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef EVALUATION_SCALAR_HPP
#define EVALUATION_SCALAR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionScalar.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for a scalar constant
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a scalar constant during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
class EvaluationScalar
{
public:

    //! @name Public Types
    //@{

    //! Type returned by this class
    typedef Real return_Type;

    //@}


    //! @name Static constants
    //@{

    //! Flag for the global current FE
    const static flag_Type S_globalUpdateFlag;

    //! Flag for the test current FE
    const static flag_Type S_testUpdateFlag;

    //! Flag for the solution current FE
    const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    EvaluationScalar()
    {}

    //! Copy constructor
    EvaluationScalar (const EvaluationScalar& evaluation)
        : M_value (evaluation.M_value)
    {}

    //! Expression-based constructor
    explicit EvaluationScalar (const ExpressionScalar& expression)
        : M_value (expression.value() )
    {}

    //! Destructor
    ~EvaluationScalar()
    {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/)
    {}

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "scalar";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/)
    {}

    //! Do nothing setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* /*testCFE*/)
    {}

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* /*solutionCFE*/)
    {}

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&)
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& /*q*/) const
    {
        return M_value;
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& /*q*/, const UInt& /*i*/) const
    {
        return M_value;
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& /*q*/, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_value;
    }

    //@}

private:

    // Value stored
    Real M_value;
};


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
