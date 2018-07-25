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
 *   @file EducatedBasisNN.hpp
     @brief This file contains the interface of an educated basis (monodimensional) with Dirichlet conditions on both side.
            \f[ u' = 0 \f] on both side.

     @date 06/2013
     @author M. Aletti <teo.aletti@gmail.com>
     @author A. Bortolossi <andrea.bortolossi89@gmail.com>


 */
#ifndef __EDUCATEDBASISNN_HPP__
#define __EDUCATEDBASISNN_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/himod/basis/Basis1DAbstract.hpp>

//#include <boost/shared_ptr.hpp>

namespace LifeV
{

class EducatedBasisNN : public Basis1DAbstract
{
public:
    //! @name Constructors & Destructor
    //@{

    /*!
    The only available constructor, this class should be used only through the factory so it need only to be
    default-constructable.
     */
    EducatedBasisNN() : Basis1DAbstract(), M_icurrent (0) {};
    //@}

    //! @name Methods
    //@{
    //! This method provide the next EigenValue of the problem
    Real Next()
    {
        return M_icurrent++ * M_PI / M_L;
    };

    void evaluateBasis (         Basis1DAbstract::MBMatrix_type& phi,
                                 Basis1DAbstract::MBMatrix_type& dphi,
                                 const Basis1DAbstract::MBVector_type& eigenvalues,
                                 const QuadratureRule* quadrule) const;

    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const {};

    Real evalSinglePoint ( const Real& eigen, const Real& yh) const;
    //@}

    Real chi() const
    {
        return 0.0;
    };

private:

    UInt M_icurrent;
};

// ===================================================
// Macros
// ===================================================

inline
Basis1DAbstract* createNN()
{
    return new EducatedBasisNN() ;
}

namespace
{
static bool registerNN = Basis1DFactory::instance().registerProduct ( "neuneu",  &createNN);
}


} //End LifeV namespace
#endif
