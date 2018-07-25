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
 *   @file Projector.hpp
     @brief This class it is used into addrhs_highprec into himodassembler

     @date 06/2013
     @author M. Aletti <teo.aletti@gmail.com>
     @author A. Bortolossi <andrea.bortolossi@gmail.com>
 */
#ifndef __PROJECTOR__
#define __PROJECTOR__

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/himod/modalbasis/ModalSpaceAbstract.hpp>

namespace LifeV
{

class Projector
{

public:
    //! @name Pulbic Types
    //@{
    //! typedef for the returning type of the functor
    typedef Real return_Type;
    //! typedef for the type of the function to project
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
    //! typedef for the modalbasis
    typedef ModalSpaceAbstract                           modalbasis_type;
    //! typedef on the pointer to the modalbasis
    typedef boost::shared_ptr<modalbasis_type>              modalbasis_ptrType;
    //@}
    //! @name Constructor and Distructor
    //@{
    //! The constructor simply assign the pointer to the modalbasis
    Projector (modalbasis_ptrType modalbasisptr) : M_modalbasis (modalbasisptr) {};
    //! The distructor doesn't do anything
    ~Projector() {};
    //@}

    //! @name Setters
    //@{
    //! This methods has to be used to set the function one want to project
    void setFunction (const function_Type& f)
    {
        M_f = f;
    }
    //! This methods let the user to change the mode on which compute the Fourier Coefficient
    void setMode (const UInt& mode)
    {
        M_mode = mode;
    }
    //@}
    //! @name Operators
    //@{
    //! The actual implementation of the functor it is just a wrapper to FourierCoeffPointWise method in ModalBasis.hpp
    return_Type operator() (const VectorSmall<1>& position) const
    {
        return M_modalbasis->fourierCoeffPointWise (position[0], M_f, M_mode);
    }
    //@}

private:

    //! The function to project
    function_Type       M_f;
    //! The number of the base function
    UInt                M_mode;
    //! The pointer to the modal basis
    modalbasis_ptrType  M_modalbasis;
};

} // end LifeV namespace
#endif // end of __PROJECTOR__
