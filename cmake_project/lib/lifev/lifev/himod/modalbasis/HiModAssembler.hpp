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
 *   @file HiModAssembler.hpp
     @brief This file contains a class for assemble an ADR problem with the Hierarchical Model Reduction method

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
 */

#ifndef __HIMODASSEMBLER__
#define __HIMODASSEMBLER__

namespace LifeV
{

enum sectionGeometry { Rectangular, Circular }; // rectangular=0, circular=1.

template< typename mesh_type, typename matrix_type, typename vector_type, int N>
class HiModAssembler;

}


#include <lifev/himod/modalbasis/HiModAssemblerRectangular.hpp>
#include <lifev/himod/modalbasis/HiModAssemblerCircular.hpp>

#endif
