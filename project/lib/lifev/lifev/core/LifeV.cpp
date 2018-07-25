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
    @file
    @brief definition of nDimensions

    @author Simone Deparis <simone.deparis@epfl.ch>

    @date 2009-12-09

 */

#include <lifev/core/LifeV.hpp>

namespace LifeV
{
void Flag::showMe ( flag_Type const& flag, std::ostream& out )
{
    out << "Flag -- ";
    flag_Type bit = 0x01;
    for ( UInt i = 0; i < sizeof ( flag_Type ) * 8; i++ )
    {
        out << static_cast<bool> ( flag & bit );
        bit <<= 1;
    }
    out << std::endl;
}
const UInt nDimensions (NDIM);
} //end namespace LifeV

