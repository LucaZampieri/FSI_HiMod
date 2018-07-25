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
    @brief Vector for elementary assembly

    @contributor Matteo Astorino <matteo.astorino@epfl.ch>
    @mantainer Matteo Astorino <matteo.astorino@epfl.ch>

 */

#include <lifev/core/array/VectorElemental.hpp>

namespace LifeV
{

VectorElemental::VectorElemental ( int nNode1, int nbr1 ) :
    super ( nNode1* nbr1 )
{
    _nBlockRow = nbr1;
    _nRow.resize ( _nBlockRow );
    _firstRow.resize ( _nBlockRow );
    int first = 0, n;
    for ( n = 0; n < nbr1; n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
}

VectorElemental::VectorElemental ( int nNode1, int nbr1,
                                   int nNode2, int nbr2 ) :
    super ( nNode1 * nbr1 + nNode2* nbr2 )
{
    _nBlockRow = nbr1 + nbr2;
    _nRow.resize ( _nBlockRow );
    _firstRow.resize ( _nBlockRow );
    int first = 0, n;
    for ( n = 0; n < nbr1; n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
    for ( n = nbr1; n < nbr1 + nbr2; n++ )
    {
        _nRow[ n ] = nNode2;
        _firstRow[ n ] = first;
        first += nNode2;
    }
}

VectorElemental::VectorElemental ( int nNode1, int nbr1,
                                   int nNode2, int nbr2,
                                   int nNode3, int nbr3 ) :
    super ( nNode1 * nbr1 + nNode2 * nbr2 + nNode3* nbr3 )
{
    _nBlockRow = nbr1 + nbr2 + nbr3;
    _nRow.resize ( _nBlockRow );
    _firstRow.resize ( _nBlockRow );
    int first = 0, n;
    for ( n = 0; n < nbr1; n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
    for ( n = nbr1; n < nbr1 + nbr2; n++ )
    {
        _nRow[ n ] = nNode2;
        _firstRow[ n ] = first;
        first += nNode2;
    }
    for ( n = nbr1 + nbr2; n < nbr1 + nbr2 + nbr3; n++ )
    {
        _nRow[ n ] = nNode3;
        _firstRow[ n ] = first;
        first += nNode3;
    }
}


void VectorElemental::showMe ( std::ostream& c )
{
    for ( int i = 0; i < _nBlockRow; i++ )
    {
        c << "Block (" << i << "), " << block ( i ) << std::endl;
    }
}
}
