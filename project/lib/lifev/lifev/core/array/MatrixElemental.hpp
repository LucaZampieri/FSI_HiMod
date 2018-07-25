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
    @brief Matrix for elementary assembly

    @contributor Matteo Astorino <matteo.astorino@epfl.ch>
    @mantainer Matteo Astorino <matteo.astorino@epfl.ch>

 */

#ifndef _ELEMMAT_H_INCLUDED
#define _ELEMMAT_H_INCLUDED

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/RNM.hpp>

namespace LifeV
{
class MatrixElemental
{

public:
    typedef KNM<Real> matrix_type;
    typedef KNM_<Real> matrix_view;
    //typedef Tab2d matrix_type;

    ~MatrixElemental();
    MatrixElemental ( UInt nNode1, UInt nbr1, UInt nbc1 ); // constructor for 1 finite element

    //! This is the constructor for the local matrix involving 2 finite elements
    /*!
      There are 6 arguments for this constructor. The arguements nbNode1 and nbNode2
      are simply the number of basis function for each finite element.

      The 4 others are used to describe the number of blocks on want in the matrix.

      Warning: all the blocks do not have the same size!

      Indeed, nbr1 describes the number of blocks that have a "height" of nbNode1,
      nbr2 the number of blocks that have a "height" of nbNode2, and the same with
      the columns for nbc1 and nbc2 (Think about it as a "local Stokes matrix": nbr1
      and nbc1 are 3 (components for the velocity), nbr2 and nbc2 are 1 (pressure
      is scalar)).
     */
    MatrixElemental ( UInt nNode1, UInt nbr1, UInt nbc1,
                      UInt nNode2, UInt nbr2, UInt nbc2 );
    MatrixElemental ( UInt nNode1, UInt nbr1, UInt nbc1,
                      UInt nNode2, UInt nbr2, UInt nbc2,
                      UInt nNode3, UInt nbr3, UInt nbc3 ); // constructor for 3 finite elements
    matrix_type& mat()
    {
        return _mat;
    }
    UInt nBlockRow() const
    {
        return _nBlockRow;
    }
    UInt nBlockCol() const
    {
        return _nBlockCol;
    }

    //Tab2dView block( UInt i, UInt j )
    matrix_view block ( UInt i, UInt j )
    {
        return _mat ( SubArray ( _nRow[ i ], _firstRow[ i ] ),
                      SubArray ( _nCol[ j ], _firstCol[ j ] ) );
        //Tab2dView __mr (_mat, TabRange(_firstRow[i], _nRow[i]), TabRange(_firstCol[j], _nCol[j]));
        //return __mr;
    }
    void zero()
    {
        //_mat = ZeroMatrix( _mat.size1(), _mat.size2() );
        _mat = 0.0;
    };
    void showMe ( std::ostream& c = std::cout );

    void   operator *= (Real coef)
    {
        this->_mat *= coef;
    }
private:

    matrix_type _mat;

    UInt _nBlockRow; // number of block rows
    UInt _nBlockCol; // number of block columns
    std::vector<UInt> _nRow; // _nRow[i]=nb of rows in the i-th block row
    std::vector<UInt> _firstRow; //_firstRow[i]=index of first row of i-th block row
    std::vector<UInt> _nCol; // _nCol[i]=nb of col in the i-th block col
    std::vector<UInt> _firstCol; //_firstCol[i]=index of first col of i-th block col
};
}
#endif


