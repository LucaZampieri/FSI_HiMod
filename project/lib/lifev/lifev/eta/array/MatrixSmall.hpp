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
 *      @file
         @brief This file contains a simple matrix class

        @contributor Ivan Kuraj <ivan.kuraj@epfl.ch>
 */

#ifndef _MATRIXSMALL_H_
#define _MATRIXSMALL_H_ 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshVertex.hpp>
#include <lifev/core/array/RNM.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace std;

// macros for defining actions on an out-of-bound reference
#define MATRIX_SMALL_DIMENSION_CHECK_NO_CHECK 0
#define MATRIX_SMALL_DIMENSION_CHECK_ASSERT 1
#define MATRIX_SMALL_DIMENSION_CHECK_EXCEPTION 2
#define MATRIX_SMALL_DIMENSION_CHECK MATRIX_SMALL_DIMENSION_CHECK_NO_CHECK

// LifeV namespace.
namespace LifeV
{
//! class MtrixSmall   This class implements a simple matrix

/*!
  @author Ivan Kuraj <ivan.kuraj@epfl.ch>

  This class implements a simple matrix.
  <br>
  It allows all kind of geometric operations on the node,
  such as summation, multiplication by scalar, scalar product,
  cross product, etc.
  The implementation is oriented to best perform with small (less than 30) n

*/

template <UInt Dim1, UInt Dim2>
class MatrixSmall
{

private:
    typedef Real& (MatrixSmall::*DereferenceMethod) ( UInt const& i );
    typedef const Real& (MatrixSmall::*ConstDereferenceMethod) ( UInt const& i ) const;

    class Dereferencer
    {
    public:
        Dereferencer (UInt iRow, MatrixSmall<Dim1, Dim2> const& iMatrix) : row (iRow), matrix (iMatrix) {}
        Real& operator[] ( UInt const& i )
        {
#if MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_ASSERT
            ASSERT ( i < Dim2, "trying to access an index that exceeds the second dimension of the matrix" );
#elif MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_EXCEPTION
            if (i >= Dim2)
            {
                stringstream messageStream (stringstream::out);
                messageStream << "Trying to access an index [" << row << "][" << i << "] that exceeds the second dimension of the matrix " << Dim2;
                throw range_error (messageStream.str() );
            }
#endif
            return (const_cast<Real&> (matrix.M_coords [ row ] [ i ]) );
        }

        const Real& operator[] ( UInt const& i ) const
        {
            return ( (const_cast<Dereferencer&> (*this) ) [i]);
        }

        UInt row;
        const MatrixSmall& matrix;
    };

#if MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_NO_CHECK
    typedef Real* OpIndexReturnType;
#else
    typedef Dereferencer OpIndexReturnType;
#endif

    void copyFrom ( MatrixSmall<Dim1, Dim2> const& matrix )
    {
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                M_coords[i][j] = matrix.M_coords[i][j];
            }
        //reinterpret_cast<Real*>(M_coords)[ i ]  = reinterpret_cast<const Real*>(matrix.M_coords)[ i ];
    }

    bool realEquality (const Real& r1, const Real& r2) const
    {
        // TODO change if it is not appropriate
        //return (fabs(r1 - r2) <= numeric_limits<double>::epsilon());
        //return (r1==r2);
        return (fabs (r1 - r2) <= 0.1);
    }

public:

    //! @name Constructors and destructors
    //@{

    //! Empty constructor (all components are set to zero)
    MatrixSmall()
    {
        for ( UInt i = 0; i < Dim1 * Dim2; i++ )
        {
            reinterpret_cast<Real*> (M_coords) [ i ] = 0.;
        }
    }

    //! Copy constructor
    MatrixSmall ( MatrixSmall<Dim1, Dim2> const& matrix )
    {
        copyFrom (matrix);
    }

    //! Import from a vector
    MatrixSmall (const vector< vector<Real> >& matrix)
    {
        // check if dimensions are correct
        bool isDim2Appropriate = true;
        for (unsigned i = 0; i < matrix.size(); i++)
            if (matrix[i].size() != Dim2)
            {
                isDim2Appropriate = false;
                break;
            }
        ASSERT ( matrix.size() == Dim1 &&  isDim2Appropriate,
                 " Non matching dimension for the construction of a fixed size matrix via vector ");

        for (unsigned int i = 0; i < Dim1; ++i)
        {
            for (unsigned int j = 0; j < Dim2; ++j)
            {
                M_coords[i][j] = matrix[i][j];
            }
        }
    }

    //@}

    //! @name Overloaded operators
    //@{

    //! Operator ==
    bool operator== ( MatrixSmall<Dim1, Dim2> const& matrix ) const
    {
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
                if (!realEquality (M_coords[i][j], matrix.M_coords[i][j]) )
                {
                    return false;
                }
        return true;
    }

    //! Operator !=
    bool operator!= ( MatrixSmall<Dim1, Dim2> const& matrix ) const
    {
        return ! (*this == matrix);
    }

    //! Assignment operator
    MatrixSmall<Dim1, Dim2>& operator= ( MatrixSmall<Dim1, Dim2> const& matrix )
    {
        // avoid this check for fastest common case
        //if (this != &matrix)
        copyFrom (matrix);
        return (*this);
    }

    //! Operator +=
    MatrixSmall<Dim1, Dim2>& operator+= ( MatrixSmall<Dim1, Dim2> const& matrix )
    {
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                M_coords[ i ][ j ] += matrix.M_coords[ i ][ j ];
            }
        return (*this);
    }

    //! Operator +
    MatrixSmall<Dim1, Dim2> operator+ ( MatrixSmall<Dim1, Dim2> const& matrix ) const
    {
        MatrixSmall<Dim1, Dim2> tmp ( *this );
        return (tmp += matrix);
    }

    //! Operator -=
    MatrixSmall<Dim1, Dim2>& operator-= ( MatrixSmall<Dim1, Dim2> const& matrix )
    {
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                M_coords[ i ][ j ] -= matrix.M_coords[ i ][ j ];
            }
        return (*this);
    }

    //! Operator -
    MatrixSmall<Dim1, Dim2> operator- ( MatrixSmall<Dim1, Dim2> const& matrix ) const
    {
        MatrixSmall<Dim1, Dim2> tmp ( *this );
        return (tmp -= matrix);
    }

    //! Operator *= (multiplication by scalar)
    MatrixSmall<Dim1, Dim2>&   operator*= ( Real const& factor )
    {
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                M_coords[ i ][ j ] *= factor;
            }
        return (*this);
    }

    //! Operator /= (division by scalar)
    MatrixSmall<Dim1, Dim2>& operator/= ( Real const& factor )
    {
        ASSERT ( factor != 0. , "Division by zero!" );
        *this *= 1. / factor;
        return (*this);
    }

    //! Operator / (division by scalar)
    MatrixSmall<Dim1, Dim2> operator/ ( Real const& factor ) const
    {
        MatrixSmall<Dim1, Dim2> tmp ( *this );
        return (tmp /= factor);
    }

    //! Operator / (division by scalar)
    MatrixSmall<Dim1, Dim2> operator* ( Real const& factor ) const
    {
        MatrixSmall<Dim1, Dim2> tmp ( *this );
        return (tmp *= factor);
    }

    //! Operator * (multiplication by a vector)
    VectorSmall<Dim1>  operator* ( VectorSmall<Dim2> const& vector ) const
    {
        VectorSmall<Dim1> resultVector;
        for ( UInt i = 0; i < Dim1; i++ )
        {
            resultVector[i] = 0;
            for ( UInt j = 0; j < Dim2; j++ )
            {
                resultVector[i] += vector[j] * (*this) [i][j];
            }
        }
        return (resultVector);
    }

    //! Operator []
    //const OpIndexReturnType operator[] ( UInt const & i ) const
    OpIndexReturnType operator[] ( UInt const& i ) const
    {
        return ( (const_cast<MatrixSmall<Dim1, Dim2>&> (*this) ) [i]);
    }

    //! Operator []
    OpIndexReturnType operator[] ( UInt const& i )
    {
#if MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_ASSERT
        ASSERT ( i < Dim1, "trying to access an index that exceeds the first dimension of the matrix" );
#elif MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_EXCEPTION
        stringstream messageStream (stringstream::out);
        messageStream << "Trying to access an index[" << i << "] that exceeds the second dimension of the matrix " << Dim1;
        if (i >= Dim1)
        {
            throw range_error (messageStream.str() );
        }
#endif
#if MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_NO_CHECK
        return (M_coords [ i ]);
#else
        return (Dereferencer (i, *this) );
#endif
    }

    //! Operator ()
    Real const& operator() ( UInt const& i, UInt const& j ) const
    {
#if MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_ASSERT
        ASSERT ( i < Dim1 && j < Dim2, "trying to access an index that exceeds the first dimension of the matrix" );
#elif MATRIX_SMALL_DIMENSION_CHECK  == MATRIX_SMALL_DIMENSION_CHECK_EXCEPTION
        stringstream messageStream (stringstream::out);
        messageStream << "Trying to access an index [" << i << "][" << j << "] that is out of bounds.";
        if (i >= Dim1 || j >= Dim2)
        {
            throw range_error (messageStream.str() );
        }
#endif
        return (M_coords[i][j]);
    }

    //! Operator ()
    Real& operator() ( UInt const& i, UInt const& j )
    {
        return (const_cast<Real&> (const_cast<const MatrixSmall<Dim1, Dim2>&> (*this) (i, j) ) );
    }

    //@}

    //! @name Geometric Methods
    //@{

    //! Scalar product
    /*!
    @param matrix second operand
    @return scalar product value
    */
    Real dot ( MatrixSmall<Dim1, Dim2> const& matrix ) const
    {
        Real scalarProduct = 0.;
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                scalarProduct += M_coords[ i ][ j ] * matrix.M_coords[ i ][ j ];
            }
        return scalarProduct;
    }

    //! Norm value
    /*!
    @return norm value
    */
    Real norm () const
    {
        return std::sqrt ( this->dot ( *this ) );
    }

    //! Normalize matrix
    void normalize ()
    {
        *this /= norm ();
    }

    //! Create the versor associated to this MatrixSmall
    /*!
    @return the versor associated to this MatrixSmall
    */
    MatrixSmall<Dim1, Dim2> normalized ()
    {
        return ( ( *this ) / norm () );
    }

    //@}

    //! @name Tools
    //@{

    //! function to get the size of the MatrixSmall ( for compatibility with Eigen)
    /*!
    @return the fixed size of the MatrixSmall
    */
    // we avoid the size method for matrix
    //    static UInt size() { return Dim1 * Dim2;}

    //@}

    //! @name Output stream operator overload
    //@{
    friend ostream& operator<< ( ostream& out, MatrixSmall<Dim1, Dim2> const& matrix )
    {
        out << "(" << std::endl ;
        for ( UInt i = 0; i < Dim1; i++ )
        {
            for ( UInt j = 0; j < Dim2; j++ )
            {
                out << matrix.M_coords[i][j] << " ";
            }
            out << std::endl;
        }
        out << ")";
        return (out);
    }
    //@}

    //! @name Conversion free-functions
    //@{

    //! Conversion of an array (std::vector, KN, etc. if applicable) to a MatrixSmall
    /*!
    @param coords matrix of point coordinates with operator[][] available
    @return the matrixSmall that corresponds to the input
    */
    template <typename Matrix>
    friend inline MatrixSmall<Dim1, Dim2> castToMatrixSmall ( Matrix const& coords )
    {
        MatrixSmall<Dim1, Dim2> tmp;
        for ( UInt i = 0; i < Dim1; i++ )
            for ( UInt j = 0; j < Dim2; j++ )
            {
                tmp.M_coords[ i ][ j ] = coords[ i ][ j ];
            }
        return (tmp);
    }
    //@}
    //
private:

    //! @name Data
    //@{

    //! Data storage
    Real M_coords[ Dim1 ] [ Dim2 ];

    //@}
};

//! @name External overloaded operators
//@{

//! Operator * (multiplication by scalar on the left)
template <UInt Dim1, UInt Dim2>
inline MatrixSmall<Dim1, Dim2> operator* ( Real const& factor, MatrixSmall<Dim1, Dim2> const& matrix )
{
    MatrixSmall<Dim1, Dim2> tmp ( matrix );
    return (tmp *= factor);
}

//! Operator * (multiplication by vector on the left)
template <UInt Dim1, UInt Dim2>
inline VectorSmall<Dim1> operator* ( VectorSmall<Dim2> const& vector, MatrixSmall<Dim1, Dim2> const& matrix )
{
    return (matrix * vector);
}

//@}

} // namespace LifeV


#endif //_MATRIXSMALL_H_

// -*- mode: c++ -*-
