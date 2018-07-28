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
 *   @file HiModPostProcessor.hpp
     @brief This file contains the utilities for the post-processing of the HiMod solution.

     @date 02/2017
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
 */

#ifndef __HIMODPOSTPROCESSOR_HPP__
#define __HIMODPOSTPROCESSOR_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/himod/tools/ReferenceMap.hpp>
#include <lifev/himod/modalbasis/NSHiModAssembler.hpp>

namespace LifeV
{

template< typename mesh_type, typename matrix_type, typename vector_type, int N>
class HiModPostProcessor
{
public:
   
    typedef    boost::shared_ptr<vector_type>    vector_ptrType;
    typedef MapEpetra                                       map_type;
    typedef typename map_type::comm_ptrtype                 commPtr_Type;
    typedef    NSHiModAssembler< mesh_type , matrix_type , vector_type, N > HMA_type;
    typedef    std::vector< std::vector<Real> >    stdMatrixType;
    typedef    std::vector<Real>                   stdVectorType;
    typedef    boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
        Default
     */
    HiModPostProcessor( const HMA_type& HMAsin,
                        const HMA_type& HMAcos,
                        const map_type& Map, commPtr_Type& Comm  ) :
                        M_HMAsin( HMAsin, Comm ), M_HMAcos( HMAcos, Comm ),
                        M_solSin( new vector_type ( Map, Unique ) ),
                        M_solCos( new vector_type ( Map, Unique ) ){};
    
    HiModPostProcessor( const boost::shared_ptr<vector_type>& solSin,
                        const boost::shared_ptr<vector_type>& solCos,
                        const HMA_type& HMAsin,
                        const HMA_type& HMAcos,
                        commPtr_Type& Comm  ) :
                        M_HMAsin( HMAsin, Comm ), M_HMAcos( HMAcos, Comm ),
                        M_solSin( solSin ), M_solCos( solCos ) {};
    
    ~HiModPostProcessor(){};
    //@}
    
    //! @name Methods
    //@{
    // Evaluate HiMod solution on the HiMod grid (quadrature nodes and FE nodes)
    vector_type evaluateSolutionHMgrid( const map_type& Map_3D )
    {
        vector_type solution_3D( Map_3D,Unique );
        solution_3D = M_HMAsin.evaluateBase3DGrid( *M_solSin );
        solution_3D += M_HMAcos.evaluateBase3DGrid( *M_solCos );
        return solution_3D;
    }; 

    // Evaluate HiMod solution on given points
    stdVectorType evaluateSolution3Dpoints( const GetPot& datafile, std::string& filename );

    // Read file
    void readMatrix( const std::string fileName, const UInt& nR, const UInt& nC, stdMatrixType& M );

    // Convert polar components to cartesian components
    stdVectorType PolarVec2CartesianVec( const stdVectorType& polarEvals, 
                                         const stdVectorType& x, const stdVectorType& y, const stdVectorType& z );
    //@}
    
    // SET METHODS
    void setHMsolution( const vector_ptrType& sinSol, const vector_ptrType& cosSol )
    {
        M_solSin = sinSol;
        M_solCos = cosSol;
    };
    
private:

    vector_ptrType           M_solSin;
    vector_ptrType           M_solCos;
    HMA_type                 M_HMAsin;
    HMA_type                 M_HMAcos;
        
};

template< typename mesh_type, typename matrix_type, typename vector_type, int N>
typename HiModPostProcessor<mesh_type, matrix_type, vector_type, N>::stdVectorType 
HiModPostProcessor<mesh_type, matrix_type, vector_type, N>::
evaluateSolution3Dpoints( const GetPot& data, std::string& filename )
{
    std::string coordFile( data("mesh/points_file", "points.txt" ) );
    UInt Npoints( data( "mesh/Npoints", 0 ) );
    stdMatrixType evalCoords;

    readMatrix( coordFile, 3, Npoints, evalCoords );

    stdVectorType x( evalCoords[0] );
    stdVectorType y( evalCoords[1] );
    stdVectorType z( evalCoords[2] );

    stdVectorType evals( 4*Npoints, 0 );
    stdVectorType evalsCos( 4*Npoints, 0 );
    stdVectorType evalsSin( 4*Npoints, 0 );
    evalsSin = M_HMAsin.evaluateHiModFunc( M_solSin, x, y, z );
    evalsCos = M_HMAcos.evaluateHiModFunc( M_solCos, x, y, z );
    
    std::transform( evalsCos.begin(), evalsCos.end(), evalsSin.begin(), evals.begin(), std::plus<Real>() );
    stdVectorType cartEvals = PolarVec2CartesianVec( evals, x, y, z );
    std::ofstream out( filename.c_str() );
    std::ostream_iterator<Real> it( out, "\n" );
    std::copy( cartEvals.begin(), cartEvals.end(), it );

    return evals;
}

template< typename mesh_type, typename matrix_type, typename vector_type, int N>
typename HiModPostProcessor<mesh_type, matrix_type, vector_type, N>::stdVectorType 
HiModPostProcessor<mesh_type, matrix_type, vector_type, N>::
PolarVec2CartesianVec( const stdVectorType& polarEvals, 
                       const stdVectorType& x, const stdVectorType& y, const stdVectorType& z )
{
    UInt Npts( x.size() );
    stdVectorType cartEvals( 4*Npts, 0 );

    for( UInt i(0); i<Npts; ++i )
    {
        Real theta = std::atan( z[i]/y[i] ); // theta \in [-pi/2,pi/2]
        if( y[i] < 0 && z[i] >= 0 )
        {
            theta = theta + M_PI;
        }
        else if( y[i] <= 0 && z[i] < 0 )
        {
            theta = theta + M_PI;
        }
        else if( y[i] > 0 && z[i] < 0 )
        {
            theta = theta + 2*M_PI;
        }
        cartEvals[i]        = polarEvals[i];
        cartEvals[Npts+i]   = polarEvals[Npts+i]*cos(theta) - polarEvals[2*Npts+i]*sin(theta);
        cartEvals[2*Npts+i] = polarEvals[Npts+i]*sin(theta) + polarEvals[2*Npts+i]*cos(theta);
        cartEvals[3*Npts+i] = polarEvals[3*Npts+i];
    }
    return cartEvals;
}

template< typename mesh_type, typename matrix_type, typename vector_type, int N> void
HiModPostProcessor<mesh_type, matrix_type, vector_type, N>::
readMatrix( const std::string fileName, const UInt& nR, const UInt& nC, stdMatrixType& M )
{
    // Define file and resize containers
    std::ifstream file( fileName.c_str() );
    std::vector<Real> Mvec;
    M.resize( nR );
        
    for( UInt i( 0 ); i != nR; ++ i )
    {
        M[i].resize( nC );
    }
        
    // Read matrix in vector form
    std::copy( std::istream_iterator<Real> (file),
               std::istream_iterator<Real> (),
               std::back_inserter< std::vector<Real> > (Mvec) );
               
    for( UInt k(0); k != nR*nC; ++k )
    {
        UInt col( k%nC );
        UInt row( k/nC );
            
        M[row][col] = Mvec[k];
    }
        
    return;
}

} // end namespace


#endif
