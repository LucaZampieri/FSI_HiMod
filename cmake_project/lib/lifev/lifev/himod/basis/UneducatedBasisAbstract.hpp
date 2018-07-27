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
 *   @file UneducatedBasisAbstract.hpp
     @brief This is an abstract class that provides the generic interface
     that a 2-dimensional uneducated basis should have. It derives from Basis2DAbstract.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>


 */
#ifndef __UNEDUCATEDBASISABSTRACT_HPP__
#define __UNEDUCATEDBASISABSTRACT_HPP__

#include <lifev/himod/basis/Basis2DAbstract.hpp>

namespace LifeV
{

class UneducatedBasisAbstract: public Basis2DAbstract
{
public:
    
    //! @name Constructor and Destructor
    //@{
    //! Constructor
    // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
    // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
    UneducatedBasisAbstract(): Basis2DAbstract(), M_normTheta( sqrt(2.) ) {};

    //! Destructor
    virtual ~UneducatedBasisAbstract() {};
    //@}
    
    // setMethods
     // Inherited from Basis2DAbstract:
     // setRho
     // setTheta
     // setNumberModes

    //! @name Get Methods
    //@{ 
    
    //! chi
    //! This method returns one of the two parameters associated to Robin BC, if any.
    virtual Real chi() const 
    {
        return 0;
    };
    
    //! eigenValues
    //! This method returns the vector containing M_mtot eigenvalues.
    virtual EigenContainer eigenValues() const
    {
        return M_eigenvalues;
    };
    
    virtual UInt mr1D() const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        return mr1D;
    };
    
    virtual UInt mtheta1D() const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );

        return mtheta1D;
    };
    
       virtual boost::shared_ptr<matrix_Type> massMatrix() const
    {
        return M_massMatrix;
    };

    //@}
    
    //! @name Methods
    //@{
    
    //! computeEigenvalues
    //! This method computes the eigenvalues of the bidimensional problem for boundary conditions.
    //! Eigenvalues depend directly on the kind of boundary condition, in fact the implementation
    //! is delegated to derived classes.
    virtual void computeEigenvalues() {}; 
    
    //! evaluateBasis
    //! This is method evaluates the (two-dimensional) modalbasis defined on 
    //! quadrature nodes, as well as its derivative.
    /*! @param    phirho            A reference to a MBMatrix_type.
                                phirho[p][q] stores the evaluation of the Bessel function associated to
                                the p-th eigenvalue on q-th quadrature node. As it is filled in runtime
                                by the method, it is actually an output.
        @param    dphirho            A reference to a MBMatrix_type.
                                dphirho[p][q] stores the evaluation of the derivative of the derivative of
                                the Bessel function associated to the p-th eigenvalue on q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    phitheta        A reference to a MBMatrix_type.
                                phitheta[p][q] stores the evaluation of the angular function (which is independent
                                from the eigenvalue) in the q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    dphitheta        A reference to a MBMatrix_type.
                                dphitheta[p][q] stores the evaluation of the derivative of the angular function
                                in the q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    eigenvalues        A reference to MBVector_type.
                                Vector of (two-dimensional) eigenvalues.
        @param    quadrulerho        A pointer to a QuadratureRule for rho-integration.
        @param    quadruletheta    A pointer to a QuadratureRule for theta-integration.
         @return                    void
    */
    void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const EigenContainer& eigenvalues,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta ) const {};

    virtual void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta
                        ) = 0;
                       
    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const = 0;

    virtual void
    gramSchmidtOrthogonalization( MBMatrix_type& phirho,
                                   MBMatrix_type& dphirho,
                                   const QuadratureRule* quadrulerho,
                                   const UInt mr1d,
                                   const UInt mtheta1d )
    {
        Real projection( 0 );
        Real norm( 0 );
    
        for( UInt p( 0 ); p != mtheta1d; ++p )
        {
            for( UInt k( 1 ); k != mr1d; ++k )
            {
                for( UInt pr( 0 ); pr != k; ++pr )
                {
                    // compute the projection
                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        projection += phirho[mr1d * p + k][n] *
                                        phirho[mr1d * p + pr][n] *
                                        quadrulerho->quadPointCoor( n, 0 ) *
                                        quadrulerho->weight( n );
                    }

                    // compute orthogonal function (not normalized yet)
                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        phirho[mr1d*p+k][n] -= phirho[mr1d*p+pr][n]*projection;
                        dphirho[mr1d*p+k][n] -= dphirho[mr1d*p+pr][n]*projection;
                    }
                    projection = 0;
                } // end pr-for
                
                // compute norm
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    norm += phirho[mr1d*p+k][n] * phirho[mr1d*p+k][n] *
                            quadrulerho->quadPointCoor( n, 0 ) *
                            quadrulerho->weight( n );
                }
                    
                // normalize
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    phirho[mr1d*p+k][n] /= sqrt( norm );
                    dphirho[mr1d*p+k][n] /= sqrt( norm );
                }
                    
                norm = 0;
            } // end k-for
        } // end p-for
};
    
    virtual void
    checkOrthogonality( const MBMatrix_type& phirho,
                        const QuadratureRule* quadrulerho,
                        const UInt mr1d,
                        const UInt mtheta1d, const bool verbose = 0 ) const
    {
        Real mass( 0 );
        Real res( 0 );
        
        if( verbose )
        {
            std::cout << " -------------------------------------------------------- " << std::endl;
            std::cout << "                     CHECKING ORTHOGONALITY                     " <<std::endl;
            std::cout<<std::endl;
        }
        
        for( UInt p( 0 ); p != mtheta1d; ++p )
        {
            for( UInt k( 0 ); k != mr1d; ++k )
            {
                for( UInt j( 0 ); j != mr1d; ++j )
                {
                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        mass += phirho[mr1d*p+k][n] * phirho[mr1d*p+j][n] *
                                quadrulerho->quadPointCoor( n, 0 ) *
                                quadrulerho->weight( n );
                    }
                    
                    // Cumulate error
                    if( mr1d*p+k == mr1d*p+j )
                    {
                        res += abs( mass - 1 );
                    }
                    else
                    {
                        res += abs( mass );
                    }
                    
                    mass = 0;
                }
            }
        } // end p-for
        
        if( verbose )
        {
            if( res <1e-10 )
            {
                std::cout << "====> TEST PASSED : sum of squared errors = " << res << std::endl;
            }
            else
            {
                std::cout << "====> TEST FAILED : sum of squared errors = " << res << std::endl;
            }
        
        
            std::cout<<std::endl;        
            std::cout << "                             END of TEST                      " <<std::endl;
            std::cout << " -------------------------------------------------------- " << std::endl;
            std::cout<<std::endl;
        }
    };
    
    //! evalSinglePoint
    //! This method evaluates the two-dimensional modal basis on an arbitrary point p.
    //! As the eigenvalue depends on the boundary condition, it is forced to be implemented in derived classes.
    /*!    @param    eigen        A reference to Real.
                            The eigenvalue associated to the modal function you want to compute.
        @param    index        A const UInt.
                            The index of the eigenvalue (0:M_mtot-1). It is needed beacause the coefficients 
                            of the modal function depend on it.
        @param    p            A reference to a const MBVector_type.
                            Polar coordinates (rho, theta) in (0,1)x(0,1) of the point of evaluation on the reference domain.
        @return                Real value of the evaluation of the (two-dimensional) 'index'-th modal function in point 'p'.
    */
    Real
    evalSinglePoint( const Real& eigen, const UInt index, const MBVector_type& p ) const {};


    virtual Real
    evalSinglePoint( const UInt index, const MBVector_type& p,
                       const QuadratureRule* quadrulerho, const MBMatrix_type& phirho
                        ) const
    {
        const int trigonometricBasis(0);
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        // Note: remind that when the basis is "domesticated" at boundary the 0-th index is dropped!
        UInt thetaIndex = static_cast<UInt> ( index / mr1D );
        UInt rIndex = index - mr1D * thetaIndex;

        Real y( p[0] );
        Real z( p[1] );

        Real rhoh = sqrt( y * y + z * z );
        Real thetah = atan( z / y );
        
        Real interp;
        
        // Interpolate the value
        UInt n( 0 );
        
        while( n != quadrulerho->nbQuadPt() && quadrulerho->quadPointCoor( n, 0 ) < rhoh )
        {
            ++n;
        }

        // N.B. La funzione è divisa per Rho per compensare la moltiplicazione  in evaluateBasis
        if( n == 0 )
        {
            interp = phirho[rIndex][0] / M_Rho;
        }
        else
        {
            interp = ( phirho[rIndex][n] / M_Rho - phirho[rIndex][n-1] / M_Rho ) /
                            ( quadrulerho->quadPointCoor( n, 0 ) - quadrulerho->quadPointCoor( n - 1, 0 ) ) *
                            ( rhoh - quadrulerho->quadPointCoor( n - 1, 0 ) ) +
                            phirho[rIndex][n-1] / M_Rho;
        }
        
        Real norm( 1. );
        switch( trigonometricBasis )
        {
            case 0:
                        if( thetaIndex == 0 )
                        {
                            norm = 1 / sqrt( 2 * M_PI );
                        }
                        else
                        {
                            norm = 1 / sqrt( M_PI );
                        }
                        return interp * norm * cos( thetaIndex * thetah );
            case 1:
                        return interp / sqrt( M_PI ) * sin( thetaIndex * thetah );
            default:
                    return 0;
        }
        return 0;
    }
    //@}
    
    virtual void rectangularTruncation( const UInt& mtot, UInt& mr1D, UInt& mtheta1D, const bool verbose = 0 ) const
    {
        mtheta1D = static_cast<UInt>( sqrt( M_mtot ) );
        mr1D = static_cast<UInt>( mtot / mtheta1D );
    
        while( mtheta1D > 1 && mr1D * mtheta1D != mtot )
        {
            mtheta1D -= 1;
            mr1D = static_cast<UInt>( mtot / mtheta1D );
        }
    
        if( verbose )
        {
            std::cout<<"Mtot = Mr * Mtheta : "<<mtot<<" = "<<mr1D<<" * "<<mtheta1D<<std::endl;
        }
        //UInt tmp = mr1D;
        //mr1D = mtheta1D;
        //mtheta1D = tmp;
        //mr1D=mtot; mtheta1D=1;

        return;
    };

    protected:

        //! normalization coefficient for goniometric basis
        Real M_normTheta;
        
};

typedef FactorySingleton< Factory < UneducatedBasisAbstract,  std::string> > UneducatedBasisFactory;

} // namespace LifeV

// Here the includes of object in the factory, we put them down here because the other options
// would have obliged the user to include all this header, that is too much complicated.
#include <lifev/himod/basis/Robert.hpp>
#include <lifev/himod/basis/ChebyshevLinearShift.hpp>
#include <lifev/himod/basis/ChebyshevQuadraticShift.hpp>
#include <lifev/himod/basis/ChebyshevQuadShift.hpp>
#include <lifev/himod/basis/ChebyshevQuadraticShiftNatural.hpp>
#include <lifev/himod/basis/ChebyshevQuadShiftNatural.hpp>
#include <lifev/himod/basis/Bessel.hpp>
#include <lifev/himod/basis/Zernike.hpp>
#include <lifev/himod/basis/ZernikeNatural.hpp>

#endif
