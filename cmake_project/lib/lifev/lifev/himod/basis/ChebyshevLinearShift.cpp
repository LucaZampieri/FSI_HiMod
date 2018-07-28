#include <lifev/himod/basis/ChebyshevLinearShift.hpp>
#include <cmath>

namespace LifeV
{
    void
    ChebyshevLinearShift::
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta )
    {
        int trigonometricBasis(0);
            UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        M_eigenvalues.resize( M_mtot );
        phirho.resize( M_mtot );
        dphirho.resize( M_mtot );
        phitheta.resize( M_mtot );
        dphitheta.resize( M_mtot );

        for ( UInt i = 0; i != M_mtot; ++i )
        {
            phirho[i].resize( quadrulerho->nbQuadPt() );
            dphirho[i].resize( quadrulerho->nbQuadPt() );
            phitheta[i].resize( quadruletheta->nbQuadPt() );
            dphitheta[i].resize( quadruletheta->nbQuadPt() );

        }
        
        MBVector_type normR( mr1D, 0 );
        for( UInt i = 0; i != mr1D; ++i )
        {
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real x = 2 * quadrulerho->quadPointCoor( n, 0 ) - 1;
                normR[i] += cos( i * acos( x ) ) * cos( i * acos( x ) ) * M_Rho * quadrulerho->quadPointCoor( n, 0 ) * 
                                M_Rho * quadrulerho->weight( n );
            }
        }
        
        for ( UInt pt = 0; pt != mtheta1D; ++pt )
        {
            for ( UInt pr = 0; pr != mr1D; ++pr )
            {
                M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                M_eigenvalues[pt*mr1D+pr].order = pt;
                M_eigenvalues[pt*mr1D+pr].index = pr;
                    
                for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                {
                        Real x = 2 * quadrulerho->quadPointCoor( n, 0 ) - 1;
                    // NOTA BENE: LE FUNZIONI SONO MOLTIPLICATE PER UN M_Rho IN PIÙ PER COMPENSARE IL normrho CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                        phirho[mr1D*pt + pr][n] = M_Rho *
                                                    cos( pr * acos( x ) ) / sqrt( normR[pr] );
                                                    
                        dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                                2 * pr * sin( pr * acos( x ) ) /
                                                sqrt( 1 - x * x ) ) * 1 / sqrt( normR[pr] );
                } // end n-for
            
                switch( trigonometricBasis )
                {
                    case 0:
                            for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                            {
                                    // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                                    // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
                                    phitheta[mr1D*pt + pr][n]  = M_normTheta * ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                    dphitheta[mr1D*pt + pr][n] = M_normTheta * pt *
                                                                      ( - sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                            }        
                            break;                            
                    case 1:
                            for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                            {
                                phitheta[mr1D*pt + pr][n]  = M_normTheta * ( sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                dphitheta[mr1D*pt + pr][n] = M_normTheta * pt *
                                                                  ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                            }
                            break;
                    default:
                            break;
                }
            } // end p-for
        }
        
        return;
    }
    
    /*Real
    ChebyshevLinearShift::
    evalSinglePoint( const UInt index, const MBVector_type& p,
                       const QuadratureRule* quadrulerho, const int trigonometricBasis ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        UInt rIndex = remainder( index, mr1D );
        UInt thetaIndex = static_cast<UInt> ( index / mr1D );
        
        Real rhoh = p[0];
        Real thetah = p[1];
        
        switch( trigonometricBasis )
        {
            case 0:
                    return cos( rIndex * acos( 2 * rhoh - 1 ) ) * cos( thetaIndex * 2 * M_PI * thetah );
            case 1:
                    return cos( rIndex * acos( 2 * rhoh - 1 ) ) * sin( thetaIndex * 2 * M_PI * thetah );
            default:
                    return 0;
        }
    }*/
} //End LifeV namespace

