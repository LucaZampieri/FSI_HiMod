#include <lifev/himod/basis/ChebyshevQuadraticShiftNatural.hpp>
#include <cmath>

namespace LifeV
{
    void
    ChebyshevQuadraticShiftNatural::
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta )
        {
            const int trigonometricBasis(0);
            UInt mr1D( 0 ), mtheta1D( 0 );
            rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
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
            MBVector_type normRr( mr1D, 0 );
        
            UInt i( 0 );
            for( UInt pr = 0; pr != mr1D; ++pr )
            {
                i = pr;
                for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                {
                    Real x = 2 * quadrulerho->quadPointCoor( n, 0 ) * quadrulerho->quadPointCoor( n, 0 ) - 1;
                    normR[pr] += ( cos( i * acos( x ) ) * cos( i * acos( x ) ) 
                                 ) *
                                 M_Rho * quadrulerho->quadPointCoor( n, 0 ) * 
                                 M_Rho * quadrulerho->weight( n );
                    normRr[pr] += ( quadrulerho->quadPointCoor( n, 0 ) * quadrulerho->quadPointCoor( n, 0 ) *
                                    cos( i * acos( x ) ) * cos( i * acos( x ) ) 
                                  ) * 
                                  M_Rho * quadrulerho->quadPointCoor( n, 0 ) * 
                                  M_Rho * quadrulerho->weight( n );
                }
            } // end norm computation
        
            for ( UInt pt = 0; pt != mtheta1D; ++pt )
            {
                if( remainder( pt, 2 ) == 0 )
                {
                    for ( UInt pr = 0; pr != mr1D; ++pr )
                    {
                        i = pr;
                        for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                        {
                            Real x = 2 * quadrulerho->quadPointCoor( n, 0 ) * quadrulerho->quadPointCoor( n, 0 ) - 1;
                
                            // NOTA BENE: LE FUNZIONI SONO MOLTIPLICATE PER UN M_Rho IN PIÙ PER COMPENSARE IL normrho
                            //CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                            phirho[mr1D*pt + pr][n] = M_Rho * (
                                                      cos( i * acos( x ) ) 
                                                      ) / sqrt( normR[pr] );
                                                    
                            dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                                       i * sin( i * acos( x ) ) /
                                                       sqrt( 1 - x * x ) *
                                                       4 * quadrulerho->quadPointCoor( n, 0 )
                                                       ) / sqrt( normR[pr] );
                        }
                    }
                } // endif( pt is even )
                else
                {
                    for ( UInt pr = 0; pr != mr1D; ++pr )
                    {
                        i = pr;
                        for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                        {
                            Real x = 2 * quadrulerho->quadPointCoor( n, 0 ) * quadrulerho->quadPointCoor( n, 0 ) - 1;
                            
                            phirho[mr1D*pt + pr][n] = M_Rho * (
                                                      quadrulerho->quadPointCoor( n, 0 ) *
                                                      cos( i * acos( x ) ) 
                                                      )/ sqrt( normRr[pr] );
                                                    
                            dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                                       quadrulerho->quadPointCoor( n, 0 ) *
                                                       i * sin( i * acos( x ) ) /
                                                       sqrt( 1 - x * x )  *
                                                       4 * quadrulerho->quadPointCoor( n, 0 ) +

                                                       cos( i * acos( x ) ) 
                                                       ) / sqrt( normRr[pr] );
                            
                        }// end n-for
                    }// end pr-for
                }// endif( pt is odd )
            } // end pt-for
    
    gramSchmidtOrthogonalization( phirho, dphirho,
                                      quadrulerho,
                                      mr1D,
                                      mtheta1D );
                                        
    checkOrthogonality( phirho, quadrulerho, mr1D, mtheta1D, 0 );
    
    switch( trigonometricBasis )
    {
            case 0:
                      for ( UInt pt = 0; pt != mtheta1D; ++pt )
                      {
                          Real correction( 1 );
                          if( pt==0 )
                          {
                              correction = 1. / sqrt( 2 );
                          }
    
                          for ( UInt pr = 0; pr != mr1D; ++pr )
                          {
                              for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                              {
                                  // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE 
                                  //NEGLI INTEGRALI DI ModalSpace e HMA.
                                  // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
                                  phitheta[mr1D*pt + pr][n]  = M_normTheta * correction *
                                                               ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                  dphitheta[mr1D*pt + pr][n] = M_normTheta * correction * pt *
                                                               ( - sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                               } // end n-for
                           } // end pr-for
                       } // end pt-for
                
                       break;
            
                case 1:            

                       for ( UInt pt = 0; pt != mtheta1D; ++pt )
                       {
                           for ( UInt pr = 0; pr != mr1D; ++pr )
                           {
                               for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                               {
                                   phitheta[mr1D*pt + pr][n]  = M_normTheta *
                                                                ( sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                   dphitheta[mr1D*pt + pr][n] = M_normTheta * pt *
                                                                ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                               } // end n-for
                            } // end pr-for
                        }

                        break;
                        
        default:
                        break;
        } // end switch
        
        return;
    }
    
    /*Real
    ChebyshevQuadraticShiftNatural::
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
                    if( remainder( mtheta1D, 2 )==0 )
                    {
                        return cos( rIndex * acos( 2 * rhoh * rhoh - 1 ) ) * cos( thetaIndex * 2 * M_PI * thetah );
                    }
                    else
                    {
                        return rhoh * cos( rIndex * acos( 2 * rhoh * rhoh - 1 ) ) * cos( thetaIndex * 2 * M_PI * thetah );
                    }
            case 1:
                    if( remainder( mtheta1D, 2 )==0 )
                    {
                        return cos( rIndex * acos( 2 * rhoh * rhoh - 1 ) ) * sin( thetaIndex * 2 * M_PI * thetah );
                    }
                    else
                    {
                        return rhoh * cos( rIndex * acos( 2 * rhoh * rhoh - 1 ) ) * sin( thetaIndex * 2 * M_PI * thetah );
                    }
            default:
                    return 0;
        }
    }*/
} //End LifeV namespace

