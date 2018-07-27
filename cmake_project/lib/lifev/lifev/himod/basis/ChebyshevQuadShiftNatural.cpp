#include <lifev/himod/basis/ChebyshevQuadShiftNatural.hpp>
#include <cmath>

namespace LifeV
{
    Real
    ChebyshevQuadShiftNatural::
    radialFunction( const UInt& pr, const UInt& pt, const Real& rStar,
                   const QuadratureRule* quadrulerho ) const 
    {

        Real norm(0);
        Real xStar( 2 * rStar * rStar - 1);
        Real phirho(0);

        UInt i( pr );
        if( remainder( pt, 2 ) == 0 )
        {
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                norm += ( T(i,x) * T(i,x)
                         ) *
                         M_Rho * r * 
                         M_Rho * dr;
            }

            phirho = T(i,xStar) / sqrt( norm );

        }
        else
        {
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                norm += ( r * r *
                          T(i,x) * T(i,x)
                          ) * 
                          M_Rho * r * 
                          M_Rho * dr;
            }

            phirho = (
                     rStar * T(i,xStar)
                     )/ sqrt( norm );
        }


       return phirho; 
    }

    Real
    ChebyshevQuadShiftNatural::
    angularFunction( const UInt& pt, const Real& tStar, const UInt& trigFunction ) const
    {
        Real phitheta( 0 );

        Real normTheta( 1. / sqrt( M_PI ) );
        switch( trigFunction )
        {
            case 0: {
                    if( pt==0 )
                    {
                        normTheta = 1. / sqrt( 2*M_PI );
                    }
                    else
                    {
                        normTheta = 1. / sqrt( M_PI );
                    }
                    phitheta = normTheta *
                               ( cos ( pt * tStar ) );

                    }break;
            
            case 1: {           
                    phitheta = normTheta *
                               ( sin ( pt * tStar ) );
                    
                    }break;
                        
            default:
                    break;
        }

        return phitheta; 
    }

    Real 
    ChebyshevQuadShiftNatural::
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar,
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 ), offset( M_mtot/2 );
        rectangularTruncation( M_mtot/2, mr1D, mtheta1D );

        UInt kInd( k<offset? k : k-offset );
        UInt trigFun( k<offset? 0 : 1 );
        UInt pt( kInd/mr1D );
        UInt pr( kInd - mr1D*pt );

        Real phir( xi( pr, pt, rStar, quadrulerho, phirho ) );
        Real phit( angularFunction( pt, tStar, trigFun ) );
        return phir*phit;
    }

    Real
    ChebyshevQuadShiftNatural::
    xi( const UInt& kr, const UInt& kt, const Real& rStar, const QuadratureRule* quadrulerho, const MBMatrix_type& phirho ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot/2, mr1D, mtheta1D );
        std::vector<Real> B( kr+1, 0 );
        std::vector<Real> Bquad( quadrulerho->nbQuadPt(), 0 );
        std::vector<std::vector<Real> > Bkquad( kr+1, Bquad );

        for( UInt h(0); h<kr+1; ++h )
        {
            B[h] = radialFunction( h, kt, rStar, quadrulerho );
            // For the normalization
            for( UInt n(0); n<quadrulerho->nbQuadPt(); ++n )
                Bkquad[h][n] = radialFunction( h, kt, quadrulerho->quadPointCoor(n,0), quadrulerho );
        }

        Real projection( 0 );
        Real norm( 0 );
    
        for( UInt k( 1 ); k < kr+1; ++k )
        {
            for( UInt pr( 0 ); pr < k; ++pr )
            {
                    // compute the projection
                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        projection += Bkquad[k][n] *
                                      Bkquad[pr][n] *
                                      //phirho[mr1D*kt + pr][n] *
                                      quadrulerho->quadPointCoor( n, 0 ) *
                                      quadrulerho->weight( n );
                    }

                    // compute orthogonal function (not normalized yet)
                    B[k] -= B[pr]*projection;

                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        Bkquad[k][n] -= Bkquad[pr][n]*projection;
                    }
                    
                    projection = 0;
                } // end pr-for
                
                // compute norm
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    norm += Bkquad[k][n] * Bkquad[k][n] *
                            quadrulerho->quadPointCoor( n, 0 ) *
                            quadrulerho->weight( n );
                }
                    
                // normalize
                B[k] /= sqrt( norm );
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    Bkquad[k][n] /= sqrt(norm);
                }
                norm = 0;
                
            } // end k-for

        return B[ kr ];

    }

   void
    ChebyshevQuadShiftNatural::
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta )
    {
        assert( M_mtot%2==0 && "The number of modes must be even!" );
        UInt mr1D( 0 ), mtheta1D( 0 ), offset( M_mtot/2 );
        rectangularTruncation( M_mtot/2, mr1D, mtheta1D );
    
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
        MBVector_type normRr( mr1D, 0 );

        UInt i( 0 );
        for( UInt pr = 0; pr != mr1D; ++pr )
        {
            i = pr;
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                normR[pr] += ( T(i,x) * T(i,x) 
                             ) *
                             M_Rho * r * 
                             M_Rho * dr;
                normRr[pr] += ( r * r *
                                T(i,x) * T(i,x) 
                              ) * 
                              M_Rho * r * 
                              M_Rho * dr;
            }
        } // end norm computation
    
        for ( UInt pt = 0; pt != mtheta1D; ++pt )
        {
            if( remainder( pt, 2 ) == 0 )
            {
                for ( UInt pr = 0; pr != mr1D; ++pr )
                {
                    i = pr;

                    M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                    M_eigenvalues[offset+pt*mr1D+pr].lambda = pt*mr1D+pr;
                    //M_eigenvalues[pt*mr1D+pr].order = pt;
                    M_eigenvalues[pt*mr1D+pr].index = pr;
                    M_eigenvalues[offset+pt*mr1D+pr].index = pr;

                    M_eigenvalues[pt*mr1D+pr].order = pt; 
                    M_eigenvalues[offset+pt*mr1D+pr].order = pt+1; 

                    for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                    {
                        Real r( quadrulerho->quadPointCoor( n, 0 ) );
                        Real x = 2 * r * r - 1;

                        phirho[mr1D*pt + pr][n] = T(i,x) / sqrt( normR[pr] );
                        //phirho[offset+mr1D*pt + pr][n] = phirho[mr1D*pt + pr][n];
                        phirho[offset+mr1D*pt + pr][n] = (
                                                  r * T(i,x)
                                                  )/ sqrt( normRr[pr] );
                                                
                        dphirho[mr1D*pt + pr][n] = ( dT(i,x) *
                                                 4 * r
                                                   ) / sqrt( normR[pr] );
                        //dphirho[offset+mr1D*pt + pr][n] = dphirho[mr1D*pt + pr][n];
                        dphirho[offset+mr1D*pt + pr][n] = ( 
                                                   r *
                                                   dT(i,x)  *
                                                   4 * r +

                                                   T(i,x) 
                                                   ) / sqrt( normRr[pr] );
                    }
                }
            } // endif( pt is even )
            else
            {
                for ( UInt pr = 0; pr != mr1D; ++pr )
                {
                    M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                    M_eigenvalues[offset+pt*mr1D+pr].lambda = pt*mr1D+pr;

                    M_eigenvalues[pt*mr1D+pr].order = pt;
                    M_eigenvalues[offset+pt*mr1D+pr].order = pt+1; 

                    M_eigenvalues[pt*mr1D+pr].index = pr;
                    M_eigenvalues[offset+pt*mr1D+pr].index = pr;

                    i = pr;
                    for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                    {
                        Real r( quadrulerho->quadPointCoor( n, 0 ) );
                        Real x = 2 * r * r - 1;
                        
                        phirho[mr1D*pt + pr][n] = (
                                                  r * T(i,x)
                                                  )/ sqrt( normRr[pr] );
                        //phirho[offset + mr1D*pt + pr][n] = phirho[mr1D*pt + pr][n];
                        phirho[offset+mr1D*pt + pr][n] = T(i,x) / sqrt( normR[pr] );
                                                
                        dphirho[mr1D*pt + pr][n] = ( 
                                                   r *
                                                   dT(i,x)  *
                                                   4 * r +

                                                   T(i,x) 
                                                   ) / sqrt( normRr[pr] );
                        //dphirho[offset + mr1D*pt + pr][n] = dphirho[mr1D*pt + pr][n];
                        dphirho[offset+mr1D*pt + pr][n] = ( dT(i,x) *
                                                 4 * r
                                                   ) / sqrt( normR[pr] );
                        
                    }// end n-for
                }// end pr-for
            }// endif( pt is odd )
        } // end pt-for
    
        gramSchmidtOrthogonalization( phirho, dphirho,
                                      quadrulerho,
                                      mr1D,
                                      mtheta1D );
                                        
        //checkOrthogonality( phirho, quadrulerho, mr1D, mtheta1D, 0 );
    
        Real normTheta( 1. / sqrt( 2*M_PI ) );
        Real theta( 0 );

        for ( UInt pt = 0; pt != mtheta1D; ++pt )
        {
             UInt ptSin( pt+1 );
            for ( UInt pr = 0; pr != mr1D; ++pr )
            {
                for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                {
                    theta = 2 * M_PI * quadruletheta->quadPointCoor( n, 0 );

                    phitheta[mr1D*pt + pr][n]  = normTheta *
                                                 ( cos ( pt * theta ) );
                    phitheta[offset + mr1D*pt + pr][n]  = normTheta *
                                                  ( sin ( ptSin * theta ) );

                    dphitheta[mr1D*pt + pr][n] = normTheta * pt *
                                                 ( - sin ( pt * theta ) );
                    dphitheta[offset+ mr1D*pt + pr][n] = normTheta * ptSin *
                                                  ( cos ( ptSin * theta ) );
                 } // end n-for
             } // end pr-for
               
            // adjust the norm after pt=0
            // The norm of cos(0), indeed, is sqrt(2*pi), while the one of sin or cos is sqrt(pi).
            normTheta = 1. / sqrt( M_PI );
         } // end pt-for
                
        return;
    }
    
} //End LifeV namespace

