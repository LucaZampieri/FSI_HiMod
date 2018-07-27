#include <lifev/himod/basis/ChebyshevQuadraticShift.hpp>
#include <cmath>

namespace LifeV
{
    Real
    ChebyshevQuadraticShift::
    radialFunction( const UInt& pr, const UInt& pt, const Real& rStar,
                   const QuadratureRule* quadrulerho ) const 
    {

        Real norm(0);
        Real xStar( 2 * rStar * rStar - 1);
        Real phirho(0);

        UInt i( pr + 1 );
        if( remainder( pt, 2 ) == 0 )
        {
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                norm += ( T(i,x) * T(i,x) -
                           2 * T(i,x) +
                           1
                         ) *
                         M_Rho * r * 
                         M_Rho * dr;
            }

            phirho = M_Rho * (
                     T(i,xStar) - 1
                     ) / sqrt( norm );


        }
        else
        {
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                norm += ( r * r *
                            T(i,x) * T(i,x) -
                            2 * r * r *
                            T(i,x)  +
                            r * r
                          ) * 
                          M_Rho * r * 
                          M_Rho * dr;
            }

            phirho = M_Rho * (
                     rStar * T(i,xStar) - rStar
                     )/ sqrt( norm );
        }


       return phirho; 
    }

    Real
    ChebyshevQuadraticShift::
    angularFunction( const UInt& pt, const Real& tStar ) const
    {
        Real phitheta( 0 );

        int trigFunction(0);
        switch( trigFunction )
        {
            case 0: {
                    Real correction( 1 );
                    if( pt==0 )
                    {
                        correction = 1. / sqrt( 2 );
                    }
                                // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE NEGLI INTEGRALI
                                // DI ModalSpace e HMA.
                                // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
                    phitheta = M_normTheta * correction *
                               ( cos ( pt * tStar ) );

                    }break;
            
            case 1: {           
                    phitheta = M_normTheta *
                               ( sin ( pt * tStar ) );
                    
                    }break;
                        
            default:
                    break;
        }

        return phitheta; 
    }

    Real 
    ChebyshevQuadraticShift::
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar,
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        UInt pt( k/mr1D );
        UInt pr( k - mr1D*pt );

        Real phir( xi( pr+1, pt, rStar, quadrulerho, phirho ) );
        Real phit( angularFunction( pt, tStar ) );
        return phir*phit;
    }

    Real
    ChebyshevQuadraticShift::
    xi( const UInt& kr, const UInt& kt, const Real& rStar, const QuadratureRule* quadrulerho, const MBMatrix_type& phirho ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        std::vector<Real> B( mr1D*kt+kr, 0 );
        std::vector<Real> Bquad( quadrulerho->nbQuadPt(), 0 );

        for( UInt pt(0); pt<=kt; ++pt )
        {
            UInt mr( pt==kt ? kr : mr1D );
            for( UInt h(0); h<mr; ++h )
                B[pt*mr1D+h] = radialFunction( h, pt, rStar, quadrulerho );
        }

        for( UInt n(0); n<=quadrulerho->nbQuadPt(); ++n )
            Bquad[n] = radialFunction( kr, kt, quadrulerho->quadPointCoor(n,0), quadrulerho );

        Real projection( 0 );
        //Real norm( 0 );
    
        for( UInt p( 0 ); p <= kt; ++p )
        {
            UInt mr1d( p==kt ? kr : mr1D );
            for( UInt k( 1 ); k != mr1d; ++k )
            {
                for( UInt pr( 0 ); pr != k; ++pr )
                {
                    // compute the projection
                    for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                    {
                        projection += Bquad[n] *
                                      phirho[mr1D * p + pr][n] *
                                      quadrulerho->quadPointCoor( n, 0 ) *
                                      quadrulerho->weight( n );
                    }

                    // compute orthogonal function (not normalized yet)
                    B[mr1D*p+k] -= B[mr1D*p+pr]*projection;
                    
                    projection = 0;
                } // end pr-for
                
                // This should be unnecessary because phirho is already normalized.
                // compute norm
                /* for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    norm += phirho[mr1D*p+k][n] * phirho[mr1D*p+k][n] *
                            quadrulerho->quadPointCoor( n, 0 ) *
                            quadrulerho->weight( n );
                }
                */    
                // normalize
                // B[mr1D*p+k] /= sqrt( norm );
                
            } // end k-for
        } // end p-for

        return B[ mr1D*kt+kr ];

    }


    void
    ChebyshevQuadraticShift::
    evaluateBasis( MBMatrix_type& phirho,
                   MBMatrix_type& dphirho,
                   MBMatrix_type& phitheta,
                   MBMatrix_type& dphitheta,
                   const QuadratureRule* quadrulerho,     
                   const QuadratureRule* quadruletheta )
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        M_eigenvalues.resize( M_mtot );
        phirho.resize( M_mtot );
        dphirho.resize( M_mtot );
        phitheta.resize( M_mtot );
        dphitheta.resize( M_mtot );
        
        /*boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
        MapEpetra Map( M_mtot, Comm );

        //Using the map it is possible to define the system matrix
        M_massMatrix = static_cast< boost::shared_ptr<matrix_Type>  > ( new matrix_Type( Map ) );

        std::vector<UInt> block_row( 1, M_mtot ); // one block with mtot cells
        std::vector<UInt> block_col( 1, M_mtot );

        M_massMatrix->setBlockStructure( block_row, block_col );
        *M_massMatrix *= 0.0;
        */
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
            i = pr + 1;
            for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                Real r( quadrulerho->quadPointCoor( n, 0 ) );
                Real dr( quadrulerho->weight( n ) );
                Real x = 2 * r * r - 1;

                normR[pr] += ( T(i,x) * T(i,x) -
                
                               2 * T(i,x) +
                            
                               1
                             ) *
                             M_Rho * r * 
                             M_Rho * dr;
                normRr[pr] += ( r * r *
                                T(i,x) * T(i,x) -
                                
                                2 * r * r *
                                T(i,x)  +
                                
                                r * r
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
                    M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                    M_eigenvalues[pt*mr1D+pr].order = pt;
                    M_eigenvalues[pt*mr1D+pr].index = pr;
                    
                    i = pr + 1;
                    for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                    {
                        Real r( quadrulerho->quadPointCoor( n, 0 ) );
                        Real x = 2 * r * r - 1;
                            
                        // NOTA BENE: LE FUNZIONI SONO MOLTIPLICATE PER UN M_Rho IN PIÙ PER COMPENSARE IL normrho
                        //CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                        phirho[mr1D*pt + pr][n] = M_Rho * (
                                                  T(i,x) - 1
                                                  ) / sqrt( normR[pr] );
                                                    
                        dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                                   dT(i,x) *
                                                   4 * r
                                                   ) / sqrt( normR[pr] );
                    }
                }
            } // endif( pt is even )
            else
            {
                for ( UInt pr = 0; pr != mr1D; ++pr )
                {
                    M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                    M_eigenvalues[pt*mr1D+pr].order = pt;
                    M_eigenvalues[pt*mr1D+pr].index = pr;
                    
                    i = pr + 1;
                    for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                    {
                        Real r( quadrulerho->quadPointCoor( n, 0 ) );
                        Real x = 2 * r * r - 1;
                            
                        phirho[mr1D*pt + pr][n] = M_Rho * (
                                                  r * T(i,x) - r
                                                  )/ sqrt( normRr[pr] );
                                                    
                        dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                                   r *
                                                   dT(i,x)  *
                                                   4 * r +

                                                   T(i,x) -
                                                    
                                                   1
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
    
        int trigonometricBasis(0);
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
                                // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE NEGLI INTEGRALI
                                // DI ModalSpace e HMA.
                                // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
                                phitheta[mr1D*pt + pr][n]  = M_normTheta * correction *
                                                             ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                dphitheta[mr1D*pt + pr][n] = M_normTheta * correction * pt *
                                                             ( - sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                            } // end n-for
                        } // end pr-for
                    }
                
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

} //End LifeV namespace

