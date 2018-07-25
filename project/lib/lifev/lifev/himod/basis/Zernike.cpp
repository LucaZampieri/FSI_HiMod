#include <lifev/himod/basis/Zernike.hpp>
#include <cmath>

namespace LifeV
{

    void
    Zernike::
    triangularTruncation( const UInt& mtot, UInt& mr1D, UInt& mtheta1D, const bool verbose ) const
    {
        std::vector<UInt> n;
        std::vector<UInt> m;

        for( UInt k(0); k!=mtot; ++k )
        {
           n.push_back( std::ceil( 0.5*(-3+std::sqrt( 9+8*k )) ) );
           m.push_back( std::abs(2*k-n[k]*(n[k]+2)) );
        }

        std::vector<UInt> n_unique( *(std::unique( n.begin(), n.end() )) );
        std::vector<UInt> m_unique( *(std::unique( m.begin(), m.end() )) );
        mr1D = n_unique.size();
        mtheta1D = m_unique.size();
    }
    
    UInt
    Zernike::
    factorial( const UInt& k ) const
    {
        UInt fact( 1 );
        for( UInt i(2); i<=k; ++i )
             fact *= i;

        return fact;
    }

    Real
    Zernike::
    Rmn( const UInt& m_, const UInt& n_, const Real& rStar ) const 
    {
        // UInt --> Real conversions
        Real m( (Real)(m_) );
        Real n( (Real)(n_) );

        Real z(0);

        UInt d( 0.5*(n-m)+1 );

        std::vector<Real> rPowern( d, std::pow(rStar,m) );
        for( UInt i(1); i < d; ++i )
        {
            rPowern[i] = rPowern[i-1]*rStar*rStar;
        }
       
        for( UInt k(0); k<d; ++k  ) 
        {
            //UInt k( kk-1 );
            Real p( ( 1.-2.*(k%2) ) * factorial( n-k ) / 
                    ( factorial( k ) *  
                      factorial( (UInt)(0.5*(n-m)) - k ) *
                      factorial( (UInt)(0.5*(n+m)) - k ) ) );
            z += p*rPowern[d-1-k];
        }

        return z; 
    }

    Real
    Zernike::
    radialFunction( const UInt& m, const UInt& n, const Real& rStar ) const
    {
        Real phirho( Rmn(m,n,rStar)*(1-rStar*rStar)  );
        return phirho;
    }

    Real
    Zernike::
    derRadialFunction( const UInt& m, const UInt& n, const Real& rStar ) const
    {
        Real dphirho( derRmn(m,n,rStar)*(1-rStar*rStar) + Rmn(m,n,rStar)*(-2*rStar)  );
        return dphirho;
    }

    Real
    Zernike::
    derRmn( const UInt& m_, const UInt& n_, const Real& r ) const 
    {
        // UInt --> Real conversions
        Real m( (Real)(m_) );
        Real n( (Real)(n_) );

        Real der( 0. );

        if( m_==0 && n_==0 )
        {
           der = 0.;
        }
        else if( m_==1 && n_==1 )
        {
           der = 1.;
        }
        else
        {
           der = ( ( n+2 )*r*r + m )/( r*(1-r*r) )*Rmn( m_,n_,r ) - (m+n+2)/(1-r*r)*Rmn(m_+1,n_+1,r);
        }

        return der;
    }

    Real
    Zernike::
    angularFunction( const int& m, const Real& tStar ) const
    {
        Real phitheta( 0 );
        Real m_abs( std::fabs( (Real) (m) ) );

        Real normTheta( 1. / sqrt( M_PI ) );
        phitheta = m>=0? std::cos( m_abs*tStar )*normTheta : std::sin( m_abs*tStar )*normTheta;
        phitheta = m==0? phitheta/std::sqrt(2) : phitheta;

        return phitheta; 
    }

    Real
    Zernike::
    derAngularFunction( const int& m, const Real& tStar ) const
    {
        Real dphitheta( 0 );
        Real m_abs( std::fabs( (Real) (m) ) );

        Real normTheta( 1. / sqrt( M_PI ) );
        dphitheta = m>=0? -m_abs*std::sin( m_abs*tStar )*normTheta : m_abs*std::cos( m_abs*tStar )*normTheta;
        dphitheta = m==0? dphitheta/std::sqrt(2) : dphitheta;

        return dphitheta; 
    }


    Real 
    Zernike::
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar,
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const
    {
        UInt n( std::ceil( 0.5*(-3+std::sqrt( 9+8*k )) ) );
        int m( 2*k-n*(n+2) );

        Real phir( radialFunction( std::abs(m), n, rStar ) );
        Real phit( angularFunction( m, tStar ) );
        return phir*phit;
    }

    void
    Zernike::
    evaluateBasis( MBMatrix_type& phirho,
                   MBMatrix_type& dphirho,
                   MBMatrix_type& phitheta,
                   MBMatrix_type& dphitheta,
                   const QuadratureRule* quadrulerho,     
                   const QuadratureRule* quadruletheta )
    {
        UInt mr1D(0);
        UInt mtheta1D(0);
        triangularTruncation( M_mtot, mr1D, mtheta1D ); 

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
            phitheta[i].resize( quadruletheta->nbQuadPt() );
            dphitheta[i].resize( quadruletheta->nbQuadPt() );
        }
        
        // Initialize eigenvalues 
        for( UInt p = 0; p != M_mtot; ++p )
        {
            UInt n( std::ceil( 0.5*(-3+std::sqrt( 9+8*p )) ) );
            int m( 2*p-n*(n+2) );

            M_eigenvalues[p].lambda = p;
            M_eigenvalues[p].index = n;
            M_eigenvalues[p].order = fabs(m);
            
        }
      
        // Compute basis functions
        // Maximum degree of the radial functions (which is also the maximum angular frequency) 
        UInt nMax( std::ceil( 0.5*(-3+std::sqrt( 9+8*(M_mtot-1) )) ) );
        for ( UInt m = 0; m <= nMax; ++m )
        {
            //UInt index( orderedIndices[p] );
            UInt nRmn( 0.5*(nMax-m) + 1 );
            std::vector<std::vector<Real> > xi( nRmn, std::vector<Real> ( quadrulerho->nbQuadPt(), 0 ) );
            std::vector<std::vector<Real> > dxi( nRmn, std::vector<Real> ( quadrulerho->nbQuadPt(), 0 ) );

            MBVector_type normR( nRmn, 0 );
        
            for( UInt i(0); i<nRmn; ++i )
            {
                UInt n( m+2*i );

                for( UInt j = 0; j != quadrulerho->nbQuadPt(); ++j )
                {
                    Real r( quadrulerho->quadPointCoor( j, 0 ) );
                    Real dr( quadrulerho->weight( j ) );
                    Real z( radialFunction( m,n,r ) );

                    normR[i] += ( z * z
                                 ) *
                                 r * 
                                 dr;
                }
                for ( UInt j = 0; j != quadrulerho->nbQuadPt(); ++j )
                {
                    Real r( quadrulerho->quadPointCoor( j, 0 ) );

                    xi[i][j] = radialFunction(m,n,r) / sqrt( normR[i] );
                                                    
                    dxi[i][j] = derRadialFunction(m,n,r) / sqrt( normR[i] );
                }
            }

            orthogonalize( xi, dxi, quadrulerho );

            for( UInt i(0); i<nRmn; ++i )
            {
                UInt n( m+2*i );

                UInt pCos( 0.5*( m+n*(n+2) ) );
                UInt pSin( 0.5*( -m+n*(n+2) ) );

                if( pCos < M_mtot )
                {
                    phirho[pCos] = xi[i];
                    dphirho[pCos] = dxi[i];
                }
                if( pSin < M_mtot )
                {
                    phirho[pSin] = xi[i];
                    dphirho[pSin] = dxi[i];
                }
            }
        } // end m-for
    
        Real theta( 0 );
        for ( UInt p = 0; p != M_mtot; ++p )
        {
             UInt n( std::ceil( 0.5*(-3+std::sqrt( 9+8*p )) ) );
             int m( 2*p-n*(n+2) );

             //UInt index( orderedIndices[p] );
             UInt index( p );

             for ( UInt j = 0; j != quadruletheta->nbQuadPt(); ++j )
             {
                 theta = 2 * M_PI * quadruletheta->quadPointCoor( j, 0 );

                 phitheta[index][j]  = angularFunction(m, theta );
                 dphitheta[index][j] = derAngularFunction(m, theta );
             } // end n-for
             
         }
                
         return;
    }

    void
    Zernike::
    orthogonalize( MBMatrix_type& phirho,
                   MBMatrix_type& dphirho,
                   const QuadratureRule* quadrulerho )
    {
        Real projection( 0 );
        Real norm( 0 );
        UInt M( phirho.size() ); // number of radial functions
    
        for( UInt k( 1 ); k != M; ++k )
        {
            for( UInt pr( 0 ); pr != k; ++pr )
            {
                // compute the projection
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    projection += phirho[k][n] *
                                    phirho[pr][n] *
                                    quadrulerho->quadPointCoor( n, 0 ) *
                                    quadrulerho->weight( n );
                }

                // compute orthogonal function (not normalized yet)
                for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
                {
                    phirho[k][n] -= phirho[pr][n]*projection;
                    dphirho[k][n] -= dphirho[pr][n]*projection;
                }
                projection = 0;
            } // end pr-for
            
            // compute norm
            for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
            {
                norm += phirho[k][n] * phirho[k][n] *
                        quadrulerho->quadPointCoor( n, 0 ) *
                        quadrulerho->weight( n );
            }
                
            // normalize
            for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
            {
                phirho[k][n] /= sqrt( norm );
                dphirho[k][n] /= sqrt( norm );
            }
                
            norm = 0;
        } // end k-for

/*for(UInt i(0); i<M; ++i)
{
    for(UInt j(0); j<M; ++j)
    {
        Real prod( 0 );
            for( UInt n( 0 ); n != quadrulerho->nbQuadPt(); ++n )
            {
                prod += phirho[i][n] * phirho[j][n] *
                        quadrulerho->quadPointCoor( n, 0 ) *
                        quadrulerho->weight( n );
            }
 std::cout << "i=" << i << " j=" << j << " prod=" << prod << std::endl;
    }
}*/
}
 
} //End LifeV namespace

