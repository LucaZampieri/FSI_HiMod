#include <lifev/himod/basis/Bessel.hpp>
#include <cmath>

namespace LifeV
{

void Bessel::
computeEigenvalues()
{

    // besselRoots.txt contains the first 110 roots of the first 101 besselj functions, s.t.
    // besselRoots(i,j) is the i-th root of the besselj function of order j.
/*  std::ifstream besselRoots( "besselRoots.txt" );
    MBMatrix_type roots;
    std::vector<Real> allRoots;
    std::set< EigenMap2D, Comparison2D > subeigen;
    Real  j0(0), j1(0), y0(0), y1(0), j0p(0), j1p(0), y0p(0), y1p(0);
    int  nm;

    M_eigenvalues.resize( M_mtot );
    M_normRJ.resize( M_mtot );
    
    std::copy( std::istream_iterator<Real>( besselRoots ),
               std::istream_iterator<Real>(),
               std::back_inserter<std::vector<Real> >( allRoots ) );

    roots.resize( M_mtot + 1 );
    for( UInt i = 0; i != M_mtot + 1; ++i )
    {
        roots[i].resize( M_mtot + 1 );
    }

    for( UInt i = 0; i != M_mtot + 1; ++i )
    {
        for( UInt j = 0; j != M_mtot + 1; ++j )
        {
            roots[i][j] = allRoots[i*101+j];
        }
    }
*/    

    MBMatrix_type roots, droots;
    std::set< EigenMap2D, Comparison2D > subeigen;
    Real  j0(0), j1(0), y0(0), y1(0), j0p(0), j1p(0), y0p(0), y1p(0);
    int  nm;
    
    M_eigenvalues.resize( M_mtot );
    M_normRJ.resize( M_mtot );

    bessel::besselzeros( 2, 2, roots, droots );
    M_eigenvalues[0].lambda = roots[0][0];
    M_eigenvalues[0].order = 0;
    M_eigenvalues[0].index = 0;

    bessel::bessjy01b( M_eigenvalues[0].lambda, j0, j1, y0, y1, j0p, j1p, y0p, y1p );

    M_normRJ[0] = sqrt(2) / fabs( j1 );

    // Insert the next root for the same order
    subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index + 1][M_eigenvalues[0].order],
                     M_eigenvalues[0].index + 1, M_eigenvalues[0].order ) );
    // Insert the first root of the next order
    subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index][M_eigenvalues[0].order + 1],
                     M_eigenvalues[0].index, M_eigenvalues[0].order + 1 ) );
        
    deleteBesselRoots( roots, droots );

    for( UInt j( 1 ); j != M_mtot; ++j )
    {
        // extract the least eigenvalue among the following ones
        M_eigenvalues[j].lambda = subeigen.begin()->lambda;
        M_eigenvalues[j].order  = subeigen.begin()->order;
        M_eigenvalues[j].index  = subeigen.begin()->index;
                    
        // erase the extract element        
        subeigen.erase( subeigen.begin() );

        // compute the norm
        if( M_eigenvalues[j].order == 0 )
        {
            bessel::bessjy01b( M_eigenvalues[j].lambda, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = sqrt(2) / fabs( j1 ); 
        }  
        else
        {
            UInt dim = M_eigenvalues[j].order+1;
            Real* jn = new Real[dim+1];
            Real* yn = new Real[dim+1];
            Real* jnp = new Real[dim+1];
            Real* ynp = new Real[dim+1];

            bessel::bessjyna( M_eigenvalues[j].order+1, M_eigenvalues[j].lambda, nm, jn, yn, jnp, ynp );
            assert( nm == M_eigenvalues[j].order+1 );
            M_normRJ[j] = sqrt( 2 ) / fabs( jn[nm] );

            delete[] jn;
            delete[] yn;
            delete[] jnp;
            delete[] ynp;
         }

         // the first increment is because of the index of the matrix, the further one because
         // I need one-past the current eigenvalue rightward and bottomward
         bessel::besselzeros( M_eigenvalues[j].order + 2, M_eigenvalues[j].index + 2, roots, droots );

         // insert next eigenvalue, same order
         subeigen.insert( EigenMap2D::make_eigenmap2D( roots[(M_eigenvalues[j]).index + 1][(M_eigenvalues[j]).order],
                          M_eigenvalues[j].index + 1, M_eigenvalues[j].order  ) );

         // insert eigenvalue same index, next order
         subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[j].index][M_eigenvalues[j].order + 1],
                          M_eigenvalues[j].index, M_eigenvalues[j].order + 1 ) );

         deleteBesselRoots( roots, droots );

    }
    return;
}

    void
    Bessel::
    evaluateBasis( MBMatrix_type& phirho,
                   MBMatrix_type& dphirho,
                   MBMatrix_type& phitheta,
                   MBMatrix_type& dphitheta,
                   const QuadratureRule* quadrulerho,     
                   const QuadratureRule* quadruletheta )
    {
        const int trigonometricBasis(0);
        phirho.resize( M_mtot );
        dphirho.resize( M_mtot );
        phitheta.resize( M_mtot );
        dphitheta.resize( M_mtot );
        computeEigenvalues();
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
        
        // Evaluate phirho, dphirho da Educated basis
        Real  j0, j1, y0, y1, j0p, j1p, y0p, y1p;
        
        Real    x, jQuadraturePoint, djQuadraturePoint, yQuadraturePoint, dyQuadraturePoint;
     
        for ( UInt m( 0 ); m != M_mtot; ++m )
        {
            for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                x = M_eigenvalues[m].lambda * quadrulerho->quadPointCoor( n, 0 );

                if( M_eigenvalues[m].order == 0 )
                {
                    bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
                    jQuadraturePoint = j0;
                    djQuadraturePoint = j0p;
                    yQuadraturePoint = y0;
                    dyQuadraturePoint = y0p;
                }
                else if( M_eigenvalues[m].order == 1 )
                {
                    bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
                    jQuadraturePoint = j1;
                    djQuadraturePoint = j1p;
                    yQuadraturePoint = y1;
                    dyQuadraturePoint = y1p;
                }
                else
                {
                    int nm;
                    UInt dim = M_eigenvalues[m].order;
                    Real* jn = new Real[dim+1];
                    Real* yn = new Real[dim+1];
                    Real* jnp = new Real[dim+1];
                    Real* ynp = new Real[dim+1];
    
                    bessel::bessjyna( M_eigenvalues[m].order, x, nm, jn, yn, jnp, ynp );
                    assert( nm == M_eigenvalues[m].order );
                
                    jQuadraturePoint = jn[nm];
                    djQuadraturePoint = jnp[nm];
                    yQuadraturePoint = yn[nm];
                    dyQuadraturePoint = ynp[nm];
                 
                    delete[] jn;
                    delete[] yn;
                    delete[] jnp;
                    delete[] ynp;
                }

                phirho[m][n]  = M_normRJ[m] * jQuadraturePoint;
                dphirho[m][n] = M_eigenvalues[m].lambda * M_normRJ[m] * djQuadraturePoint;

            } // end n-for
        } // end m-for

        Real normTheta( 1. / sqrt( 2*M_PI ) );
        switch( trigonometricBasis )
        {
            case 0:
                    normTheta = 1. / sqrt( 2*M_PI );
                    for ( UInt m = 0; m != M_mtot; ++m )
                    {
                        UInt pt = M_eigenvalues[m].order;
                        normTheta = ( pt==0?  1. / sqrt( 2*M_PI ): 1. / sqrt( M_PI ) );
                        for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                        {
                            phitheta[m][n]  = normTheta *
                                              ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                            dphitheta[m][n] = normTheta * pt *
                                              ( - sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                        } // end n-for
                     } // end m-for
            
                     break;
            
            case 1:            

                    normTheta = 1. / sqrt( M_PI );
                    for ( UInt m = 0; m != M_mtot; ++m )
                    {
                        UInt pt = M_eigenvalues[m].order;
                        for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                        {
                            phitheta[m][n]  = normTheta *
                                              ( sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                            dphitheta[m][n] = normTheta * pt *
                                              ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                        } // end n-for
                    } // end m-for
                
                    break;
                        
               default:
                       break;
           } // end switch
        
        return;
    }

} //End LifeV namespace

