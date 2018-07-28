#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include <lifev/himod/basis/EducatedBasisD.hpp>
#include <lifev/navier_stokes/function/bessel/bessel.hpp>
#include <fstream>

namespace LifeV
{

void EducatedBasisD::
computeEigenvalues()
{

		// besselRoots.txt contains the first 110 roots of the first 101 besselj functions, s.t.
		// besselRoots(i,j) is the i-th root of the besselj function of order j.
/*		std::ifstream besselRoots( "besselRoots.txt" );
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
	
		M_eigenvalues[0].lambda = roots[0][0] * roots[0][0] / ( M_Rho * M_Rho );
		M_eigenvalues[0].order = 0;
		M_eigenvalues[0].index = 0;

		bessel::bessjy01b( sqrt( M_eigenvalues[0].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );

		M_normRJ[0] = 1. / ( 1. / sqrt( 2 ) * fabs( j0p ) );

		// Insert the next root for the same order
		subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index + 1][M_eigenvalues[0].order],
														M_eigenvalues[0].index + 1, M_eigenvalues[0].order ) );
		// Insert the first root of the next order
		subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index][M_eigenvalues[0].order + 1],
														M_eigenvalues[0].index, M_eigenvalues[0].order + 1 ) );
		
		deleteBesselRoots( roots, droots );

		for( int j = 1; j != M_mtot; ++j )
		{
			// extract the least eigenvalue among the following ones
		    M_eigenvalues[j].lambda = ( subeigen.begin()->lambda ) * ( subeigen.begin()->lambda ) / ( M_Rho * M_Rho );
		    M_eigenvalues[j].order  = subeigen.begin()->order;
		    M_eigenvalues[j].index  = subeigen.begin()->index;
		            
		    // erase the extract element        
		    subeigen.erase( subeigen.begin() );

		    // compute the norm
		    if( M_eigenvalues[j].order == 0 )
			{
		        bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
		        M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( j0p ) ); 
		    }  
			else if( M_eigenvalues[j].order == 1 )
		    {
		        bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
		        M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( j1p ) );
		    } 
			else
		    {
		    	UInt dim = M_eigenvalues[j].order;
		   		Real* jn = new Real[dim+1];
		   		Real* yn = new Real[dim+1];
		   		Real* jnp = new Real[dim+1];
		   		Real* ynp = new Real[dim+1];

		        bessel::bessjyna( M_eigenvalues[j].order, sqrt( M_eigenvalues[j].lambda ) * M_Rho, nm, jn, yn, jnp, ynp );
		        assert( nm == M_eigenvalues[j].order );
		        M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( jnp[nm] ) );

		        delete[] jn;
		        delete[] yn;
		        delete[] jnp;
		        delete[] ynp;
		    }

		    // the first increment is because of the index of the matrix, the further one because I need one-past the current eigenvalue rightward and bottomward
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

/*void EducatedBasisD::
computeEigenvalues()
{
	MBMatrix_type roots, droots;
    std::set< EigenMap2D, Comparison2D > subeigen;
    Real  j0(0), j1(0), y0(0), y1(0), j0p(0), j1p(0), y0p(0), y1p(0);
    int  nm;
	
	M_eigenvalues.resize( M_mtot );
	M_normRJ.resize( M_mtot );

	bessel::besselzeros( 2, 2, roots, droots );
	
	M_eigenvalues[0].lambda = roots[0][0] * roots[0][0] / ( M_Rho * M_Rho );
	M_eigenvalues[0].order = 0;
	M_eigenvalues[0].index = 0;

	bessel::bessjy01b( sqrt( M_eigenvalues[0].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );

	M_normRJ[0] = 1. / ( 1. / sqrt( 2 ) * fabs( j0p ) );

	// Insert the next root for the same order
    subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index + 1][M_eigenvalues[0].order], M_eigenvalues[0].index + 1, M_eigenvalues[0].order ) );
    // Insert the first root of the next order
    subeigen.insert( EigenMap2D::make_eigenmap2D( roots[M_eigenvalues[0].index][M_eigenvalues[0].order + 1], M_eigenvalues[0].index, M_eigenvalues[0].order + 1 ) );
    
	deleteBesselRoots( roots, droots );

    for( int j = 1; j != M_mtot; ++j )
    {
    	// extract the least eigenvalue among the following ones
        M_eigenvalues[j].lambda = ( subeigen.begin()->lambda ) * ( subeigen.begin()->lambda ) / ( M_Rho * M_Rho );
        M_eigenvalues[j].order  = subeigen.begin()->order;
        M_eigenvalues[j].index  = subeigen.begin()->index;
                
        // erase the extract element        
        subeigen.erase( subeigen.begin() );

        // compute the norm
        if( M_eigenvalues[j].order == 0 )
	    {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( j0p ) ); 
        }  
	    else if( M_eigenvalues[j].order == 1 )
        {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( j1p ) );
        } 
	    else
        {
        	UInt dim = M_eigenvalues[j].order;
       		Real* jn = new Real[dim+1];
       		Real* yn = new Real[dim+1];
       		Real* jnp = new Real[dim+1];
       		Real* ynp = new Real[dim+1];

            bessel::bessjyna( M_eigenvalues[j].order, sqrt( M_eigenvalues[j].lambda ) * M_Rho, nm, jn, yn, jnp, ynp );
            assert( nm == M_eigenvalues[j].order );
            M_normRJ[j] = 1. / ( 1. / sqrt( 2 ) * fabs( jnp[nm] ) );

            delete[] jn;
            delete[] yn;
            delete[] jnp;
            delete[] ynp;
        }

        // the first increment is because of the index of the matrix, the further one because I need one-past the current eigenvalue rightward and bottomward
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
*/

}// end namespace lifev
