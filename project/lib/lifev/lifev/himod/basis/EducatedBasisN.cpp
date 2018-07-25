#pragma GCC diagnostic ignored "-Wsign-compare"

#include <lifev/himod/basis/EducatedBasisN.hpp>
#include <iomanip>
namespace LifeV
{

void EducatedBasisN::
computeEigenvalues()
{
	MBMatrix_type roots, droots;
    std::set< EigenMap2D, Comparison2D > subeigen;
    Real  j0(0), j1(0), y0(0), y1(0), j0p(0), j1p(0), y0p(0), y1p(0);
    int  nm;
	
	M_eigenvalues.resize( M_mtot );
	M_normRJ.resize( M_mtot );

	bessel::besselzeros( 2, 2, roots, droots );

    // Compute the first M_mtot roots of the derivative of the first M_mtot orders
/*    for ( UInt i( 0 ); i != M_mtot; ++i )
        for ( UInt j( 0 ); j != M_mtot; ++j )
        {
            subeigen.insert( EigenMap2D::make_eigenmap2D( droots[j][i], j, i ) );
//            std::cout << "i=" << i << ", j=" << j << ", root=" << droots[j][i] << std::endl;
        }
*/
    // Extract first root (lambda=0, order 0) and compute the norm.
    M_eigenvalues[0].lambda = droots[0][0] * droots[0][0] / ( M_Rho * M_Rho );
	M_eigenvalues[0].order = 0;
	M_eigenvalues[0].index = 0;

//    std::cout << "first root: " << M_eigenvalues[0].lambda << std::endl;

	bessel::bessjy01b( sqrt( M_eigenvalues[0].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );

	M_normRJ[0] = 1. / ( 1. / sqrt( 2 ) * fabs( j0 ) );
//    subeigen.erase( subeigen.begin() );

	// Insert the next root for the same order
    subeigen.insert( EigenMap2D::make_eigenmap2D( droots[M_eigenvalues[0].index + 1][M_eigenvalues[0].order], M_eigenvalues[0].index + 1, M_eigenvalues[0].order ) );
    
    // Insert the first root of the next order
    subeigen.insert( EigenMap2D::make_eigenmap2D( droots[M_eigenvalues[0].index][M_eigenvalues[0].order + 1], M_eigenvalues[0].index, M_eigenvalues[0].order + 1 ) );
    
	deleteBesselRoots( roots, droots );


    // Extract the next eigenvalues
    for( int j = 1; j != M_mtot; ++j )
    {

    	// extract the least eigenvalue among the following ones
        M_eigenvalues[j].lambda = ( subeigen.begin()->lambda ) * ( subeigen.begin()->lambda ) / ( M_Rho * M_Rho );
        M_eigenvalues[j].order  = subeigen.begin()->order;
        M_eigenvalues[j].index  = subeigen.begin()->index;
// std::cout << "Eval order " << M_eigenvalues[j].order << " index " << M_eigenvalues[j].index << " lambda " << M_eigenvalues[j].lambda << std::endl;
        // erase the extract element        
        subeigen.erase( subeigen.begin() );

        // compute the norm
        if( M_eigenvalues[j].order == 0 )
        {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / ( sqrt( 2 ) * sqrt( M_eigenvalues[j].lambda ) * M_Rho ) *
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho - M_eigenvalues[j].order * M_eigenvalues[j].order ) *
                          fabs( j0 ) ); 
        }  
        else if( M_eigenvalues[j].order == 1 )
        {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / ( sqrt( 2 ) * sqrt( M_eigenvalues[j].lambda ) * M_Rho ) *
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho - M_eigenvalues[j].order * M_eigenvalues[j].order ) *
                          fabs( j1 ) ); 
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

            M_normRJ[j] = 1. / ( 1. / ( sqrt( 2 ) * sqrt( M_eigenvalues[j].lambda ) * M_Rho ) *
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho - M_eigenvalues[j].order * M_eigenvalues[j].order ) *
                          fabs( jn[nm] ) ); 
			
            delete[] jn;
            delete[] yn;
            delete[] jnp;
            delete[] ynp;
        }

        // the first increment is because of the index of the matrix, the further one because I need one-past the current eigenvalue rightward and bottomward
       	bessel::besselzeros( M_eigenvalues[j].order + 2, M_eigenvalues[j].index + 2, roots, droots );

        // add next eigenvalue, same order
       	subeigen.insert( EigenMap2D::make_eigenmap2D( droots[(M_eigenvalues[j]).index + 1][M_eigenvalues[j].order],
                         M_eigenvalues[j].index + 1, M_eigenvalues[j].order  ) );
        
        // add eigenvalue same index, next order
        subeigen.insert( EigenMap2D::make_eigenmap2D( droots[M_eigenvalues[j].index][(M_eigenvalues[j]).order + 1],
                         M_eigenvalues[j].index, M_eigenvalues[j].order + 1 ) );
		
        deleteBesselRoots( roots, droots );

	}
	
    return;
}

} // end namespace bessel
