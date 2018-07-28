#pragma GCC diagnostic ignored "-Wsign-compare"

#include <lifev/himod/basis/EducatedBasisR.hpp>
#include <lifev/core/algorithm/NonLinearBrent.hpp>

namespace LifeV
{

Real EducatedBasisR::
Bisection( functor_ptrType f, Real a, Real b, const Real tol, const int maxIt )
{
	Real fa, fb, fx( 100 ), x;
	
	fa=( *f )( a );
	if (fa == 0)
	{
		return a;
	}
	
	fb=( *f )( b );
	if ( fb == 0 )
	{
		 return b;
	}

	assert( b > a );
	assert( fa * fb < 0 );
	UInt i=0;
	for ( i = 0; i != maxIt && fabs( fx ) > tol; ++i )
	{
		x=( a + b ) / 2.;
		fx=( *f )( x );
		if ( fx == 0 )
			return x;
		if ( fa * fx < 0 )
		{
			b = x;
		}
		else
		{
			a = x;
		}
	}
	return x;
}

void EducatedBasisR::
computeEigenvalues()
{	

	MBMatrix_type	 roots, edges, dedges;
    std::set< EigenMap2D, Comparison2D > subeigen;
    Real  j0, j1, y0, y1, j0p, j1p, y0p, y1p;
    Real omega, a, b;
    int  nm;
    Real eps1 = 1e-12;
    
    M_eigenvalues.resize( M_mtot );
    M_normRJ.resize( M_mtot );
    roots.resize( M_mtot );
    for( UInt i = 0; i != M_mtot; ++i )
    	roots[i].resize( M_mtot );
 
    bessel::besselzeros( M_mtot, M_mtot, edges, dedges );    
     
    // The first root of the functor is between 0 and the first root of the akin Bessel function
    for( int i = 0; i != M_mtot; ++i ) //order
    {
        M_ptrFunctor = static_cast<functor_ptrType> ( new  EducatedBasisFunctorR( M_mu, M_chi, M_Rho, i ) );
	    
       	a = eps1;
       	b = edges[0][i];
        omega = Bisection( M_ptrFunctor, a, b, eps1 / 10., 10000 );
//      std::cout<<"root order "<<i<<": "<<omega<<std::endl;
        roots[0][i] = omega; // the root of the functor is the eigenvalue on the reference domain
        subeigen.insert( EigenMap2D::make_eigenmap2D( omega, 0, i ) );
//      std::cout<<"root order "<< i <<" index "<< 0 <<": "<<std::setprecision(16)<<omega<<std::endl;
    }
   
	
    // Find the following roots for all the orders
    for( int i = 0; i != M_mtot; ++i ) //order
    {
        M_ptrFunctor = static_cast<functor_ptrType> ( new  EducatedBasisFunctorR( M_mu, M_chi, M_Rho, i ) );
	    
        for( int j = 1; j != M_mtot; ++j ) //index
        {
        	a = edges[j-1][i];
        	b = edges[j][i];
                omega = Bisection( M_ptrFunctor, a, b, eps1 / 10., 10000 );
//		std::cout<<"root order "<<i<<" index "<<j<<": "<<std::setprecision(16)<<omega<<std::endl;
	        roots[j][i] = omega; // the root of the functor is the eigenvalue on the reference domain
                subeigen.insert( EigenMap2D::make_eigenmap2D( omega, j, i ) );
        }
    }
 
    for( int j = 0; j != M_mtot; ++j )
    {
        M_eigenvalues[j].lambda=( *subeigen.begin() ).lambda * ( *subeigen.begin() ).lambda / ( M_Rho * M_Rho );
        M_eigenvalues[j].order=( *subeigen.begin() ).order;
        M_eigenvalues[j].index=( *subeigen.begin() ).index;
//        std::cout<<"eigenvalue "<<j<<": "<< std::setprecision(16)<< sqrt(M_eigenvalues[j].lambda)<<" order "<<M_eigenvalues[j].order<<std::endl;
        subeigen.erase( subeigen.begin() );
        
        // compute the norm
        if( M_eigenvalues[j].order == 0 )
        {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / ( sqrt( 2 ) * sqrt( M_eigenvalues[j].lambda ) * M_Rho ) *
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho + M_chi / M_mu * M_chi / M_mu * M_Rho * M_Rho -
                          M_eigenvalues[j].order * M_eigenvalues[j].order ) *
                          fabs( j0 ) );
        }  
        else if( M_eigenvalues[j].order == 1 )
        {
            bessel::bessjy01b( sqrt( M_eigenvalues[j].lambda ) * M_Rho, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
            M_normRJ[j] = 1. / ( 1. / ( sqrt( 2 ) * sqrt( M_eigenvalues[j].lambda ) * M_Rho ) *
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho + M_chi / M_mu * M_chi / M_mu * M_Rho * M_Rho -
                          M_eigenvalues[j].order * M_eigenvalues[j].order ) *
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
                          sqrt( M_eigenvalues[j].lambda * M_Rho * M_Rho + M_chi / M_mu * M_chi / M_mu * M_Rho * M_Rho -
                          M_eigenvalues[j].order * M_eigenvalues[j].order ) *
                          fabs( jn[nm] ) ); 
            delete[] jn;
            delete[] yn;
            delete[] jnp;
            delete[] ynp;
        }

    }
    return;
}

} // namespace lifev
