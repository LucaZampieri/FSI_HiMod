#ifndef __UTILITYUNSTEADYPRESSUREFUNCTIONS_HPP__
#define __UTILITYUNSTEADYPRESSUREFUNCTIONS_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/navier_stokes/function/bessel/bessel.hpp>
#include <cmath>
#include <complex>
#include <lifev/core/filter/GetPot.hpp>

GetPot dataFile( "data" );

// We define a set of data to assemble a problem with known solution
static LifeV::Real R = dataFile( "himod/rho", 1. );
static LifeV::Real dR = dataFile( "himod/drho", 0 );
static LifeV::Real Theta = dataFile( "himod/theta", 2. * M_PI );
static LifeV::Real Lx = dataFile ("mesh/lx", 1.);
static LifeV::Real density = dataFile( "fluid/physics/density", 1. );
static LifeV::Real nu = dataFile( "fluid/physics/nu", 0.1 );
static LifeV::Real mu = nu * density;
static LifeV::Real omega = dataFile( "fluid/problem/T", 2. * M_PI );
static LifeV::Real A = dataFile( "fluid/problem/A", 1. );
// The force term
LifeV::Real fx( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{   
	return 0;
}

LifeV::Real fxrho( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return fx( t, x, rho, theta, 0 ) * rho;
}

LifeV::Real fr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return 0;
}

LifeV::Real frrho( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return fr( t, x, rho, theta, 0 ) * rho;
}

LifeV::Real ftheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return 0;
}

LifeV::Real fthetarho( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return ftheta( t, x, rho, theta, 0 ) * rho;
}

// The exact solution
LifeV::Real x_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& c )
{
    
    std::complex<LifeV::Real>	one( 1., 0. );
    std::complex<LifeV::Real>	i( 0., 1. );
    std::complex<LifeV::Real>	p( std::real( r * sqrt( -i*omega / nu ) ),
                                   std::imag( r * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>	p0( std::real( R * sqrt( -i*omega / nu ) ),
                                    std::imag( R * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>  j0( 0, 0 ), j0p0( 0, 0 ),
    					j1( 0, 0 ), y0( 0, 0 ), y1( 0, 0 ), j0p( 0, 0 ), j1p( 0, 0 ), y0p( 0, 0 ), y1p( 0, 0 );
	
	bessel::cbessjy01( p, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
	bessel::cbessjy01( p0, j0p0, j1, y0, y1, j0p, j1p, y0p, y1p );
    
	return std::real( A * exp( i * omega * t ) / ( i * omega ) * ( one - j0 / j0p0 ) );
}

LifeV::Real r_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return 0;
}

LifeV::Real theta_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return 0;
}

LifeV::Real p_es( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    std::complex<LifeV::Real>	i( 0., 1. );
    
	return std::real( A * exp( i * omega * t ) * ( Lx - x ) );
}

// -----------------------------		Map		-----------------------------------

LifeV::Real Jr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 1. / R;
}

LifeV::Real Jtheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{   
    return 1. / R;
}

LifeV::Real Dr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return - dR / R * r;
}

LifeV::Real Dtheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real Drtheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real Dthetar( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

// IMPORTANT: This is actually the inverse of the Jacobian, i.e. the shape factor which shows up in the integrals!
LifeV::Real Jacobian( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{  
    return fabs( 1. / ( Jr( t, x, r, theta, 0 ) * Jtheta( t, x, r, theta, 0 ) - Drtheta( t, x, r, theta, 0 ) * Dthetar( t, x, r, theta, 0 ) ) );
}

LifeV::Real inverseRhat( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
	return R * r;
}

// ---------------------	Boundary Data		------------------------------

LifeV::Real xInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
// Dirichlet BC	
//    return x_ues( t, 0., r, theta, 0 );
// Neumann BC	
    return p_es( t, 0., r, theta, 0 );
}

LifeV::Real rInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
// Dirichlet
    return 0;
// Neumann
/*    std::complex<LifeV::Real>  one( 1., 0. );
    std::complex<LifeV::Real>  i( 0., 1. );
    std::complex<LifeV::Real>  p( std::real( r * sqrt( -i*omega / nu ) ),
                                   std::imag( r * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>  p0( std::real( R * sqrt( -i*omega / nu ) ),
                                    std::imag( R * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>  j0( 0, 0 ), j0p0( 0, 0 ),
                               j1( 0, 0 ), y0( 0, 0 ), y1( 0, 0 ), j0p( 0, 0 ), j1p( 0, 0 ), y0p( 0, 0 ), y1p( 0, 0 );
	
    bessel::cbessjy01( p, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
    bessel::cbessjy01( p0, j0p0, j1, y0, y1, j0p, j1p, y0p, y1p );
    
// Note: J0'(z) = -J1(z)
    return std::real( nu * A * exp( i * omega * t ) / ( i * omega ) * ( sqrt( -i*omega / nu ) * j1 / j0p0 ) );
*/
}

LifeV::Real thetaInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real xOut( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real rOut( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real thetaOut( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

#endif
