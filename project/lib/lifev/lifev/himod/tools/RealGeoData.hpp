#ifndef __REALGEODATA_HPP__
#define __REALGEODATA_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/navier_stokes/function/bessel/bessel.hpp>
#include <cmath>
#include <complex>
#include <lifev/core/filter/GetPot.hpp>

GetPot dataFile( "data" );

// We define a set of data to assemble a problem with known solution
static LifeV::Real Theta = dataFile( "himod/theta", 2. * M_PI );
static LifeV::Real Lx = dataFile ("mesh/lx", 1.);
static LifeV::Real density = dataFile( "fluid/physics/density", 1. );
static LifeV::Real nu = dataFile( "fluid/physics/nu", 0.1 );
static LifeV::Real mu = nu * density;
static LifeV::Real deltaP = dataFile( "fluid/physics/deltaP", 5. );
static LifeV::Real omega = dataFile( "fluid/physics/T", 2. * M_PI );
static LifeV::Real A = dataFile( "fluid/physics/A", 1. );

static LifeV::Real Rin = dataFile( "fluid/physics/Rin", 0.2 );

// The force term
LifeV::Real fx( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{   
    return 0;
}

LifeV::Real fr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

LifeV::Real ftheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& rho, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
}

// The exact solution
LifeV::Real x_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& id )
{
    // Poiseuille
    // return ( Rin * Rin - r * r ) / ( 4 * mu ) * deltaP / Lx;
    
    // Womersley
    /*LifeV::Real Rx( Rin );
    std::complex<LifeV::Real>    one( 1., 0. );
    std::complex<LifeV::Real>    i( 0., 1. );
    std::complex<LifeV::Real>    p( std::real( r * sqrt( -i*omega / nu ) ),
                                   std::imag( r * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>    p0( std::real( Rx * sqrt( -i*omega / nu ) ),
                                    std::imag( Rx * sqrt( -i*omega / nu ) ) );
    std::complex<LifeV::Real>  j0( 0, 0 ), j0p0( 0, 0 ),
                        j1( 0, 0 ), y0( 0, 0 ), y1( 0, 0 ), j0p( 0, 0 ), j1p( 0, 0 ), y0p( 0, 0 ), y1p( 0, 0 );
    
    bessel::cbessjy01( p, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
    bessel::cbessjy01( p0, j0p0, j1, y0, y1, j0p, j1p, y0p, y1p );
    
    return std::real( A * exp( i * omega * t ) / ( i * omega ) * ( one - j0 / j0p0 ) );
    */
    LifeV::Real Rx( Rin );
    //LifeV::Real Qt( -23.1754*t*t*t*t+39.0591*t*t*t-18.5291*t*t+2.6454*t+0.45 );
    LifeV::Real Qt( -15.8504*t*t*t*t+28.8767*t*t*t-15.8001*t*t+2.7738*t+0.45 );
    LifeV::Real A( Qt*2 );
    return A * ( 1 - r*r/(Rx*Rx) );

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
    //return deltaP / Lx * ( Lx - x );

    std::complex<LifeV::Real> i( 0., 1. );
    return std::real( A * exp(i*omega*t)*( Lx - x ) );
}

// ---------------------    Boundary Data        ------------------------------

LifeV::Real xInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Dirichlet
     return x_ues( t, 0., r, theta, 0 );
    
    // Neumann
    //return p_es( t, 0., r, theta, 0 );
}

LifeV::Real rInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Dirichlet
    return 0;
    //
    // Neumann
    /*std::complex<LifeV::Real>  one( 1., 0. );
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
