#ifndef __UTILITYFUNCTIONSPIPE_HPP__
#define __UTILITYFUNCTIONSPIPE_HPP__

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

static LifeV::Real Rin = dataFile( "fluid/physics/Rin", 1. );
static LifeV::Real Rout = dataFile( "fluid/physics/Rout", 0.5 );
static LifeV::Real z0 = dataFile( "fluid/physics/z0", 1. );
static LifeV::Real delta0 = dataFile( "fluid/physics/delta0", 0.2 );
static LifeV::Real Z0( z0*std::sin( std::acos( (z0-delta0)/z0 ) ) );
//static LifeV::Real Z0 = dataFile( "physics/Z0", 4*Rin );

static LifeV::Real occlusion = dataFile( "fluid/physics/occlusion", 0.5 );
static LifeV::Real delta( occlusion );

// -----------------------------        Map        -----------------------------------
static LifeV::Real psi( const LifeV::Real& x )
{
    LifeV::Real xShift(x-Lx/2);

    if( abs(xShift)>Z0 )
        return std::asin( Z0/z0 );
    else
        return std::asin( abs(xShift)/z0 );
    
}

static LifeV::Real fDelta( const LifeV::Real& x )
{
    LifeV::Real xShift(x-Lx/2);

    if( abs(xShift)>Z0 )
        return 0;
    else
        return z0 * std::cos( psi(x) ) - (z0-delta0);
    
}

static LifeV::Real phi( const LifeV::Real& x )
{
    LifeV::Real xShift(x-Lx/2);

    if( abs(xShift)>Z0 )
        return M_PI/2.;
    else
        return std::asin(1-fDelta(x)/Rin);
    
}

static LifeV::Real alpha( const LifeV::Real& x, const LifeV::Real& theta )
{
    LifeV::Real phix( phi(x) );
    if( theta > phix && theta < 3/2*M_PI )
        return theta-M_PI;
    else if( theta >= 3/2*M_PI && theta < 2*M_PI-phix )
        return 2*M_PI-theta;
    else
        return 0;
}

static LifeV::Real dalpha( const LifeV::Real& x, const LifeV::Real& theta )
{
    LifeV::Real phix( phi(x) );
    if( theta > phix && theta < 3/2*M_PI )
        return 1;
    else if( theta >= 3/2*M_PI && theta < 2*M_PI-phix )
        return -1;
    else
        return 0;
}

LifeV::Real fR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Cylinder
     return R;
    
    // Cone
    // return ( Rout - Rin ) / Lx * x + Rin;
    
    // Stenosis
    // CASE (1) - exponential
    // return 1 - occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );
    
    // CASE (2) - axisymmetric
    /* LifeV::Real xShift(x-Lx/2);
    if( xShift>-Z0 && xShift<Z0 )
    {
        return Rin * ( 1 - delta / (2*Rin) * ( 1 + std::cos(M_PI*xShift/Z0) ) );
    }
    else
    {
        return Rin;
    }
    */

    // CASE (3) - non-axisymmetric
    /*LifeV::Real xShift( x-Lx/2 );
    LifeV::Real deltax( fDelta(x) );
    LifeV::Real phix( phi(x) );
    
    if( abs(xShift)>Z0 )
    {
        return Rin;
    }
    else
    {
        if( theta>2*M_PI-phix || theta < M_PI+phix )
        {
            return Rin;
        }
        else
        {
            return ( Rin - deltax )/std::sin( theta-M_PI );
        }
    }
*/
    // Aneurysm
    // return 0.5 + occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );
}
LifeV::Real fdR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Cylinder
     return 0;
    
    // Cone
    // return ( Rout - Rin ) / Lx;
    
    // Stenosis
    // CASE (1) - exponential occlusion
    // return 2 * ( x - Lx / 2 ) * occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );

    // CASE (2) - axisymmetric
    /* LifeV::Real xShift(x-Lx/2);
    if( xShift>-Z0 && xShift<Z0 )
    {
        return delta/2 * M_PI/Z0 * std::sin(M_PI*xShift/Z0);
    }
    else
    {
        return 0.;
    }
    */

    // CASE (3) - non-axisymmetric
    /*LifeV::Real xShift(x-Lx/2);
    LifeV::Real phix( phi(x) );
    LifeV::Real ddeltax( -xShift / ( z0*sqrt( 1 - xShift*xShift/(z0*z0) ) ) );
    
    if( abs(xShift)>Z0 )
    {
        return 0.;
    }
    else
    {
        if( theta>2*M_PI-phix || theta < M_PI+phix )
        {
            return 0.;
        }
        else
        {
            return -ddeltax/std::sin( theta-M_PI );
        }
    }
*/
    // Aneurysm
    // return - 2 * ( x - Lx / 2 ) * occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );
}

LifeV::Real Jr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 1. / fR( t, x, r, theta, 0 );
}

LifeV::Real Jtheta( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{   
    return 1. / fR( t, x, r, theta, 0 );
}

LifeV::Real Dr( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return - fdR( t, x, r, theta, 0 ) / fR( t, x, r, theta, 0 );
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

    // CASE (3) - non-axisymmetric stenosis
    /*LifeV::Real xShift( x-Lx/2 );
    LifeV::Real phix( phi(x) );
    LifeV::Real deltax( fDelta(x) );
    LifeV::Real Rx( fR(t,x,r,theta,0) );

    if( abs(xShift)>Z0 )
    {
        return 0.;
    }
    else
    {
        if( theta>=2*M_PI-phix || theta <= M_PI+phix )
        {
            return 0.;
        }
        else
        {
            // -cot * csc = -1 / (tan*sin)
            // Remark. Drtheta = -1/(R*R)*dR/dTheta
            return -( -1./( Rx*Rx*std::tan(alpha(x,theta))*std::sin(alpha(x,theta)) )*(Rin-deltax)*dalpha(x,theta) );
        }
    }
*/
}

// IMPORTANT: This is actually the inverse of the Jacobian, i.e. the shape factor which shows up in the integrals!
LifeV::Real Jacobian( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{  
    return fabs( 1. / ( Jr( t, x, r, theta, 0 ) * Jtheta( t, x, r, theta, 0 ) - Drtheta( t, x, r, theta, 0 ) * Dthetar( t, x, r, theta, 0 ) ) );
}

LifeV::Real inverseRhat( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return fR( t, x, r, theta, 0 ) * r;
}

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
LifeV::Real x_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& c )
{
    LifeV::Real Rx( fR( t, x, r, theta, 0 ) );
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
    
    /*LifeV::Real Rx( fR( t, x, r, theta, 0 ) );
    LifeV::Real Qt( t*t );
    //LifeV::Real Qt( -23.1754*t*t*t*t+39.0591*t*t*t-18.5291*t*t+2.6454*t+0.45 );
    LifeV::Real A( Qt*2 );
    return A * ( 1 - r*r/(Rx*Rx) );
 */   
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
    std::complex<LifeV::Real>    i( 0., 1. );
    
    return std::real( A * exp( i * omega * t ) * ( Lx - x ) );
}

// ---------------------    Boundary Data        ------------------------------

LifeV::Real xInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
// Dirichlet BC    
    return x_ues( t, 0., r, theta, 0 );
    //return (1-r*r)*2;
// Neumann BC    
//    return p_es( t, 0., r, theta, 0 );
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
