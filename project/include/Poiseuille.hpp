#ifndef __POISEUILLE_HPP__
#define __POISEUILLE_HPP__

#include <lifev/core/LifeV.hpp>
#include <cmath>
#include <lifev/core/filter/GetPot.hpp>

GetPot dataFile( "include/dataPoiseuille" );

// We define a set of data to assemble a problem with known solution
static LifeV::Real R = dataFile( "himod/rho", 1. );
static LifeV::Real a = dataFile( "fluid/physics/a", 1. );
static LifeV::Real b = dataFile( "fluid/physics/b", 0.1 );
static LifeV::Real dR = dataFile( "himod/drho", 0 );
static LifeV::Real Theta = dataFile( "himod/theta", 2. * M_PI );
static LifeV::Real Lx = dataFile ("mesh/lx", 1.);
static LifeV::Real density = dataFile( "fluid/physics/density", 1. );
static LifeV::Real nu = dataFile( "fluid/physics/nu", 0.1 );
static LifeV::Real mu = nu * density;
static LifeV::Real deltaP = dataFile( "fluid/physics/deltaP", 5. );
static LifeV::Real A = dataFile( "fluid/physics/A", 1. );
static LifeV::Real omega = 2 * M_PI;

static LifeV::Real Rin = dataFile( "fluid/physics/Rin", 0.2 );
static LifeV::Real Rout = dataFile( "fluid/physics/Rout", 0.2 );
static LifeV::Real z0 = dataFile( "fluid/physics/z0", 0.5 );
static LifeV::Real delta0 = dataFile( "fluid/physics/delta0", 0.15 );
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

static LifeV::Real fR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Cylinder
     return R;
    // Ellipse (for a = b the ellipse is a disk)
    //return a*b/std::sqrt( b*b* std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ) );

    // CASE (3) - non-axisymmetric: elliptical section, exponential x
    //LifeV::Real S( a*b/std::sqrt( b*b* std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ) ) );
    //return S*( 1 - 0.5*exp( -( x - Lx/2 )*( x - Lx/2 )/0.1 ) );

    // Cone
    // return ( Rout - Rin ) / Lx * x + Rin;

    // Stenosis
    // CASE (1) - exponential
    // return 1 - occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );

    // CASE (2) - axisymmetric
    /*LifeV::Real xShift(x-Lx/3);
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
static LifeV::Real fdR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    // Cylinder or Ellipse
     return 0;

    // CASE (3) - non-axisymmetric: elliptical section, exponential x
    //LifeV::Real S( a*b/std::sqrt( b*b* std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ) ) );
    //return S*( 0.5*exp( -( x - Lx/2 )*( x - Lx/2 )/0.1 )*2*( x - Lx/2 )/0.1 );

    // Cone
    // return ( Rout - Rin ) / Lx;

    // Stenosis
    // CASE (1) - exponential occlusion
    // return 2 * ( x - Lx / 2 ) * occlusion * exp( - ( x - Lx / 2 ) * ( x - Lx / 2 ) );

    // CASE (2) - axisymmetric
    /* LifeV::Real xShift(x-Lx/3);
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
    // REMARK. Dr is (-R'/R*r), but the r factor is included in modal space, so to factor the integrals.
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
    //LifeV::Real dRdt( a*b* ( b*b - a*a ) * std::cos( theta ) * std::sin( theta ) / pow( b*b*std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ), 3./2. ) );
    //return - dRdt/( fR( t, x, r, theta, 0 )*fR( t, x, r, theta, 0 ) );

    // Elliptical section, exponential x
    //LifeV::Real S( a*b/std::sqrt( b*b* std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ) ) );
    //LifeV::Real Sx( S*( 1 - 0.5*exp( -( x - Lx/2 )*( x - Lx/2 )/0.1 ) ) );
    //LifeV::Real dSdt( a*b* ( b*b - a*a ) * std::cos( theta ) * std::sin( theta ) / pow( b*b*std::cos( theta ) *std::cos( theta ) + a*a*std::sin( theta ) *std::sin( theta ), 3./2. )*( 1 - 0.5*exp( -( x - Lx/2 )*( x - Lx/2 )/0.1 ) ) );
    //return - dSdt/( fR( t, x, r, theta, 0 )*fR( t, x, r, theta, 0 ) );

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
        if( theta>2*M_PI-phix || theta < M_PI+phix )
        {
            return 0.;
        }
        else
        {
            // -cot * csc = -1 / (tan*sin)
            // Remark. Drtheta = -1/(R*R)*dR/dTheta
            return -( -1./( Rx*Rx*std::tan(theta-M_PI)*std::sin(theta-M_PI) )*(Rin-deltax) );
        }
    }
    */
}

// IMPORTANT: This is actually the inverse of the Jacobian, i.e. the shape factor which shows up in the integrals!
LifeV::Real Jacobian( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return fabs( 1. / ( Jr( t, x, r, theta, 0 ) * Jtheta( t, x, r, theta, 0 ) - Drtheta( t, x, r, theta, 0 ) * Dthetar( t, x, r, theta, 0 ) ) );
}

LifeV::Real JacobianWall( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return fR( t, x, r, theta, 0 );
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
LifeV::Real x_ues( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& i )
{
    return ( fR( t, x, r, theta, 0 ) * fR( t, x, r, theta, 0 ) - r * r ) / ( 4 * mu ) * deltaP / Lx;
    //return ( 1- r*r / ( fR( t, x, r, theta, 0 ) * fR( t, x, r, theta, 0 ) ) );
    //LifeV::Real u( 1/(2*mu)*A*a*a*b*b/(a*a+b*b) );
    //return u*( 1 - r*std::cos(theta)/a*r*std::cos(theta)/a - r*std::sin(theta)/b*r*std::sin(theta)/b );
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
    return A * ( Lx - x );
}

// ---------------------    Boundary Data        ------------------------------

LifeV::Real xInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return x_ues( t, 0., r, theta, 0 );
}

LifeV::Real rInf( const LifeV::Real& t, const LifeV::Real& /*x*/, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
    return 0;
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
