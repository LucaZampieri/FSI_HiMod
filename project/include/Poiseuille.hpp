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
static LifeV::Real mu2 = nu * density; // Luca
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

static LifeV::Real fR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
     return R;
}
static LifeV::Real fdR( const LifeV::Real& t, const LifeV::Real& x, const LifeV::Real& r, const LifeV::Real& theta, const LifeV::ID& /*i*/ )
{
     return 0;
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




#endif
