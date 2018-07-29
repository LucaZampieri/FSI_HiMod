#include <include/FSIData.hpp>

#include <iostream>
#include <iomanip> // for nicer output of print all

namespace LifeV
{

FSIData::FSIData( GetPot dataFile ) :
         // HiMod
         D_mx( dataFile( "himod/mx", 10 ) ),
         D_mr( dataFile( "himod/mr", 10 ) ),
         D_mtheta( dataFile( "himod/mtheta", 10 ) ),
         D_mp( dataFile( "himod/mp", 8 ) ),
         D_Nelements( dataFile( "mesh/num_elements", 10 ) ),

         // time
         D_t0( dataFile( "time/t0", 0. ) ),
         D_dt( dataFile( "time/dt", 0.1 ) ),
         D_alpha( 1./D_dt ),
         D_T( dataFile( "time/T", 1. ) ),

         D_L(        dataFile( "fluid/structure/L", 5. ) ),
         D_rho_s(    dataFile( "fluid/structure/rho_s", 1. ) ),
         D_h_s(      dataFile( "fluid/structure/h_s", 0.1 ) ),
         D_e(        dataFile( "fluid/structure/e", 4.e5 ) ),
         D_occlusion(dataFile( "fluid/structure/occlusion", 0.45 ) ),
         D_nu(       dataFile( "fluid/physics/nu", 1. ) ),

         // Structure
         D_Rin(       dataFile( "fluid/structure/Rin" , 1. ) ),
         D_Rout(      dataFile( "fluid/structure/Rout", 1. ) ),
         D_Z0(        dataFile( "fluid/structure/Z0"  , 1. ) ),
         D_delta0(    dataFile( "fluid/structure/delta0", 0.45 ) ),

         D_ux0_str(     dataFile( "functions/ux0",     "0" ) ),
         D_ur0_str(     dataFile( "functions/ur0",     "0" ) ),
         D_utheta0_str( dataFile( "functions/utheta0", "0" ) ),
         D_p1_str(      dataFile( "functions/p1",      "0" ) ),
         D_p2_str(      dataFile( "functions/p2",      "0" ) ),
         D_fx_str(      dataFile( "functions/fx",      "0" ) ),
         D_fr_str(      dataFile( "functions/fr",      "0" ) ),
         D_ftheta_str(  dataFile( "functions/ftheta",  "0" ) ),
         D_Radius_str(  dataFile( "functions/Radius" , "0" ) ),
         D_dRadius_str( dataFile( "functions/dRadius", "0" ) ),


         D_ux0( muparser_function( D_ux0_str ) ),
         D_ur0( muparser_function( D_ur0_str ) ),
         D_utheta0( muparser_function( D_utheta0_str ) ),
         D_p1( muparser_timeFunction( D_p1_str ) ),
         D_p2( muparser_timeFunction( D_p2_str ) ),
         D_fx( muparser_function( D_fx_str ) ),
         D_fr( muparser_function( D_fr_str ) ),
         D_ftheta(  muparser_function( D_ftheta_str ) ),
         D_theta( 2*PI ),

         D_Jr( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 1. / D_Radius(t,x,r,theta,0); } ),
         D_Jtheta( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 1. / D_Radius(t,x,r,theta,0); } ),
         D_Dr( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return - D_dRadius(t,x,r,theta,0) / D_Radius(t,x,r,theta,0); } ),
         D_Dtheta( [] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 0; } ),
         D_Drtheta( [] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 0; } ),
         D_Dthetar( [] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 0; } ),
         D_Jacobian( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return fabs( 1. / ( D_Jr(t,x,r,theta,0) * D_Jtheta(t,x,r,theta,0) - D_Drtheta(t,x,r,theta,0) * D_Dthetar(t,x,r,theta,0) ) ); } ),
         D_JacobianWall( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return D_Radius(t,x,r,theta,0); } ),
         D_inverseRhat( [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return D_Radius(t,x,r,theta,0)*r; } )
{
  std::cout <<"Constructed FSIData" <<std::endl;

    int case_radius = dataFile( "fluid/structure/case_radius", 0 );
    switch (case_radius)
    {
        case 0:
            D_Radius = muparser_function( D_Radius_str );
            D_dRadius = muparser_function( D_dRadius_str );
        break;
        case 1: // stenosis : Exponential Occlusion
            D_Radius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 1 - /*occlusion*/0.45 * exp( - ( x - D_L / 2 ) * ( x - D_L / 2 ) ); };
            D_dRadius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 2 * ( x - D_L / 2 ) * /*occlusion*/0.45 * exp( - ( x - D_L / 2 ) * ( x - D_L / 2 ) ); };
        break;
        case 2:  // Aneurysm
            D_Radius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return 0.5 + /*occlusion*/0.45 * exp( - ( x - D_L / 2 ) * ( x - D_L / 2 ) ); };
            D_dRadius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return - 2 * ( x - D_L / 2 ) * /*occlusion*/0.45 * exp( - ( x - D_L / 2 ) * ( x - D_L / 2 ) ); };
        break;
        case 3: // Cone
            D_Radius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return ( D_Rout - D_Rin ) / D_L * x + D_Rin; };
            D_dRadius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ ) { return ( D_Rout - D_Rin ) / D_L; };
        break;
        case 4: // axisymmetric
            D_Radius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ )
            {
              Real xShift(x-D_L/3);
              if( xShift>-D_Z0 && xShift<D_Z0 )
              {
                  return D_Rin * ( 1 - D_delta0 / (2*D_Rin) * ( 1 + std::cos(M_PI*xShift/D_Z0) ) );
              }
              else
              {
                  return D_Rin;
              }
            };
            D_dRadius = [this] ( const Real& t, const Real& x, const Real& r, const Real& theta, const ID& /*i*/ )
            {
              Real xShift(x-D_L/3);
             if( xShift>-D_Z0 && xShift<D_Z0 )
             {
                 return D_delta0/2 * M_PI/D_Z0 * std::sin(M_PI*xShift/D_Z0);
             }
             else
             {
                 return 0.;
             }
            };
        break;
    }
}

void FSIData::printAll() const
{
  std::cout <<"The data of the fluid structure interaction problem are:" <<std::endl;
  std::cout << "Number of modes: \n";
  std::cout << std::setw(5)  << std::left     <<"mx: "     <<D_mx;
  std::cout << std::setw(10) << std::internal <<"mr: "     <<D_mr;
  std::cout << std::setw(10) << std::internal <<"mtheta: " <<D_mtheta;
  std::cout << std::setw(5)  << std::right    <<"mp: "     <<D_mp;

  std::cout <<"Nelements: " <<D_Nelements <<std::endl;

  std::cout <<"t0: " <<D_t0 <<std::endl;
  std::cout <<"dt: " <<D_dt <<std::endl;
  std::cout <<"T: "  <<D_T  <<std::endl;

  std::cout << "Radius "  << D_Radius_str << std::endl;
  std::cout << "dRadius " << D_dRadius_str << std::endl;


  std::cout <<"L: " <<D_L <<std::endl;
  std::cout <<"rho_s: " <<D_rho_s <<std::endl;
  std::cout <<"h_s: " <<D_h_s <<std::endl;
  std::cout <<"e: " <<D_e <<std::endl;

  std::cout <<"nu: " <<D_nu <<std::endl;

  std::cout <<"ux0 = " <<D_ux0_str <<std::endl;
  std::cout <<"ur0 = " <<D_ur0_str <<std::endl;
  std::cout <<"utheta0 = " <<D_utheta0_str <<std::endl;

  std::cout <<"p1 = " <<D_p1_str <<std::endl;
  std::cout <<"p2 = " <<D_p2_str <<std::endl;

  std::cout <<"fx = " <<D_fx_str <<std::endl;
  std::cout <<"fr = " <<D_fr_str <<std::endl;
  std::cout <<"ftheta = " <<D_ftheta_str <<std::endl;
}

}
