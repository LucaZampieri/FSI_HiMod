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
         D_t0( dataFile( "time/t0", 0. ) ),
         D_dt( dataFile( "time/dt", 0.1 ) ),
         D_alpha( 1./D_dt ),
         D_T( dataFile( "time/T", 1. ) ),

         D_theta( 2*PI ),
         D_L( dataFile( "fluid/structure/L", 5. ) ),
         D_rho_s( dataFile( "fluid/structure/rho_s", 1. ) ),
         D_h_s( dataFile( "fluid/structure/h_s", 0.1 ) ),
         D_e( dataFile( "fluid/structure/e", 4.e5 ) ),
         D_csi( dataFile( "fluid/structure/csi", 1. ) ),
         D_rho_f( dataFile( "fluid/physics/rho_f", 1. ) ),
         D_nu( dataFile( "fluid/physics/nu", 1. ) ),
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
         //D_e( D_E * D_h_s / ((1 - D_csi*D_csi) * D_R*D_R)),
>>>>>>> f65a88a19c8aa250af643994d791a1f56cab1715
>>>>>>> 24c6253a84befbfd29040e8c40b4ff0024900469

         D_ux0_str( dataFile( "functions/ux0", "0" ) ),
         D_ur0_str( dataFile( "functions/ur0", "0" ) ),
         D_utheta0_str( dataFile( "functions/utheta0", "0" ) ),
         D_p1_str( dataFile( "functions/p1", "0" ) ),
         D_p2_str( dataFile( "functions/p2", "0" ) ),
         D_fx_str( dataFile( "functions/fx", "0" ) ),
         D_fr_str( dataFile( "functions/fr", "0" ) ),
         D_ftheta_str( dataFile( "functions/ftheta", "0" ) ),
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
         D_Radius(  muparser_function( D_Radius_str ) ),
         D_dRadius( muparser_function( D_dRadius_str ) )
{
  std::cout <<"Constructed FSIData" <<std::endl;
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

  std::cout << "radius and dradius blahblahblah"<<std::endl;

  std::cout <<"theta: " <<D_theta <<std::endl;
  std::cout <<"L: " <<D_L <<std::endl;
  std::cout <<"rho_s: " <<D_rho_s <<std::endl;
  std::cout <<"h_s: " <<D_h_s <<std::endl;
  std::cout <<"E: " <<D_E <<std::endl;
  std::cout <<"csi: " <<D_csi <<std::endl;

  std::cout <<"rho_f: " <<D_rho_f <<std::endl;
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
