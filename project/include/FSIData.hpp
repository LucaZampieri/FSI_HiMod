#ifndef FSI_DATA_HPP
#define FSI_DATA_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <cmath>
#include <functional>
#include "include/muparser_function.hpp"
#include "include/muparser_timeFunction.hpp"
#include <string>

namespace LifeV
{

/* Class which collect all the data to solve a fluid structure interaction problem */

class FSIData
{
  public:

    static constexpr Real PI = std::atan(1.)*4; // pi greco
    typedef std::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& )> function_Type;
    typedef std::function<Real ( const Real& )> oneDFunction_Type;

    // All the functions that return the data
    const UInt& mx() const
    {return D_mx;}

    const UInt& mr() const
    {return D_mr;}

    const UInt& mtheta() const{
      return D_mtheta;
    }

    const UInt& mp() const{
      return D_mp;
    }

    const UInt& Nelements() const
    {
      return D_Nelements;
    }

    const Real& t0() const
    {
      return D_t0;
    }
    const Real& dt() const
    {
      return D_dt;
    }
    const Real& T() const{
      return D_T;
    }

    const Real& theta() const
    {
      return D_theta;
    }

    const Real& L() const
    {
      return D_L;
    }

    const Real& rho_s() const
    {
      return D_rho_s;
    }
    const Real& h_s() const
    {
      return D_h_s;
    }

    const Real& nu() const
    {
      return D_nu;
    }

    const function_Type& ux0() const
    {
      return D_ux0;
    }
    const function_Type& ur0() const
    {
      return D_ur0;
    }
    const function_Type& utheta0() const
    {
      return D_utheta0;
    }

    const oneDFunction_Type& p1() const
    {
      return D_p1;
    }
    const oneDFunction_Type& p2() const
    {
      return D_p2;
    }

    const function_Type& fx() const
    {
      return D_fx;
    }

    const function_Type& fr() const
    {
      return D_fr;
    }

    const function_Type& ftheta() const
    {
      return D_ftheta;
    }

    const function_Type& Radius() const
    {
      return D_Radius;
    }

    const function_Type& dRadius() const
    {
      return D_dRadius;
    }

    const Real& alpha() const
    {
      return D_alpha;
    }

    const Real& e() const
    {
      return D_e;
    }



    // Constructor
    FSIData( GetPot dataFile );

    // Function which print all the data
    void printAll() const;

  private:

    UInt D_mx;                // number of modal functions for the axial velocity
    UInt D_mr;                // number of modal functions for the radial velocity
    UInt D_mtheta;            // number of modal functions for the angular velocity
    UInt D_mp;                // number of modal functions for the pressure

    UInt D_Nelements;         // number of finite elements along the axial direction

    Real D_t0;                // initial time
    Real D_dt;                // time step
    Real D_alpha;              // 1/dt
    Real D_T;                 // final time

    Real D_theta;             // 2*pigreco
    Real D_L;                 // length of the cylinder
    Real D_rho_s;             // density of the structure
    Real D_h_s;               // thickness of the structure
    Real D_e;                 // Coefficient structure

    Real D_nu;                // dynamical viscosity of the fluid

    std::string D_ux0_str;      // string of the initial axial velocity
    std::string D_ur0_str;      // string of the initial radial velocity
    std::string D_utheta0_str;  // string of the initial angular velocity

    std::string D_p1_str;   // string of the pressure at the inlet
    std::string D_p2_str;   // string of the pressure at the outlet

    std::string D_fx_str;       // string of the volumetric axial force
    std::string D_fr_str;       // string of the volumetric radial force
    std::string D_ftheta_str;   // string of the volumetric angular

    std::string D_Radius_str;   // string of the radius
    std::string D_dRadius_str;  // string of the derivative of the radius

    function_Type D_ux0;      // initial axial velocity
    function_Type D_ur0;      // initial radial velocity
    function_Type D_utheta0;  // initial angular velocity

    oneDFunction_Type D_p1;        // pressure at the inlet
    oneDFunction_Type D_p2;        // pressure at the outlet

    function_Type D_Radius;   // radius (fct of x)
    function_Type D_dRadius;  // Derivative of the radius wrt x

    function_Type D_fx;       // volumetric axial force
    function_Type D_fr;       // volumetric radial force
    function_Type D_ftheta;   // volumetric angular force


};

}

#endif
