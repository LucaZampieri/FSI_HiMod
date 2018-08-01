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

    const Real& occlusion() const
    {
      return D_occlusion;
    }

    const Real& Rin() const
    {
      return D_Rin;
    }

    const Real& Rout() const
    {
      return D_Rout;
    }

    const Real& Z0() const
    {
      return D_Z0;
    }

    const Real& delta0() const
    {
      return D_delta0;
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

    function_Type Jr() const
    {
      return D_Jr;
    }

    function_Type Jtheta() const
    {
      return D_Jtheta;
    }

    function_Type Dr() const
    {
      return D_Dr;
    }

    function_Type Dtheta() const
    {
      return D_Dtheta;
    }

    function_Type Drtheta() const
    {
      return D_Drtheta;
    }

    function_Type Dthetar() const
    {
      return D_Dthetar;
    }

    function_Type Jacobian() const
    {
      return D_Jacobian;
    }

    function_Type JacobianWall() const
    {
      return D_JacobianWall;
    }

    function_Type inverseRhat() const
    {
      return D_inverseRhat;
    }

    std::string polyTypeVelocity() const
    {
      return D_polyTypeVelocity;
    }

    std::string polyTypePressure() const
    {
      return D_polyTypePressure;
    }


    // Constructor
    FSIData( GetPot dataFile );

    // Function which print all the data
    void printAll() const;

  private:

    bool D_dataFileHasBeenRead; // Boolean to check if the datafile has been read

    UInt D_mx;                // number of modal functions for the axial velocity
    UInt D_mr;                // number of modal functions for the radial velocity
    UInt D_mtheta;            // number of modal functions for the angular velocity
    UInt D_mp;                // number of modal functions for the pressure

    UInt D_Nelements;         // number of finite elements along the axial direction

    UInt D_case_radius;        // Case for the geometry

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
    Real D_occlusion;         // Occlusion of the cylinder (used for cases 1&2)
    Real D_Rin;               // input Radius for cases 3&4
    Real D_Rout;              // output radius for cases 3&4
    Real D_Z0;                // Mioddle section for case 4
    Real D_delta0;                // Mioddle section for case 4

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

    std::string D_polyTypeVelocity;
    std::string D_polyTypePressure;

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

    // Functions for the reference map;
    function_Type D_Jr;
    function_Type D_Jtheta;
    function_Type D_Dr;
    function_Type D_Dtheta;
    function_Type D_Drtheta;
    function_Type D_Dthetar;
    function_Type D_Jacobian;
    function_Type D_JacobianWall;
    function_Type D_inverseRhat;

};

}

#endif
