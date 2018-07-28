//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file CaseTest.hpp
     @brief This file contains a list of different test cases.

     @date 06/2013
     @author A. Bortolossi <andrea.bortolossi@gmail.com>
     @author M. Aletti <teo.aletti@gmail.com>
 */

#ifndef __CASETEST_HPP__
#define __CASETEST_HPP__ 1

/*
    CASE TEST:
*/

#include <lifev/himod/util/CaseTestStruct.hpp>
#include <lifev/himod/util/CaseTestStructCirc.hpp>

namespace LifeV
{
namespace CaseTest
{

/* CaseTestStruct(
                        const std::string& casename,
                        const Real& lx,
                        const Real& ly,
                        const Real& lz,
                        const Real& Mu,
                        const Real& beta1,

                        const Real& beta2,

                        const Real& beta3,
                        const Real& Sigma,
                        const Real& ChiY,
                        const Real& ChiZ,
                        const std::string& up,
                        const std::string& down,
                        const std::string& left,
                        const std::string& right,
                        const function_type& F,
                        const function_type& G,
                        const function_type& Ues
                        )
*/
Real g (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
Real f (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
Real ues (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStruct DDDD ("DDDD", 0.2, 0.1, 0.1, 1.0, 0.0, 0.0, 0.0, 0.0, 1, 1, "dir", "dir", "dir", "dir", f, g, ues);

/* Case test with full dirichlet boundary conditions and an advection field not zero.*/
Real g_adv (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
Real f_adv (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
Real ues_adv (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStruct DDDD_Adv ("DDDD_Adv", 0.2, 0.1, 0.1, 1.0, 5.0, 1.0, 1.0, 0.0, 1., 1., "dir", "dir", "dir", "dir", f_adv, g_adv, ues_adv);

Real g_adr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
Real f_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
Real ues_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStruct DDDD_ADR ("DDDD_ADR", 0.2, 0.1, 0.1, 1.0, 5.0, 1.0, 1.0, 3.0, 1.0, 1.0, "dir", "dir", "dir", "dir", f_adr, g_adr, ues_adr);

Real g_drdr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
Real f_drdr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
Real ues_drdr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStruct DRDR ("DRDR", 0.1, 0.1, 0.1, 1.0, 0.0, 0.0, 0.0, 0.0, 3.345, 0, "dir", "dir", "rob", "rob", f_drdr, g_drdr, ues_drdr);

Real g_rrrr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
Real f_rrrr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
Real ues_rrrr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStruct RRRR ("RRRR", 0.1, 0.1, 0.1, 1.0, 0.0, 0.0, 0.0, 0.0, 4.456, 4.456, "rob", "rob", "rob", "rob", f_rrrr, g_rrrr, ues_rrrr);


Real g_rr2d (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real f_rr2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real ues_rr2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
static const CaseTestStruct RR2D ("RR2D", /*lx*/ 1.0, /*ly*/1.0, /*lz*/1.0,
                                  /*mu*/1.0,
                                  /*beta1*/20.0, /*beta2*/ 2.0, /*beta3*/ 0.0,
                                  /*sigma*/2.0, /*chi*/ 1.0, 1.0,
                                  "basis", "fake", "rob", "rob", f_rr2d, g_rr2d, ues_rr2d);

Real g_2drr (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/);
Real f_2drr (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/);
Real ues_2drr (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/);
static const CaseTestStruct BDRR ("BDRR", /*lx*/ 1.0, /*ly*/1.0, /*lz*/1.0,
                                  /*mu*/1.0,
                                  /*beta1*/20.0, /*beta2*/ 0.0, /*beta3*/ 2.0,
                                  /*sigma*/2.0, /*chi*/ 1.0, 1.0,
                                  "rob", "rob", "fake", "basis", f_2drr, g_2drr, ues_2drr);

Real g_2ddd (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/);
Real f_2ddd (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/);
Real ues_2ddd (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/);
static const CaseTestStruct BDDD ("BDDD", /*lx*/ 1.0, /*ly*/1.0, /*lz*/1.0,
                                  /*mu*/1.0,
                                  /*beta1*/20.0, /*beta2*/ 0.0, /*beta3*/ 2.0,
                                  /*sigma*/2.0, /*chi*/ 1.0, 1.0,
                                  "dir", "dir", "fake", "basis", f_2ddd, g_2ddd, ues_2ddd);

Real g_dd2d (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real f_dd2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real ues_dd2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
static const CaseTestStruct DD2D ("DD2D", /*lx*/ 1.0, /*ly*/1.0, /*lz*/1.0,
                                  /*mu*/1.0,
                                  /*beta1*/20.0, /*beta2*/ 2.0, /*beta3*/ 0.0,
                                  /*sigma*/2.0, /*chi*/ 1.0, 1.0,
                                  "basis", "fake", "dir", "dir", f_dd2d, g_dd2d, ues_dd2d);

                                  
Real g_circ1 (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real f_circ1 (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real frho_circ1 (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real ues_circ1 (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
static const CaseTestStructCirc D ("D", /*lx*/ 1.0, /*rho*/0.1, /*theta*/2*M_PI,
                                  /*mu*/1.0,
                                  /*beta1*/3.0, /*beta2*/ 1.0, /*beta3*/ 1.0,
                                  /*sigma*/0.0, /*chi*/ 1.0,
                                  "dir", f_circ1, g_circ1, ues_circ1);
                                  
//------------------- SHOW CASE TESTS ----------------------------------

Real g_camini (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real f_camini (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real ues_camini (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/);

static const CaseTestStruct Camini ("Camini", /*lx*/ 0.2, /*ly*/0.1, /*lz*/0.1,
                                    /*mu*/0.05,
                                    /*beta1*/5.0, /*beta2*/ 1.0, /*beta3*/ 0.0,
                                    /*sigma*/0.3, /*chiY*/ 1.0, /*chiZ*/ 1.0,
                                    "dir", "dir", "dir", "dir", f_camini, g_camini, ues_camini);


Real ues_stent (const Real& /*t*/, const Real& /*x*/, const Real& /*rho*/, const Real& /*theta*/, const ID& /*i*/);
Real f_stent ( const Real& /* t */, const Real& x, const Real& rho, const Real&  theta  , const ID& /* i */ );
Real g_stent (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/);
static const CaseTestStructCirc stent ("stent", /*lx*/ 2.0, /*rho*/0.1, /*theta*/2*M_PI,
                                    /*mu*/0.05,
                                    /*beta1*/5.0, /*beta2*/ 1.0, /*beta3*/ 0.0,
                                    /*sigma*/0.3, /*chi*/ 1.0, 
                                    "dir", f_stent, g_stent, ues_stent);
                                  
}
} // end of LifeV namespace
#endif // __CASETEST__
