//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file CaseTestStruct.hpp
    @brief A very simple structure to contain case test data

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi@gmail.com>

 */
#ifndef __CASETESTSTRUCT_HPP__
#define __CASETESTSTRUCT_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{

class CaseTestStruct
{
public:


    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID& ) >      function_type;


    CaseTestStruct (
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
    ) :
        CaseName (casename), Lx (lx), Ly (ly), Lz (lz), mu (Mu), sigma (Sigma), chiy (ChiY), chiz (ChiZ),
        UP (up), DOWN (down), LEFT (left), RIGHT (right), f (F), g (G), ues (Ues)
    {
        beta[0] = beta1;
        beta[1] = beta2;
        beta[2] = beta3;
    };

    CaseTestStruct() {};

    CaseTestStruct (const CaseTestStruct& C)
    {
        CaseName = C.CaseName;
        Lx = C.Lx;
        Ly = C.Ly;
        Lz = C.Lz;
        mu = C.mu;
        beta = C.beta;
        sigma = C.sigma;
        chiy = C.chiy;
        chiz = C.chiz;
        UP = C.UP;
        DOWN = C.DOWN;
        LEFT = C.LEFT;
        RIGHT = C.RIGHT;
        f = C.f;
        g = C.g;
        ues = C.ues;
    };

    std::string CaseName;
    Real Lx;
    Real Ly;
    Real Lz;
    Real mu;
    VectorSmall<3> beta;
    Real sigma;
    Real chiy;
    Real chiz;
    std::string UP;
    std::string DOWN;
    std::string LEFT;
    std::string RIGHT;
    function_type f;
    function_type g;
    function_type ues;
};

} //End lifev namespace

#endif //End __CASETESTSTRUCT_HPP__
