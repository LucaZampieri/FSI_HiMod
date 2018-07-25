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
    @file CaseTestStructCirc.hpp
    @brief A very simple structure to contain case test data

    @date 11/2013
    @author Massimiliano Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
    @author Sofia Guzzetti <sofia.guzzetti@gmail.com>

 */
#ifndef __CaseTestStructCirc_HPP__
#define __CaseTestStructCirc_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{

class CaseTestStructCirc
{
public:


    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID& ) >      function_type;


    CaseTestStructCirc (
        const std::string& casename,
        const Real& lx,
        const Real& rho,
        const Real& theta,
        const Real& Mu,
        const Real& beta1,

        const Real& beta2,

        const Real& beta3,
        const Real& Sigma,
        const Real& Chi,
        const std::string& bc,
        const function_type& F,
        const function_type& G,
        const function_type& Ues
    ) :
        CaseName (casename), Lx (lx), Rho (rho), Theta (theta), mu (Mu), sigma (Sigma), chi (Chi), bc (bc), f (F), g (G), ues (Ues)
    {
        beta[0] = beta1;
        beta[1] = beta2;
        beta[2] = beta3;
    };

    CaseTestStructCirc() {};

    CaseTestStructCirc (const CaseTestStructCirc& C)
    {
        CaseName = C.CaseName;
        Lx = C.Lx;
        Rho = C.Rho;
        Theta = C.Theta;
        mu = C.mu;
        beta = C.beta;
        sigma = C.sigma;
        chi = C.chi;
		bc = C.bc;
        f = C.f;
        g = C.g;
        ues = C.ues;
    };

    std::string CaseName;
    Real Lx;
    Real Rho;
    Real Theta;
    Real mu;
    VectorSmall<3> beta;
    Real sigma;
    Real chi;
    std::string bc;
    function_type f;
    function_type g;
    function_type ues;
};

} //End lifev namespace

#endif //End __CaseTestStructCirc_HPP__
