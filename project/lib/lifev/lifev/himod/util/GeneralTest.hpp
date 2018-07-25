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
    @file generaltest.hpp
    @brief This file contains the description of ??

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi

 */
#ifndef __GENERALTEST_HPP_
#define __GENERALTEST_HPP_ 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>

#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/himod/util/CaseTestStruct.hpp>
#include <lifev/himod/util/CaseTest.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconversion"


#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/himod/modalbasis/HiModAssembler.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wconversion"


namespace LifeV
{

class GeneralTest
{
public:
    enum ListOfCaseTests
    {
        Camini,
        DDDD,
        DDDD_Adv,
        DDDD_ADR,
        DRDR,
        RRRR,
        RR2D,
        BDRR,
        BDDD,
        DD2D
    };
    //! typedef for the function type
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID& ) >      function_type;
    //! typedef for a vector in R^3
    typedef VectorSmall<3>                      TreDvector_type;
    //! typedef for the mesh
    typedef RegionMesh<LinearLine>              mesh_type;
    //! typedef for the matrix type
    typedef MatrixEpetraStructured<Real>        matrix_type;
    //! typedef for the vector type
    typedef VectorEpetraStructured              vector_type;

    typedef Preconditioner                      basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
    typedef PreconditionerIfpack                prec_Type;
    typedef boost::shared_ptr<prec_Type>        precPtr_Type;


    // Costruttori
    GeneralTest (    const QuadratureRule* quady = &quadRuleLobSeg32pt, const QuadratureRule* quadz = &quadRuleLobSeg32pt);

    // Methods
    void run();//******
    // Set Methods:
    void setDatafile (GetPot const& DataFile);

    void setCoefficients (Real const& mu, TreDvector_type const& beta, Real const& sigma,  Real const& chiy = 1, Real const& chiz = 1);
    void setExport (UInt const& nx, UInt const& ny, UInt const& nz);
    void setBC (std::string bcUp, std::string bcDown, std::string bcLeft, std::string bcRight);
    void setFunctions (function_type const& ForceTerm, function_type const& DirichletInflow);//*****
    void setData (   UInt const& m,//*******
                     UInt const& Nel);
    void setDomain (Real const& lx, Real const& ly, Real const& lz);

    void setCaseName (const std::string& CN)
    {
        M_caseName = CN;
    };

private:
    // Name of the test
    std::string M_caseName;

    // Data of the problem
    function_type M_ForceTerm;
    function_type M_DirichletInflow;

    //Coefficients
    Real M_mu;
    TreDvector_type M_beta;
    Real M_sigma;
    Real M_chiy;
    Real M_chiz;

    //BC nature
    std::string M_bcLeft;
    std::string M_bcRight;
    std::string M_bcDown;
    std::string M_bcUp;

    //Domain parameter
    Real M_lx;
    Real M_ly;
    Real M_lz;

    // Data for the test
    UInt M_m;
    UInt M_Nel;

    // Data for export
    UInt M_nx_grid;
    UInt M_ny_grid;
    UInt M_nz_grid;

    GetPot M_DataFile;

    std::map<std::string, ListOfCaseTests> M_CaseTestSwitch;

    const QuadratureRule* M_quady;
    const QuadratureRule* M_quadz;
}; //End of the class

} //End of lifeV namespace
#endif
