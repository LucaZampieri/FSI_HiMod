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
    @file CheckErrorHandler.hpp
    @brief A class that is usefull to hold the error data in a convergence test.

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi@gmail.com>

 */
#ifndef __CONVERGENCEERRORHANDLER_HPP_
#define __CONVERGENCEERRORHANDLER_HPP_ 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>

#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
namespace LifeV
{
//  UTILITY FOR CONVERGE FILE WITH GNUPLOT ----------------------------------------
class ConvergenceErrorHandler
{
public:
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID& ) >      function_type;
    typedef VectorSmall<3>                      TreDvector_type;
    typedef RegionMesh<LinearLine>              mesh_type;
    typedef MatrixEpetraStructured<Real>        matrix_type;
    typedef VectorEpetraStructured              vector_type;

    /*!
        Make a file of output of the matrix of errors, we suppose that the matrix has this shape

       NspaceIT
        m1      err_m1_h1       err_m1_h2       err_m1_h3   ...
        m2      err_m2_h1       err_m2_h2       err_m2_h3   ...
        m3      err_m3_h1       err_m3_h2       err_m3_h3   ...
        m4      err_m4_h1       err_m4_h2       err_m4_h3   ...
        ..      ..          ..          ..
    */

    /*!
        This method initializes the matrix of the errors, according to the max number of modes, max numbero of different
        spatial discretization and the step of convergence.
    */

    ConvergenceErrorHandler (const UInt& m, const UInt& n, const UInt& mstart, const UInt& step);

    void convergeFile (std::string filename);


    /*
        Add error in the correct cell of the matrix ( mode , kind of spatial discretization )
    */
    void addError (const Real& err, const UInt& mod, const UInt& n);

    /*!
        Compute normL2_ues with higher precision (Nel elements in the 1d FESpace)
    */
    static Real computeNormUes (const UInt& Nel, const UInt& M, const function_type& ues, const Real& Lx, const Real& Ly, const Real& Lz,
                                const std::string& down,  const std::string& up,
                                const std::string& left,  const std::string& right);

private:
    std::vector< std::vector<Real> > M_error;
    UInt M_Nit;
};

} //End LifeV namespaces
#endif
