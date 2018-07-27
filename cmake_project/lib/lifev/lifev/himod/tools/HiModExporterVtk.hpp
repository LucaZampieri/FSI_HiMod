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
 *   @file HiModExporterVtk.hpp
     @brief This file contains the exporter for a HiMod solution.

     @date 05/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
 */

#ifndef __HIMODEXPORTERVTK_HPP__
#define __HIMODEXPORTERVTK_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

#include <lifev/himod/basis/Basis2DAbstract.hpp>
#include <lifev/himod/modalbasis/NSModalSpaceCircular.hpp>
#include <lifev/himod/modalbasis/NSModalSpacePipe.hpp>
#include <lifev/himod/modalbasis/ModalSpaceCircular.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

class HiModExporterVtk
{
public:
    
    typedef    std::vector< std::vector<Real> >    matrix_pointType;
    typedef    std::vector< std::vector<UInt> >    matrix_tetraType;
    typedef    std::vector<Real>                   vector_Type;
    typedef    std::vector<UInt>                   vector_UInt_Type;
    typedef    Basis2DAbstract                     basis2d_type;
    typedef    boost::shared_ptr<basis2d_type>     basis2d_ptrType;
    typedef    boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
        Constant radius
     */
    HiModExporterVtk( const NSModalSpaceCircular& MB,
                      const Real& h, const UInt& Nelements ) :
                    M_uMeshSize( h ), M_Nelements( Nelements ),
                    M_R( MB.Rho() ), M_fR( fake ),
                    M_quadrulerho( MB.qrRho() ), M_quadruletheta( MB.qrTheta() ) {};
    /*!
        Non constant radius
     */                       
    HiModExporterVtk( const NSModalSpaceCircular& MB,
                      const Real& h, const UInt& Nelements, const bool xDependent ) :
                    M_uMeshSize( h ), M_Nelements( Nelements ),
                    M_R( 1. ), M_fR( MB.fRho() ),
                    M_quadrulerho( MB.qrRho() ), M_quadruletheta( MB.qrTheta() ) {};
    
    HiModExporterVtk( const NSModalSpacePipe& MB,
                      const Real& h, const UInt& Nelements ) :
                    M_uMeshSize( h ), M_Nelements( Nelements ),
                    M_R( 1. ), M_fR( MB.fRho() ),
                       M_quadrulerho( MB.qrRho() ), M_quadruletheta( MB.qrTheta() ) {};

    /*!
        Scalar unknown (ADR)
    */
    HiModExporterVtk( const ModalSpaceCircular& MB,
                        const Real& h, const UInt& Nelements ) :
                        M_uMeshSize( h ), M_Nelements( Nelements ),
                        M_R( MB.Rho() ),
                           M_quadrulerho( new QuadratureRule( MB.qrRho() ) ), M_quadruletheta( new QuadratureRule( MB.qrTheta() ) ) {};
                                                       
    ~HiModExporterVtk(){};
    //@}
    
    //! @name Methods
    //@{
    // R constant
    void writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D, const UInt& timeIter ) const;
    // R non constant
    void writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D, const UInt& timeIter, const bool xDependent ) const;
    
    // ADR
    void writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D ) const;
    //@}
    
private:
        
    const Real                        M_uMeshSize;
    const UInt                      M_Nelements;
    
    const Real                      M_R;
    const function_Type             M_fR;
    
    const QuadratureRule*           M_quadrulerho;
    const QuadratureRule*           M_quadruletheta;
    
    static Real fake( const Real& t, const Real& x, const Real& r, const Real& th, const ID& i ) { return 0; };
    
};

} // end namespace


#endif
