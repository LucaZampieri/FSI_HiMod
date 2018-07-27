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
 *   @file HiModExporterEnsight.hpp
     @brief This file contains the exporter for a HiMod solution.

     @date 05/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
 */

#ifndef __HIMODEXPORTERENSIGHT_HPP__
#define __HIMODEXPORTERENSIGHT_HPP__

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
#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

class HiModExporterEnsight
{
public:
	
	typedef	std::vector< std::vector<Real> >	matrix_pointType;
	typedef	std::vector< std::vector<UInt> >	matrix_tetraType;
	typedef	std::vector<Real>					vector_Type;
	typedef	std::vector<UInt>					vector_UInt_Type;
    typedef Basis2DAbstract                     basis2d_type;
    typedef boost::shared_ptr<basis2d_type>     basis2d_ptrType;

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
		The geometry must be already in the folder
     */
	HiModExporterEnsight( std::string fileName,
	                        const UInt& mx, const UInt& mr, const UInt& mtheta, const UInt& mp,
							const NSModalSpaceCircular& MB,
							const Real& h, const UInt& Nelements,
							const UInt& Ntimesteps, const Real& dt );
							
	~HiModExporterEnsight(){};
    //@}
    
    //! @name Methods
    //@{
    void writeSolution( const VectorEpetraStructured& sinCoeff, const VectorEpetraStructured& cosCoeff, const UInt& timeIter, const UInt& Nelements ) const;
    //@}
    
private:
	
	static void readGeometry();
	static void writeGeoFile();
	static void writeCaseFile( std::string& fileName, const UInt& Ntimesteps, const Real& dt );
	
	void interpolate( const VectorEpetraStructured& coeff,
	                  vector_Type& interpolation ) const;
	
	void evaluateModalBasis( matrix_pointType& xbasis,
	                         matrix_pointType& rbasis,
	                         matrix_pointType& tbasis,
	                         matrix_pointType& pbasis,
                             const int& trigonometricBasis ) const;
                             
	void evaluateSolution( const vector_Type& interp,
	                        const matrix_pointType& xbasis,
	                        const matrix_pointType& rbasis,
	                        const matrix_pointType& tbasis,
	                        const matrix_pointType& pbasis,
	                        vector_Type& gridSolution ) const;
	
	static matrix_pointType				M_points;
	static matrix_tetraType				M_tetrahedra;
		
	const Real						M_uMeshSize;
	const UInt                      M_Nelements;
	
	const UInt						M_mx;
	const UInt						M_mr;
	const UInt						M_mtheta;
	const UInt						M_mp;
	
	const QuadratureRule*           M_quadrulerho;
	const QuadratureRule*           M_quadruletheta;
	matrix_pointType				M_xRadialBasis;
	matrix_pointType				M_rRadialBasis;
	matrix_pointType				M_tRadialBasis;
	matrix_pointType				M_pRadialBasis;
	
	basis2d_ptrType    				M_xBasis;
	basis2d_ptrType    				M_rBasis;
	basis2d_ptrType    				M_tBasis;
	basis2d_ptrType    				M_pBasis;
};

} // end namespace


#endif
