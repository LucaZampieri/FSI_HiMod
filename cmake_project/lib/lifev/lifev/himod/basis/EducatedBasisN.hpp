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
 *   @file EducatedBasisN.hpp
     @brief This file contains the interface of a two-dimensional educated basis
     		with Neumann boundary conditions.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
 */
 
#ifndef __EDUCATEDBASISN_HPP__
#define __EDUCATEDBASISN_HPP__

#include <lifev/himod/basis/EducatedBasisAbstract.hpp>

namespace LifeV
{

class EducatedBasisN : public EducatedBasisAbstract
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
    The only available constructor: this class should be used only through the factory so it needs only to be
    default-constructable.
    */
    EducatedBasisN() : EducatedBasisAbstract() {};
    //@}
    
    //! @name Methods
    //@{
	//! computeEigenvalues
	//! This method computes the eigenvalues of the bidimensional problem for Neumann boundary conditions.
    virtual void computeEigenvalues();

    //@}

    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const {};
	// getMethod inherited from Basis2DAbstract - override
    Real chi() const
    {
        return 0.0;
    };

/*	METHODS INHERITED FROM BASIS2DABSTRACT

	SETMETHODS
	setRho
	setTheta
	setMu
	setChi
	setNumberModes
	evaluateBasis
	computeEigenvalues
	
	GETMETHODS
	eigenvalues
*/

protected:
/* 	PROTECTED MEMBERS INHERITED FROM BASIS2DABSTRACT
	Real M_Rho;
	Real M_Theta;
    UInt M_mtot;
    MBVector_type M_eigenvalues;
    Real M_normTheta;
    MBVector_type M_normRJ;
*/
};

// ===================================================
// Macros
// ===================================================

inline
EducatedBasisAbstract* createN()
{
    return new EducatedBasisN() ;
}

namespace
{
static bool registerN = EducatedBasisFactory::instance().registerProduct( "neu",  &createN );
}

} //End LifeV namespace
#endif
