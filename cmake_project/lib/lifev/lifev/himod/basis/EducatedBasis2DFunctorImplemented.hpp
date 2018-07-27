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
    @file EducatedBasis2DFunctorImplemented.hpp
    @brief	This file contains all the different kinds of functor needed to provide all the different BCs for two-dimensional problem.
    		We put it all here beacuse they are very short and we preferred to not produce many different, very short, files.

    @date 11/2013
    @author S. Guzzetti <sofia.guzzetti@gmail.com>
    @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>

 */

#ifndef __EDUCATEDBASIS2DFUNCTORIMPLEMENTED_HPP__
#define __EDUCATEDBASIS2DFUNCTORIMPLEMENTED_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/navier_stokes/function/bessel/bessel.hpp>
#include <lifev/himod/basis/EducatedBasisFunctorAbstract.hpp>

namespace LifeV
{
class EducatedBasisFunctorR : public EducatedBasisFunctorAbstract
{
public:
	//! @name Constructor
	//@{    
    
    //! Constructor
    EducatedBasisFunctorR ( const Real& mu, const Real& chi, const Real& Rho, const UInt& n ) : EducatedBasisFunctorAbstract ( mu, chi, Rho ), M_order ( n ) {};
    //@}
    
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() ( const Real& x ) const
    {
		if ( M_order == 0 || M_order == 1 )
		{
		  double j0, j1, y0, y1, j0p, j1p, y0p, y1p;
		  
		  bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
		  
		  if ( M_order == 0 )
		  {
		  	return  x / M_L * ( j0p ) + M_chi / M_mu * ( j0 );
		  }
		  else
		  {
		  	return  x / M_L * ( j1p ) + M_chi / M_mu * ( j1 );
		  }
		}
		else
		{
		  int nm;
		  Real* jn = new Real[M_order+1];
		  Real* yn = new Real[M_order+1];
		  Real* jnp = new Real[M_order+1];
		  Real* ynp = new Real[M_order+1];
		  
		  bessel::bessjyna( M_order, x, nm, jn, yn, jnp, ynp );
//		  std::cout<<"nm = "<<nm<<", order = "<<M_order<<std::endl;
//		  assert( nm == M_order );
		  
		  Real tmp =  x / M_L * ( jnp[nm] ) + M_chi / M_mu * ( jn[nm] );
		  
		  delete[] jn;
		  delete[] yn;
		  delete[] jnp;
		  delete[] ynp;

		  return tmp;
		}
    }
    //@}
private:
   const UInt M_order;
};

/* This class is useless because the roots of Bessel functions derivatives are computed through an appropriate method.

class EducatedBasisFunctorN : public EducatedBasisFunctorAbstract
{
public:
	//! @name Constructor
	//@{
	// deletion of the constructor with missing input arguments
//    EducatedBasisFunctorN ( const Real& mu, const Real& chi, const Real& L ) = delete;    

    //! Constructor
    EducatedBasisFunctorN ( const Real& mu, const Real& chi, const Real& Rho, const UInt& n ) : EducatedBasisFunctorAbstract ( mu, chi, Rho ), M_order ( n ) {};
    //@}
    
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() ( const Real& x ) const
    {
	if ( M_order == 0 || M_order == 1 )
	{
	  double j0, j1, y0, y1, j0p, j1p, y0p, y1p;
	  
	  bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
	  
	  if ( M_order == 0 ) return j0p;
	  else return ( j1p );
	}
	else
	{
	  int nm;
	  Real* jn = new Real[M_order+1];
      Real* yn = new Real[M_order+1];
      Real* jnp = new Real[M_order+1];
      Real* ynp = new Real[M_order+1];
	  
	  bessel::bessjyna( M_order, x, nm, jn, yn, jnp, ynp );
      assert( nm == M_order );	  
	  Real tmp = ( jnp[nm] );
	  
	  delete[] jn;
      delete[] yn;
      delete[] jnp;
      delete[] ynp;
      
	  return tmp;
	}
    }
    //@}
private:
   const UInt M_order;
};
*/
/* This class is useless because the roots of Bessel functions are computed through an appropriate method.

class EducatedBasisFunctorD : public EducatedBasisFunctorAbstract
{
public:
	//! @name Constructor
	//@{
	// deletion of the constructor with missing input arguments
//    EducatedBasisFunctorD ( const Real& mu, const Real& chi, const Real& L ) = delete;   
    
    //! Constructor
    EducatedBasisFunctorD ( const Real& mu, const Real& chi, const Real& Rho, const UInt& n ) : EducatedBasisFunctorAbstract ( mu, chi, Rho ), M_order( n ) {};
    //@}
    
    //! @name Operators
    //@{
    //! The actual implementation of the functor
    Real operator() ( const Real& x ) const
    {
	if ( M_order == 0 || M_order == 1 )
	{
	  double j0, j1, y0, y1, j0p, j1p, y0p, y1p;
	  
	  bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
	  
	  if ( M_order == 0 ) return ( j0 );
	  else return ( j1 );
	}
	else
	{
	  int nm;
	  Real* jn = new Real[M_order+1];
      Real* yn = new Real[M_order+1];
      Real* jnp = new Real[M_order+1];
      Real* ynp = new Real[M_order+1];
	  
	  bessel::bessjyna( M_order, x, nm, jn, yn, jnp, ynp );
      assert( nm == M_order );
	  
	  Real tmp = ( jn[nm] );
	  
	  delete[] jn;
      delete[] yn;
      delete[] jnp;
      delete[] ynp;
      
	  return tmp;
	}
    }
    //@}
private:
   const UInt M_order;
};
*/
} // end namespace lifeV

#endif
