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
    @file ModalSpace.hpp
    @brief This file contains the description of the class that
    handle all the utilities related to the modal expansion
    on the transversal fiber for the Navier-Stokes Equations

     @date 03/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>

 */
#ifndef __NSMODALSPACEABSTRACT_HPP__
#define __NSMODALSPACEABSTRACT_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

//! ModalSpace - Class for Hi-Mod modal basis
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  Class representing the modalbasis.
 *
 */
class NSModalSpaceAbstract
{

public:

    //! @name Public Types
    //@{

    //! typedef for the matrix used to contain the computed value of the modal basis on the quadrature rule
    typedef std::vector<std::vector<Real> >              MBMatrix_type;

    //! typedef for the function type
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //!
    typedef std::vector<Real>                               MBVector_type;
    //@}

    //! @name Constructor & Destructor
    //@{
	NSModalSpaceAbstract( const UInt mx, const UInt mr, const UInt mtheta, const UInt mp ): M_mx( mx ), M_mr( mr ), M_mtheta( mtheta ), M_mp( mp ) {};

    //! Destroys the monodimensional basis generators
    virtual ~NSModalSpaceAbstract(){}
    //@}

    //! @name Methods
    //@{

    //! Calls EvaluateBasis() method own by the basis generators.
    virtual void evaluateBasis()=0;

    //! Compute The Fourier Coefficients of the function \f$g=g(y,z)\f$
    virtual    std::vector<Real> xFourierCoefficients( const function_Type& g, const Real& t ) const = 0;
    virtual    std::vector<Real> rFourierCoefficients( const function_Type& g, const Real& t ) const = 0;
    virtual    std::vector<Real> thetaFourierCoefficients( const function_Type& g, const Real& t ) const = 0;

    //! Compute The k-th Fourier coefficients of a function \f$f=f(x,y,z)\f$ in the point x
    virtual    Real xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const = 0;
    virtual    Real rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const = 0;
    virtual    Real thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const = 0;

/*    virtual Real compute_R11 (const UInt& j, const UInt& k,
                              const function_Type& mu, const Real& x) const = 0;

    virtual Real compute_R10 (const UInt& j, const UInt& k,
                              const function_Type& beta, const Real& x) const = 0;

    virtual Real compute_R00 (const UInt& j, const UInt& k,
                              const function_Type& mu, const function_Type& beta, const function_Type& sigma, const Real& x) const = 0;
*/

    //! ShowMe method
    virtual    void showMe() const = 0;

    //@}

    //! @name Get Methods
    //@{
    //! return the dimension of the modal basis
    const UInt mx() const
    {
        return M_mx;
    }
    
    const UInt mr() const
    {
        return M_mr;
    }

    const UInt mtheta() const
    {
        return M_mtheta;
    }

    const UInt mp() const
    {
        return M_mp;
    }
    //@}

protected:

    //! @name Private methods
    //@{
    //! This methods compute the list of the eigenvalues, it fills the M_eigenvalues
    virtual void eigensProvider()=0;

    //@}

    //! dimension of the modal basis
    const UInt M_mx;
    const UInt M_mr;
    const UInt M_mtheta;
    const UInt M_mp;

};
} // namespace LifeV
#endif // __MODALSPACEABSTRACT_HPP_
