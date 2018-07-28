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
 *   @file Basis2DAbstract.hpp
     @brief This is an abstract class that provides the generic interface
     that a 2-dimensional basis should have. It derives from BasisAbstract.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>
     @author M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>


 */
#ifndef __BASIS2DABSTRACT_HPP__
#define __BASIS2DABSTRACT_HPP__

#include <lifev/himod/basis/BasisAbstract.hpp>

#include <cmath>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iomanip>

namespace LifeV
{

enum TrigonometricFunction { cosine, sine }; // sine=1, cosine=0.

struct EigenMap2D
{
    Real lambda;

    int order;
    int index;

    static EigenMap2D make_eigenmap2D (const Real& _lambda, const int& _k, const int& _n)
    {
        EigenMap2D a;

        a.lambda = _lambda;
        a.order = _n;
        a.index = _k;

        return a;
    }
};

struct Comparison2D
{
        bool operator() ( EigenMap2D const& a, EigenMap2D const& b ) const
        {
            if( a.lambda != b.lambda )
            {
                return ( a.lambda < b.lambda );
                }
            return a.order > b.order;
        }
};

class Basis2DAbstract: public BasisAbstract
{
public:
    //! typedef for the utility map for the eigenvalues
    typedef std::vector<EigenMap2D>                                EigenContainer;
       typedef std::vector<std::vector<bool> >                        spy_type;

    //! @name Constructor and Destructor
    //@{
    //! Constructor
    Basis2DAbstract(): M_Theta( 2*M_PI ) {};

    //! Destructor
    virtual ~Basis2DAbstract() {};
    //@}

    virtual UInt mtheta1D() const {};
    
    virtual UInt mr1D() const {};
    
    virtual boost::shared_ptr<matrix_Type> massMatrix() const {};

    //! @name setMethods
    //@{
    
    //! With this method you can set the radius of the transversal circular section of the domain. Altough this class is meant to
    // work on the reference domain, in some cases the radius of physical domain is still needed to adjust
    // the boundary conditions coefficients.
    void setRho( const Real& Rho )
    {
        M_Rho = Rho;
    };

    //! With this method you can set the range of the angle of the transversal circular section of the domain, which is set to [0,2*M_PI] by default.
    // Altough this class is meant to work on the reference domain, in some cases the radius of physical domain is still needed to adjust
    // the boundary conditions coefficients.
    void setTheta( const Real& Theta )
    {
        M_Theta = Theta;
    };

    //! With this method you can set the number of modal functions you need to approximate the solution.
       void setNumberModes( const UInt& m )
    {
        M_mtot = m;
    }    
    //@}
    
    //! @name Get Methods
    //@{ 
   
    //! eigenValues
    //! This method returns the vector containing M_mtot eigenvalues.
    virtual EigenContainer eigenValues() const = 0;
    //@}
    
    //! @name Methods
    //@{
    
    //! computeEigenvalues
    //! This method computes the eigenvalues of the bidimensional problem for boundary conditions.
    //! Eigenvalues depend directly on the kind of boundary condition, in fact the implementation
    //! is delegated to derived classes.
    virtual void computeEigenvalues() = 0; 
    
    //! evaluateBasis
    //! This is method evaluates the (two-dimensional) modalbasis defined on 
    //! quadrature nodes, as well as its derivative.
    /*! @param    phirho            A reference to a MBMatrix_type.
                                phirho[p][q] stores the evaluation of the Bessel function associated to
                                the p-th eigenvalue on q-th quadrature node. As it is filled in runtime
                                by the method, it is actually an output.
        @param    dphirho            A reference to a MBMatrix_type.
                                dphirho[p][q] stores the evaluation of the derivative of the derivative of
                                the Bessel function associated to the p-th eigenvalue on q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    phitheta        A reference to a MBMatrix_type.
                                phitheta[p][q] stores the evaluation of the angular function (which is independent
                                from the eigenvalue) in the q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    dphitheta        A reference to a MBMatrix_type.
                                dphitheta[p][q] stores the evaluation of the derivative of the angular function
                                in the q-th quadrature node.
                                As it is filled in runtime by the method, it is actually an output.
        @param    eigenvalues        A reference to MBVector_type.
                                Vector of (two-dimensional) eigenvalues.
        @param    quadrulerho        A pointer to a QuadratureRule for rho-integration.
        @param    quadruletheta    A pointer to a QuadratureRule for theta-integration.
         @return                    void
    */
    virtual void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const EigenContainer& eigenvalues,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta ) const = 0;
                       
    virtual void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta
                        ) = 0;
                       
        virtual Real 
        basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                       const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const = 0;

    virtual Real
    evalSinglePoint( const Real& eigen, const UInt index, const MBVector_type& p ) const = 0;
    
    virtual Real
    evalSinglePoint( const UInt index, const MBVector_type& p,
                       const QuadratureRule* quadrulerho, const MBMatrix_type& phirho
                        ) const = 0;
    
protected:

    //! Vector of eigenvalues
    EigenContainer M_eigenvalues;

    //! The radius of the transversal section of the domain.
    Real M_Rho;
    
    //! The maximal angle of the section (2*M_PI by default)
    Real M_Theta;

    //! Number of modes (i.e. number of modal functions).
    UInt M_mtot;
    
    boost::shared_ptr<matrix_Type> M_massMatrix;

};


// ===================================================
// Macros
// ===================================================

//! create factory for ModalSpace
//typedef FactorySingleton< Factory < Basis2DAbstract,  std::string> > Basis2DFactory;

} // namespace LifeV

// Here the includes of object in the factory, we put them down here because the other options
// would have obliged the user to include all this header, that is too much complicated.
#include <lifev/himod/basis/EducatedBasisAbstract.hpp>
#include <lifev/himod/basis/UneducatedBasisAbstract.hpp>

#endif
