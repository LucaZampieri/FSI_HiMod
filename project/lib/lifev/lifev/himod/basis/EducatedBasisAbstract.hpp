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
 *   @file EducatedBasisAbstract.hpp
     @brief This is an abstract class that provides the generic interface
     that a 2-dimensional educated basis should have. It derives from Basis2DAbstract.

     @date 11/2013
     @author S. Guzzetti <sofia.guzzetti@gmail.com>


 */
#ifndef __EDUCATEDBASISABSTRACT_HPP__
#define __EDUCATEDBASISABSTRACT_HPP__

#include <lifev/himod/basis/Basis2DAbstract.hpp>

#include <lifev/himod/basis/EducatedBasis2DFunctorImplemented.hpp>

#include <lifev/navier_stokes/function/bessel/bessel.hpp>

namespace LifeV
{

class EducatedBasisAbstract: public Basis2DAbstract
{
public:
    //! typedef for the utility map for the eigenvalues
    typedef std::vector<EigenMap2D>                                EigenContainer;
       typedef std::vector<std::vector<bool> >                        spy_type;

    //! @name Constructor and Destructor
    //@{
    //! Constructor
    EducatedBasisAbstract(): Basis2DAbstract(), M_normTheta( 1. ) {};

    //! Destructor
    virtual ~EducatedBasisAbstract() {};
    //@}

    // setMethods
     // Inherited from Basis2DAbstract:
     // setRho
     // setTheta
     // setNumberModes

    //! @name Get Methods
    //@{ 
    
    //! chi
    //! This method returns one of the two parameters associated to Robin BC, if any.
    virtual Real chi() const = 0;
    
    //! eigenValues
    //! This method returns the vector containing M_mtot eigenvalues.
    virtual EigenContainer eigenValues() const
    {
        return M_eigenvalues;
    };
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
    void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const EigenContainer& eigenvalues,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta ) const
    {
        UInt dim = eigenvalues.size();
        phirho.resize( dim );
        dphirho.resize( dim );
        phitheta.resize( dim );
        dphitheta.resize( dim );

        for ( UInt i = 0; i != dim; ++i )
        {
            phirho[i].resize( quadrulerho->nbQuadPt() );
            dphirho[i].resize( quadrulerho->nbQuadPt() );
            phitheta[i].resize( quadruletheta->nbQuadPt() );
            dphitheta[i].resize( quadruletheta->nbQuadPt() );
        }

        Real  j0, j1, y0, y1, j0p, j1p, y0p, y1p;
        
        Real    x, jQuadraturePoint, djQuadraturePoint, yQuadraturePoint, dyQuadraturePoint;
     
        for ( UInt p = 0; p != dim; ++p )
        {
            for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
            {
                x = sqrt( M_eigenvalues[p].lambda ) * quadrulerho->quadPointCoor( n, 0 ) * M_Rho;

                if( M_eigenvalues[p].order == 0 )
                {
                    bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
                    jQuadraturePoint = j0;
                    djQuadraturePoint = j0p;
                    yQuadraturePoint = y0;
                    dyQuadraturePoint = y0p;
                }
                else if( M_eigenvalues[p].order == 1 )
                {
                    bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
                    jQuadraturePoint = j1;
                    djQuadraturePoint = j1p;
                    yQuadraturePoint = y1;
                    dyQuadraturePoint = y1p;
                }
                else
                {
                    int nm;
                    UInt dim = M_eigenvalues[p].order;
                    Real* jn = new Real[dim+1];
                    Real* yn = new Real[dim+1];
                    Real* jnp = new Real[dim+1];
                    Real* ynp = new Real[dim+1];
        
                    bessel::bessjyna( M_eigenvalues[p].order, x, nm, jn, yn, jnp, ynp );
                    assert( nm == M_eigenvalues[p].order );
                    
                    jQuadraturePoint = jn[nm];
                    djQuadraturePoint = jnp[nm];
                    yQuadraturePoint = yn[nm];
                    dyQuadraturePoint = ynp[nm];
                    
                    delete[] jn;
                    delete[] yn;
                    delete[] jnp;
                    delete[] ynp;
                }

                phirho[p][n]  = M_normRJ[p] * jQuadraturePoint;
                dphirho[p][n] = sqrt( M_eigenvalues[p].lambda ) * M_Rho * M_normRJ[p] * djQuadraturePoint;

            } // end n-for
            
            for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
            {
                phitheta[p][n]  = M_normTheta * ( sin ( M_eigenvalues[p].order * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) )
                                  +  cos ( M_eigenvalues[p].order * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                dphitheta[p][n] = M_normTheta * M_eigenvalues[p].order *
                                  ( cos ( M_eigenvalues[p].order * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) )
                                  - sin ( M_eigenvalues[p].order * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
            } // end n-for

        } // end p-for
        
        return;
}

    virtual void
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta
                        ) {}
    
    virtual Real 
    basisFunction( const UInt& k, const Real& rStar, const Real& tStar, 
                   const QuadratureRule* quadrulerho, MBMatrix_type phirho ) const = 0;

    //! evalSinglePoint
    //! This method evaluates the two-dimensional modal basis on an arbitrary point p.
    //! As the eigenvalue depends on the boundary condition, it is forced to be implemented in derived classes.
    /*!    @param    eigen        A reference to Real.
                            The eigenvalue associated to the modal function you want to compute.
        @param    index        A const UInt.
                            The index of the eigenvalue (0:M_mtot-1). It is needed beacause the coefficients 
                            of the modal function depend on it.
        @param    p            A reference to a const MBVector_type.
                            Polar coordinates (rho, theta) in (0,1)x(0,1) of the point of evaluation on the reference domain.
        @return                Real value of the evaluation of the (two-dimensional) 'index'-th modal function in point 'p'.
    */
    Real
    evalSinglePoint( const Real& eigen, const UInt index, const MBVector_type& p ) const
    {
        Real rhoh = p[0];
        Real thetah = p[1];
        Real x = sqrt( eigen ) * rhoh;
        Real  j0, j1, y0, y1, j0p, j1p, y0p, y1p;
        Real  jPoint, djPoint, yPoint, dyPoint;
     
        if( M_eigenvalues[index].order == 0 )
        {
             bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
             jPoint = j0;
             djPoint = j0p;
             yPoint = y0;
             dyPoint = y0p;
         }
         else if( M_eigenvalues[index].order == 1 )
         {
             bessel::bessjy01b( x, j0, j1, y0, y1, j0p, j1p, y0p, y1p );
             jPoint = j1;
             djPoint = j1p;
             yPoint = y1;
             dyPoint = y1p;
         }
         else
         {
              int nm;
             UInt dim = M_eigenvalues[index].order;
             Real* jn = new Real[dim+1];
             Real* yn = new Real[dim+1];
             Real* jnp = new Real[dim+1];
             Real* ynp = new Real[dim+1];
        
             bessel::bessjyna( M_eigenvalues[index].order, x, nm, jn, yn, jnp, ynp );
             assert( nm == M_eigenvalues[index].order );
                    
             jPoint = jn[nm];
             djPoint = jnp[nm];
             yPoint = yn[nm];
             dyPoint = ynp[nm];
                  
             delete[] jn;
             delete[] yn;
             delete[] jnp;
             delete[] ynp;
        }
        return M_normTheta * ( sin ( M_Theta * M_eigenvalues[index].order * thetah ) +  cos ( M_Theta * M_eigenvalues[index].order * thetah ) ) *
                ( M_normRJ[index] * ( jPoint ) );
    }


    virtual Real
    evalSinglePoint( const UInt index, const MBVector_type& p,
                       const QuadratureRule* quadrulerho, const MBMatrix_type& phirho ) const {}
    //@}
protected:

    //! normalization coefficient for goniometric basis
    Real M_normTheta;

    //! normalization coefficients for bessel function associated to J (one for each mode)
    MBVector_type M_normRJ;

    // some utilities
    void addRow( spy_type& matrix )
    {
        for( int i = 0; i != matrix.size(); ++i )
            matrix[i].push_back(0);
        return;
    }

    void addColumn( spy_type& matrix )
    {

        int dim = matrix[0].size();
        matrix.push_back( std::vector<bool> ( dim, 0 ) );

        return;
    }
    
    void deleteBesselRoots( MBMatrix_type& edges, MBMatrix_type& dedges )
    {
        for( UInt i = 0; i != edges.size(); ++i )
            edges[i].clear();
            
        for( UInt i = 0; i != dedges.size(); ++i )
            dedges[i].clear();

        return;
    }
};

typedef FactorySingleton< Factory < EducatedBasisAbstract,  std::string> > EducatedBasisFactory;

} // namespace LifeV

// Here the includes of object in the factory, we put them down here because the other options
// would have obliged the user to include all this header, that is too much complicated.
#include <lifev/himod/basis/EducatedBasisD.hpp>
#include <lifev/himod/basis/EducatedBasisR.hpp>
#include <lifev/himod/basis/EducatedBasisN.hpp>

#endif
