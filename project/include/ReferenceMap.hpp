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
    on the transversal circular fiber for the Navier-Stokes Equations

     @date 03/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>

 */

 // Search for " Luca " for part modified or added by him

#ifndef __REFERENCEMAP_HPP__
#define __REFERENCEMAP_HPP__

#include <lifev/core/LifeV.hpp>
#include "QuadratureRule.hpp"

#include <boost/shared_ptr.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MapVector.hpp>


namespace LifeV
{

class MapFunctor
{
    public:
            typedef std::vector<std::vector<Real> >              MBMatrix_type;
            typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

            // Constructor and destructor
             MapFunctor( const MBMatrix_type& M, const Real& h, const bool& rDependent=0 ) : M_data(M), M_h(h), M_scale(rDependent) {};
             MapFunctor(){};
            ~MapFunctor(){};

            // The function wrapped
            Real f( Real const& t, Real const& x, Real const& r, Real const& theta, ID const& j )
            {
                // REMARK. The input argument j is the index of the theta-node!!!
                /*if( fabs(remainder( x,M_h )) > 1e-8 )
                {
                    std::cout << "x = " << x << ", h = " << M_h << " x/h = " << x/M_h << "fmod = " << fmod(x,M_h) << std::endl;
                    std::cout << "ERROR: x is not a FEM node. ABORTED." << std::endl;
                    return 0;
                }
               */
                // If x and theta are HM nodes, just read the value
                if( fabs(remainder( x,M_h )) <= 1e-8 && fabs( theta-M_quadruleTheta->quadPointCoor(j,0)*2*M_PI ) <= 1e-8 )
                {
                    UInt dof( M_data.size() );
                    UInt index( std::floor(x/M_h) );
                    UInt i;
                    if( index%2 == 0 ) // Pressure DOF
                    {
                         i = index / 2;
                    }
                    else // Velocity DOF
                    {
                         i = (UInt)( (index + dof)/2 );
                    }
                    return M_data[i][j]*( 1 + M_scale*(r-1) );
                }
                else // otherwise interpolate
                {
                    UInt dof( M_data.size() );
                    UInt pdof( (dof-1)/2+1 );
                    UInt index( std::floor(x/M_h) );
                    UInt i;
                    UInt xIndL;
                    UInt xIndR;
                    Real xL;
                    Real xR;
                    if( index%2 == 0 ) // Pressure DOF
                    {
                         i = index / 2;
                         //xIndL = ( index==dof-1? i-1:i );
                         //xIndR = ( xIndL+1 );
                         xIndL = ( index==dof-1? dof-1:i );
                         xIndR = ( index==dof-1? pdof-1:xIndL+pdof );
                         //xL =  ( index==dof-1? (index-2)*M_h:index*M_h );
                         //xR = xL+2*M_h;
                         xL =  ( index==dof-1? (index-1)*M_h:index*M_h );
                         xR = xL+M_h;
                    }
                    else // Velocity DOF
                    {
                         i = (UInt)( (index + dof)/2 );
                         //xIndL = ( i==dof-1? i-1:i );
                         //xIndR = ( xIndL+1 );
                         xIndL = ( i );
                         xIndR = ( xIndL-pdof+1 );
                         //xL =  ( i==dof-1? (index-2)*M_h:index*M_h );
                         //xR = xL+2*M_h;
                         xL = ( index*M_h );
                         xR = xL+M_h;
                    }
                    UInt k(0);

                    while( k < M_quadruleTheta->nbQuadPt() && theta - 2*M_PI*M_quadruleTheta->quadPointCoor(k,0)>=0 )
                    {
                        ++k;
                    }
                    UInt indR( k );
                    UInt indL( k-1 );
                    if( k==M_quadruleTheta->nbQuadPt() )
                    {
                        indR = 0;
                        indL = k-1;
                    }
                    Real valueLL( M_data[xIndL][indL] );
                    Real valueLR( M_data[xIndL][indR] );
                    Real valueRL( M_data[xIndR][indL] );
                    Real valueRR( M_data[xIndR][indR] );
                    Real valueR( (x-xL)/(xR-xL)*(valueRR-valueLR) + valueLR );
                    Real valueL( (x-xL)/(xR-xL)*(valueRL-valueLL) + valueLL );
                    Real thetaL( 2*M_PI*M_quadruleTheta->quadPointCoor(indL,0) );
                    Real thetaR( 2*M_PI*M_quadruleTheta->quadPointCoor(indR,0) );

                    return (theta - thetaL)/( thetaR - thetaL )*( valueR-valueL ) + valueL;
                }
            };

            // Set methods
            void setFunction( const function_Type& func )
            {
                M_f = func;
                return;
            };
            void setData( const MBMatrix_type& data )
            {
                M_data = data;
                return;
            };
            void setH( const Real& h )
            {
                M_h = h;
                return;
            };
            void setScaling( const bool& s )
            {
                M_scale = s;
                return;
            };
            void setThetaQuadRule( const QuadratureRule& qr )
            {
                M_quadruleTheta = new QuadratureRule( qr );
                return;
            };

            // Get methods
            function_Type getFunction(){ return M_f; };

            const MBMatrix_type getData(){ return M_data; };

    protected:
            MBMatrix_type M_data;
            Real          M_h;
            function_Type       M_f;
            bool          M_scale;
            QuadratureRule*                        M_quadruleTheta;
};

//! ReferenceMap - Class representing the function which maps the physical domain into a cylinder with section of unit radius.
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  Class representing the function which maps the physical domain into a cylinder with section of unit radius.
 *
 */
class ReferenceMap
{
    public:
            // The following typedef may be useful in case of x-variant coefficients
            // typedef        vector< vector< vector< Real > > >        3dVector_type;
            //! Typedef for a function type
            typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

            //! typedef for the matrix used to contain the computed value of the modal basis on the quadrature rule
            typedef std::vector<std::vector<Real> >              MBMatrix_type;

            typedef std::vector<Real>                            MBVector_type;

            typedef VectorEpetraStructured                       vector_Type;

            typedef boost::shared_ptr<vector_Type>               vector_ptrType;

            ReferenceMap( const function_Type& jr,
                          const function_Type& jtheta,
                          const function_Type& dr,
                          const function_Type& dtheta,
                          const function_Type& drtheta,
                          const function_Type& dthetar,
                          const function_Type& jacobian,
                          const function_Type& jacobianWall, // Luca
                          const function_Type& inverseRhat,
                          const QuadratureRule* quadrho = &quadRuleSeg32pt,
                          const QuadratureRule* quadtheta = &quadRuleSeg32pt ) :
                        M_functionJr( jr ),
                        M_functionJtheta( jtheta ),
                        M_functionDr( dr ),
                        M_functionDtheta( dtheta ),
                        M_functionDrtheta( drtheta ),
                        M_functionDthetar( dthetar ),
                        M_functionJacobian( jacobian ),
                        M_functionJacobianWall( jacobianWall ), //Luca
                        M_inverseRhat( inverseRhat ),
                        M_quadruleRho( quadrho ),
                        M_quadruleTheta( quadtheta )
                        {
                            M_xJr = vector_ptrType();
                            M_xJtheta = vector_ptrType();
                            M_xDr = vector_ptrType();
                            M_xDtheta = vector_ptrType();
                            M_xDrtheta = vector_ptrType();
                            M_xDthetar = vector_ptrType();
                            M_xJacobian = vector_ptrType();
                            M_xJacobianWall = vector_ptrType(); // Luca
                            M_xR = vector_ptrType();
                            M_xdR = vector_ptrType();
                        };

            // Constructor for non-axisymmetric geometries
            ReferenceMap( const function_Type& jr, const function_Type& jtheta, const function_Type& dr, const function_Type& dtheta,
                        const function_Type& drtheta, const function_Type& dthetar, const function_Type& jacobian,
                        const function_Type& inverseRhat, const MapEpetra& epetraMap,
                        const QuadratureRule* quadrho = &quadRuleSeg32pt, const QuadratureRule* quadtheta = &quadRuleSeg32pt ) :
                        M_functionJr( jr ), M_functionJtheta( jtheta ), M_functionDr( dr ), M_functionDtheta( dtheta ),
                        M_functionDrtheta( drtheta ), M_functionDthetar( dthetar ), M_functionJacobian( jacobian ),
                        M_inverseRhat( inverseRhat ),
                        M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta )
                        {
                            M_xJr = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xJtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xDr = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xDtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xDrtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xDthetar = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xJacobian = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xR = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                            M_xdR = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
                        };

            // Constructor for patient-specific geometries
            ReferenceMap( const std::string& Rfile, const std::string& dxRfile, const std::string& dtRfile,
                          const MapEpetra& epetraMap, const UInt& femDof, const Real& h,
                          const QuadratureRule* quadrho = &quadRuleSeg32pt, const QuadratureRule* quadtheta = &quadRuleSeg32pt );

            ~ReferenceMap()
            {
                            M_xJr.reset();
                            M_xJtheta.reset();
                            M_xDr.reset();
                            M_xDtheta.reset();
                            M_xDrtheta.reset();
                            M_xDthetar.reset();
                            M_xJacobian.reset();
                            M_xJacobianWall.reset(); // Luca
                            M_xR.reset();
                            M_xdR.reset();
            };

            void evaluateMap( const Real& t );
            void evaluateAxialMap( const Real& h, const UInt& dof,
                                    const function_Type& Radius, const function_Type& dRadius );
            void evaluateFtensor( const Real& t, const Real& h, const UInt& dof );

            void readMatrix( const std::string fileName, const UInt& nR, const UInt& nC, MBMatrix_type& M );

            const MBMatrix_type Jr() const
            {
                return M_Jr;
            }

            const MBMatrix_type Jtheta() const
            {
                return M_Jtheta;
            }

            const MBMatrix_type Dr() const
            {
                return M_Dr;
            }

            const MBMatrix_type Dtheta() const
            {
                return M_Dtheta;
            }

            const MBMatrix_type Drtheta() const
            {
                return M_Drtheta;
            }

            const MBMatrix_type Dthetar() const
            {
                return M_Dthetar;
            }

            const MBMatrix_type Jacobian() const
            {
                return M_Jacobian;
            }

            const MBMatrix_type JacobianWall() const //Luca
            {
                return M_JacobianWall;
            }

            // ----------------     Axial dependence     ------------------
            const vector_Type xJr() const
            {
                return *M_xJr;
            }

            const vector_Type xJtheta() const
            {
                return *M_xJtheta;
            }

            const vector_Type xDr() const
            {
                return *M_xDr;
            }

            const vector_Type xDtheta() const
            {
                return *M_xDtheta;
            }

            const vector_Type xDrtheta() const
            {
                return *M_xDrtheta;
            }

            const vector_Type xDthetar() const
            {
                return *M_xDthetar;
            }

            const vector_Type xJacobian() const
            {
                return *M_xJacobian;
            }

            onst vector_Type xJacobianWall() const //Luca
            {
                return *M_xJacobianWall;
            }

            const vector_Type xR() const
            {
                return *M_xR;
            }

            const vector_Type xdR() const
            {
                return *M_xdR;
            }

            //! return the Quadrature rule
            const QuadratureRule&  qrRho() const
            {
                return *M_quadruleRho;
            }

            //! return the Quadrature rule
            const QuadratureRule&  qrTheta() const
            {
                return *M_quadruleTheta;
            }

            // Return functions
            const function_Type fR() const
            {
                return M_functionR;
            }

            const function_Type fDxR() const
            {
                return M_functionDxR;
            }

            const function_Type fDtR() const
            {
                return M_functionDtR;
            }

            const function_Type inverseRhat() const
            {
                return M_inverseRhat;
            }

            const function_Type fJacobian() const
            {
                return M_functionJacobian;
            }

            const function_Type fJacobianWall() const //Luca
            {
                return M_functionJacobianWall;
            }

            const function_Type fJr() const
            {
                return M_functionJr;
            }

            const function_Type fJtheta() const
            {
                return M_functionJtheta;
            }

            const function_Type fDr() const
            {
                return M_functionDr;
            }

            const function_Type fDtheta() const
            {
                return M_functionDtheta;
            }

            const function_Type fDrtheta() const
            {
                return M_functionDrtheta;
            }

            const function_Type fDthetar() const
            {
                return M_functionDthetar;
            }

    private:
            function_Type                          M_functionR;
            function_Type                          M_functionDxR;
            function_Type                          M_functionDtR;
            function_Type                          M_functionJr;
            function_Type                          M_functionJtheta;
            function_Type                          M_functionDr;
            function_Type                          M_functionDtheta;
            function_Type                          M_functionDrtheta;
            function_Type                          M_functionDthetar;
            function_Type                          M_functionJacobian;
            function_Type                          M_functionJacobianWall; // Luca

            function_Type                          M_inverseRhat;

            MapFunctor                             M_Rfunctor;
            MapFunctor                             M_RhatFunctor;
            MapFunctor                             M_Jfunctor;
            MapFunctor                             M_JfunctorWall; // Luca

            const QuadratureRule*                        M_quadruleRho;
            const QuadratureRule*                        M_quadruleTheta;

            MBMatrix_type                                M_R;
            MBMatrix_type                                M_dxR;
            MBMatrix_type                                M_dtR;
            MBMatrix_type                                M_Jr;
            MBMatrix_type                                M_Jtheta;
            MBMatrix_type                                M_Dr;
            MBMatrix_type                                M_Dtheta;
            MBMatrix_type                                M_Drtheta;
            MBMatrix_type                                M_Dthetar;
            MBMatrix_type                                M_Jacobian;
            MBMatrix_type                                M_JacobianWall; // Luca

            vector_ptrType                                M_xJr;
            vector_ptrType                                M_xJtheta;
            vector_ptrType                                M_xDr;
            vector_ptrType                                M_xDtheta;
            vector_ptrType                                M_xDrtheta;
            vector_ptrType                                M_xDthetar;
            vector_ptrType                                M_xJacobian;
            vector_ptrType                                M_xJacobianWall; // Luca
            vector_ptrType                                M_xR;
            vector_ptrType                                M_xdR;

// ---------------------------- New FSI ----------------------------------------
// All the section is Luca
     public:
       //! Evaluation of the map
       //void evaluateMapFSI( const Real& t );
// -----------------------------------------------------------------------------

};

} // namespace
#endif
