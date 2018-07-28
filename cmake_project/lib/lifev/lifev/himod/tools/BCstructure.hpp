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
    handles the enforcement of boundary conditions in in/outlet
    on the reference cylinder for the Navier-Stokes Equations

     @date 03/2014
     @author S. Guzzetti <sofia.guzzetti@gmail.com>

 */

#ifndef __BCHANDLER_HPP__
#define __BCHANDLER_HPP__

namespace LifeV
{

//! ReferenceMap - Class representing the function which maps the physical domain into a cylinder with section of unit radius.
/*!
 *  @author S. Guzzetti <sofia.guzzetti@gmail.com>
 *  This class handles the enforcement of boundary conditions in in/outlet on the reference cylinder for the Navier-Stokes Equations.
 *
 */

enum BCname { dir, neu, rob };

class BCScalarData
{
	public:
			typedef boost::function< Real (const Real& , const Real&, const Real&, const Real&, const UInt&) >	function_Type;

			BCScalarData(){
							std::cout << "Warning: empty Boundary Condition. Don't forget to add it through the setMethods!" << std::endl;
							};
			
			BCScalarData( const BCname& nameType, const function_Type& nameFunc ):
						M_BCType( nameType ), M_BCData( nameFunc ) {};
						
			~BCScalarData(){};
			
			void setBCType( const BCname& nameType );
			void setBCValue( const function_Type& nameFunc );
			
			BCname BCType() const
			{
				return M_BCType;
			}
			
			function_Type BCData() const
			{
				return M_BCData;
			}
			
	private:
			BCname		M_BCType;
			function_Type	M_BCData;
};

class BCdata
{
	public:
			typedef boost::function< Real (const Real& , const Real&, const Real&, const Real&, const UInt&) >	function_Type;

			BCdata(){
						std::cout << "Warning: empty Boundary Conditions. Don't forget to add them through the setMethods!" << std::endl;
					};
					
			BCdata( const BCScalarData& xIn, const BCScalarData& rIn, const BCScalarData& thetaIn,
					const BCScalarData& xOut, const BCScalarData& rOut, const BCScalarData& thetaOut ):
					M_xInflow( xIn ), M_rInflow( rIn ), M_thetaInflow( thetaIn ), M_xOutflow( xOut ), M_rOutflow( rOut ), M_thetaOutflow( thetaOut ) {};
			
			~BCdata(){};
			
			void setXVelocityInflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_xInflow.setBCType( nameType );
				M_xInflow.setBCValue( nameFunc );
				return;
			};
			
			void setRVelocityInflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_rInflow.setBCType( nameType );
				M_rInflow.setBCValue( nameFunc );
				return;
			};
			
			void setThetaVelocityInflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_thetaInflow.setBCType( nameType );
				M_thetaInflow.setBCValue( nameFunc );
				return;
			};
			
			void setXVelocityOutflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_xOutflow.setBCType( nameType );
				M_xOutflow.setBCValue( nameFunc );
				return;
			};
			
			void setRVelocityOutflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_rOutflow.setBCType( nameType );
				M_rOutflow.setBCValue( nameFunc );
				return;
			};
			
			void setThetaVelocityOutflow( const BCname& nameType, const function_Type& nameFunc )
			{
				M_thetaOutflow.setBCType( nameType );
				M_thetaOutflow.setBCValue( nameFunc );
				return;
			};
			
			BCScalarData xInflow() const
			{
				return M_xInflow;
			};
			
			BCScalarData rInflow() const
			{
				return M_rInflow;
			};
			
			BCScalarData thetaInflow() const
			{
				return M_thetaInflow;
			};
			
			BCScalarData xOutflow() const
			{
				return M_xOutflow;
			};
			
			BCScalarData rOutflow() const
			{
				return M_rOutflow;
			};
			
			BCScalarData thetaOutflow() const
			{
				return M_thetaOutflow;
			};
			
	private:
			BCScalarData		M_xInflow;
			BCScalarData		M_rInflow;
			BCScalarData		M_thetaInflow;
			
			BCScalarData		M_xOutflow;
			BCScalarData		M_rOutflow;
			BCScalarData		M_thetaOutflow;
};

} // namespace
#endif
