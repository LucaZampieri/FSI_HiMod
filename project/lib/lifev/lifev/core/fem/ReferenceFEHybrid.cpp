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
    @file
    @brief Reference finite element for Hdiv space.

    @author Alessio Fumagalli
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10-05-2010

    @contributor
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/ReferenceFEHybrid.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

// Costructor.
ReferenceFEHybrid::ReferenceFEHybrid ( std::string name, FE_TYPE type, ReferenceShapes shape, UInt nbDofPerVertex, UInt nbDofPerEdge,
                                       UInt nbDofPerFace, UInt nbDofPerVolume, UInt nbDof, UInt nbLocalCoor, const UInt& numberBoundaryFE,
                                       const CurrentFEManifold** boundaryFEList, const Real* refCoor, DofPatternType patternType ) :
    ReferenceFE ( name, type, shape, nbDofPerVertex, nbDofPerEdge, nbDofPerFace, nbDofPerVolume,
                  nbDof, nbLocalCoor, 1, static_cast<function_Type*> (NULL),  static_cast<function_Type*> (NULL),
                  static_cast<function_Type*> (NULL),  static_cast<function_Type*> (NULL), refCoor,
                  patternType, static_cast<ReferenceFE*> (NULL) ),
    M_numberBoundaryFE ( numberBoundaryFE ), M_boundaryFEList ( boundaryFEList )
{}

// Destructor.
ReferenceFEHybrid::~ReferenceFEHybrid()
{}



} // Namespace LifeV
