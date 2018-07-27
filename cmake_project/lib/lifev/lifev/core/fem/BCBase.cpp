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
    @brief classes to handle boundary conditions

    @author M.A. Fernandez
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 06-2002

    @date 11-2002

  This file contains the classes which may be used to store boundary
  conditions. A boundary condition object will have the following
  elements:
<ol>
  <li> a name identifying a specific BC,

  <li> a flag identifying a specific part of the mesh boundary,

  <li> a type (Essential, Natural, Robin, Flux, Resistance),

  <li> a mode of implementation (Scalar, Full, Component, Normal,
     Tangential, Resistance, Directional),

  <li> a functor holding the data function,

  <li> a bool vector  describing the components involved in this boundary condition

  <li> a list of pointers to identifiers allowing the user to know to
     which DOF the boundary condition applies.
</ol>

 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCBase.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

BCBase::BCBase()
{
}

BCBase::BCBase ( const bcName_Type& name, const bcFlag_Type& flag,
                 const bcType_Type& type, const bcMode_Type& mode,
                 BCFunctionBase& bcFunction, const bcComponentsVec_Type& components )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components ( components ),
    M_bcFunction ( bcFunction.clone() ),
    M_bcFunctionFEVectorDependent(),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( -1 ),
    M_finalized ( false )
{
    if ( M_mode != Component )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase ( const bcName_Type& name,
                 const bcFlag_Type&  flag,
                 const bcType_Type&      type,
                 const bcMode_Type&      mode,
                 BCFunctionBase&    bcFunction ) :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction ( bcFunction.clone() ),
    M_bcFunctionFEVectorDependent(),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( -1 ),
    M_finalized ( false )
{
    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
        case Scalar:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( 0 );
            break;
        case Tangential:
            numberOfComponents = nDimensions - 1;
            M_components.reserve ( numberOfComponents );
            for ( ID i = 0; i < numberOfComponents; ++i )
            {
                M_components.push_back ( i );
            }
            break;
        case Normal:
            // Normal Essential boundary condition (cf Gwenol Grandperrin Master Thesis)
            if (type == Essential)
            {
                numberOfComponents = 1;
                M_components.reserve ( numberOfComponents );
                M_components.push_back ( nDimensions - 1 );
            }
            else
            {
                numberOfComponents = nDimensions;
                M_components.reserve ( numberOfComponents );
                for ( ID i = 0; i < numberOfComponents; ++i )
                {
                    M_components.push_back ( i );
                }
            }
            break;
        case Directional:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( nDimensions - 1 );
            break;
        default:
            ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase ( const bcName_Type& name,
                 const bcFlag_Type&  flag,
                 const bcType_Type&      type,
                 const bcMode_Type&      mode,
                 BCFunctionBase&    bcFunction,
                 const UInt&        numberOfComponents )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction ( bcFunction.clone() ),
    M_bcFunctionFEVectorDependent(),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( -1 ),
    M_finalized ( false )
{
    if ( M_mode != Full )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
    M_components.reserve ( numberOfComponents );
    for ( ID i = 0; i < numberOfComponents; ++i )
    {
        M_components.push_back ( i );
    }

}


BCBase::BCBase ( const bcName_Type& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector,
                 const bcComponentsVec_Type& components )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components (components),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent(),
    M_bcVector ( bcVector.clone() ),
    M_isStored_BcVector ( true ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( -1 ),
    M_finalized ( false )
{
    if ( mode != Component )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase ( const bcName_Type& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent(),
    M_bcVector ( bcVector.clone() ),
    M_isStored_BcVector ( true ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( 0 ), // The others are initialize to -1 we should follow a common convention.
    M_finalized ( false )
{
    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
        case Scalar:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( 0 );

            break;
        case Tangential:
            numberOfComponents = nDimensions - 1;
            M_components.reserve ( numberOfComponents );
            for ( ID i = 0; i < numberOfComponents; ++i )
            {
                M_components.push_back ( i );
            }

            break;
        case Normal:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( nDimensions - 1 );

            break;
        case Directional:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( nDimensions - 1 );
            break;
        default:
            ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}


BCBase::BCBase ( const bcName_Type& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector,
                 const UInt& numberOfComponents )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent(),
    M_bcVector ( bcVector.clone() ),
    M_isStored_BcVector ( true ),
    M_isStored_BcFunctionVectorDependent (false),
    M_offset ( -1 ),
    M_finalized ( false )
{
    if ( mode != Full )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    M_components.reserve ( numberOfComponents );
    for ( ID i = 0; i < numberOfComponents; ++i )
    {
        M_components.push_back ( i );
    }

}

BCBase::BCBase ( const bcName_Type&     name,
                 const bcFlag_Type&      flag,
                 const bcType_Type&          type,
                 const bcMode_Type&          mode,
                 BCFunctionUDepBase&    bcFunctionFEVectorDependent,
                 const bcComponentsVec_Type& components ) :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components ( components ),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent ( bcFunctionFEVectorDependent.clone() ),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (true),
    M_offset ( 0 ), // The others are initialize to -1 we should follow a common convention.
    M_finalized ( false )
{
    if ( M_mode != Component )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase ( const bcName_Type&  name,
                 const bcFlag_Type&   flag,
                 const bcType_Type&       type,
                 const bcMode_Type&       mode,
                 BCFunctionUDepBase& bcFunctionFEVectorDependent) :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent ( bcFunctionFEVectorDependent.clone() ),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (true),
    M_offset ( -1 ),
    M_finalized ( false )
{

    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
        case Scalar:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( 0 );
            break;
        case Tangential:
            numberOfComponents = nDimensions - 1;
            M_components.reserve ( numberOfComponents );
            for ( ID i = 0; i < numberOfComponents; ++i )
            {
                M_components.push_back ( i );
            }
            break;
        case Normal:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( nDimensions - 1 );
            break;
        case Directional:
            numberOfComponents = 1;
            M_components.reserve ( numberOfComponents );
            M_components.push_back ( nDimensions - 1);
            break;
        default:
            ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

}

BCBase::BCBase ( const bcName_Type&  name,
                 const bcFlag_Type&   flag,
                 const bcType_Type&       type,
                 const bcMode_Type&       mode,
                 BCFunctionUDepBase& bcFunctionFEVectorDependent,
                 const UInt&         numberOfComponents )
    :
    M_name ( name ),
    M_flag ( flag ),
    M_type ( type ),
    M_mode ( mode ),
    M_components(),
    M_bcFunction(),
    M_bcFunctionFEVectorDependent ( bcFunctionFEVectorDependent.clone() ),
    M_bcVector(),
    M_isStored_BcVector ( false ),
    M_isStored_BcFunctionVectorDependent (true),
    M_offset ( -1 ),
    M_finalized ( false )
{
    if ( M_mode != Full )
    {
        ERROR_MSG ( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    M_components.reserve ( numberOfComponents );
    for ( ID i = 0; i < numberOfComponents; ++i )
    {
        M_components.push_back ( i );
    }

}

BCBase::BCBase ( const BCBase& bcBase ) :
    M_name                                  ( bcBase.M_name ),
    M_flag                                  ( bcBase.M_flag ),
    M_type                                  ( bcBase.M_type ),
    M_mode                                  ( bcBase.M_mode ),
    M_components                            ( bcBase.M_components ),
    M_bcFunction                            ( ),
    M_bcFunctionFEVectorDependent           ( ),
    M_bcVector                              ( ),
    M_isStored_BcVector                     ( bcBase.M_isStored_BcVector ),
    M_isStored_BcFunctionVectorDependent    ( bcBase.M_isStored_BcFunctionVectorDependent ),
    M_idSet                                 ( ),
    M_idVector                              ( ),
    M_offset                                ( bcBase.M_offset ),
    M_finalized                             ( bcBase.M_finalized )
{
    // If the shared_ptr is not empty we make a true copy
    if ( bcBase.M_bcFunction.get() != 0 )
    {
        M_bcFunction = bcBase.M_bcFunction->clone();
    }
    if ( bcBase.M_bcFunctionFEVectorDependent.get() != 0 )
    {
        M_bcFunctionFEVectorDependent = bcBase.M_bcFunctionFEVectorDependent->clone();
    }
    if ( bcBase.M_bcVector.get() != 0 )
    {
        M_bcVector = bcBase.M_bcVector->clone();
    }

    // Important!!: The set member M_idSet is always empty at this point, it is just
    // an auxiliary container used at the moment of the boundary update (see BCHandler::bcUpdate)

    // The list of ID's must be empty
    if ( !M_idVector.empty() || !bcBase.M_idVector.empty() )
    {
        ERROR_MSG ( "BCBase::BCBase : The BC copy constructor does not work with list of identifiers which are not empty" );
    }
}


BCBase::~BCBase()
{
}

// ===================================================
// Methods
// ===================================================

ID BCBase::component ( const ID i ) const
{
    ASSERT_BD ( i < M_components.size() );
    return M_components[ i ];
}

bool  BCBase::isRobinCoeffAVector()  const
{
    if ( M_isStored_BcVector )
    {
        return  (*M_bcVector).isRobinCoeffAVector();
    }
    else
    {
        ERROR_MSG ( "BCBase::Robin : A data vector must be specified before calling this method" );
        return 0.;
    }
}

bool BCBase::isBetaCoeffAVector()   const
{
    if ( M_isStored_BcVector )
    {

        return   (*M_bcVector).isBetaCoeffAVector();
    }
    else
    {
        ERROR_MSG ( "BCBase::beta: A data vector must be specified before calling this method" );
        return 0.;
    }
}


Real BCBase::robinCoeffVector ( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector).robinCoeffVector ( iDof, iComponent );
    }
    else
    {
        ERROR_MSG ( "BCBase::RobinVec : A data vector must be specified before calling this method" );
        return 0.;
    }
}

Real BCBase::betaCoeffVector ( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector).betaCoeffVector ( iDof, iComponent );
    }
    else
    {
        ERROR_MSG ( "BCBase::RobinVec : A data vector must be specified before calling this method" );
        return 0.;
    }

}


const BCFunctionBase* BCBase::pointerToFunctor() const
{
    return M_bcFunction.get();
}

const BCFunctionUDepBase* BCBase::pointerToFunctorUDep() const
{
    return M_bcFunctionFEVectorDependent.get();
}

const BCVectorBase* BCBase::pointerToBCVector() const
{
    return M_bcVector.get();
}

void
BCBase::addBCIdentifier ( BCIdentifierBase* identifierToAddPtr )
{
    M_idSet.insert ( boost::shared_ptr<BCIdentifierBase> ( identifierToAddPtr ) );
}


UInt
BCBase::list_size() const
{
    return M_idVector.size();
}

std::ostream&
BCBase::showMe ( bool verbose, std::ostream& out ) const
{
    out << "********************************" << std::endl;
    out << "BC Name              : " << M_name << std::endl;
    out << "Flag                 : " << M_flag << std::endl;
    out << "Type                 : " << M_type << std::endl;
    out << "Mode                 : " << M_mode << std::endl;
    out << "Number of components : " << M_components.size() << std::endl;
    out << "List of components   : ";
    for ( ID i = 0; i < M_components.size(); ++i )
    {
        out << M_components[ i ] << " ";
    }
    out << std::endl;
    out << "Offset               : " << M_offset << std::endl;
    out << "Number of stored ID's: " << M_idVector.size() << std::endl;

    if ( verbose && M_finalized )
    {
        unsigned int count ( 0 ), lines ( 10 );
        out << "IDs in list";
        for ( std::vector<boost::shared_ptr<BCIdentifierBase> >::const_iterator i = M_idVector.begin();
                i != M_idVector.end(); ++i )
        {
            if ( count++ % lines == 0 )
            {
                out << std::endl;
            }
            out << ( *i ) ->id() << "  ";
        }
        if ( count % lines != 0 )
        {
            out << std::endl;
        }
        if ( isDataAVector() )
        {
            M_bcVector->showMe ( verbose, out );
        }
    }

    out << "********************************" << std::endl;
    return out;
}


// ===================================================
// Operators
// ===================================================

BCBase& BCBase::operator= ( const BCBase& BCb )
{

    M_name = BCb.M_name;
    M_flag = BCb.M_flag;
    M_type = BCb.M_type;
    M_mode = BCb.M_mode;
    M_finalized = BCb.M_finalized;
    M_isStored_BcVector = BCb.M_isStored_BcVector;
    M_bcFunctionFEVectorDependent = BCb.M_bcFunctionFEVectorDependent;
    M_bcFunction = BCb.M_bcFunction;
    M_bcVector = BCb.M_bcVector;
    M_isStored_BcVector = BCb.M_isStored_BcVector;
    M_isStored_BcFunctionVectorDependent = BCb.M_isStored_BcFunctionVectorDependent;
    M_offset  = BCb.M_offset;
    M_finalized = BCb.M_finalized;
    M_components = BCb.M_components;

    // Important!!: The set member M_idSet is always empty at this
    // point, it is just an auxiliary container used at the moment of
    // the boundary update (see BCHandler::bcUpdate)

    // The list of ID's must be empty
    if ( !M_idVector.empty() || !BCb.M_idVector.empty() )
    {
        ERROR_MSG ( "BCBase::operator= : The BC assigment operator does not work with lists of identifiers which are not empty" );
    }

    return *this;
}

const BCIdentifierBase*
BCBase::operator[] ( const ID& i ) const
{
    ASSERT_PRE ( M_finalized, "BC List should be finalized before being accessed" );
    ASSERT_BD ( i < M_idVector.size() );
    return M_idVector[ i ].get();
}

Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& iComponent ) const
{
    return M_bcFunction->operator() ( t, x, y, z, iComponent );
}

Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& iComponent, const Real& u ) const
{
    /* is there a better way ? */
    debugStream (800) << "debug800 in BCBase::operator(6x)\n";
    return M_bcFunctionFEVectorDependent->operator() (t, x, y, z, iComponent, u);
    debugStream (800) << "debug800 out BCBase::operator(6x)\n";
}

Real BCBase::operator() ( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector ) ( iDof, iComponent );
    }
    else
    {
        ERROR_MSG ( "BCBase::operator() : A data vector must be specified before calling this method" );
        return 0.;
    }
}


// ===================================================
// Set Methods
// ===================================================

void
BCBase::setBCVector ( const BCVectorBase& bcVector )
{
    M_bcVector = boost::shared_ptr<BCVectorBase > ( bcVector.clone() );
    M_isStored_BcVector = true;
    M_isStored_BcFunctionVectorDependent = false;
}

void
BCBase::setBCFunction ( const BCFunctionBase& bcFunction )
{
    M_bcFunction = bcFunction.clone();
    M_isStored_BcVector = false;
    M_isStored_BcFunctionVectorDependent = false;
}

void
BCBase::setBCFunction ( const BCFunctionUDepBase& bcFunctionFEVectorDependent )
{
    M_bcFunctionFEVectorDependent = bcFunctionFEVectorDependent.clone();
    M_isStored_BcVector = false;
    M_isStored_BcFunctionVectorDependent = true;
}


// ===================================================
// Get Methods
// ===================================================

std::string BCBase::name() const
{
    return M_name;
}

bcFlag_Type BCBase::flag() const
{
    return M_flag;
}

bcType_Type BCBase::type() const
{
    return M_type;
}

bcMode_Type BCBase::mode() const
{
    return M_mode;
}

UInt BCBase::numberOfComponents() const
{
    return M_components.size();
}

Real BCBase::robinCoeff() const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector ).robinCoeff();
    }
    else
    {
        ERROR_MSG ( "BCBase::robinCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

Real BCBase::resistanceCoeff() const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector ).resistanceCoeff();
    }
    else
    {
        ERROR_MSG ( "BCBase::resistanceCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

Real BCBase::betaCoeff() const
{
    if ( M_isStored_BcVector )
    {
        return ( *M_bcVector ).betaCoeff();
    }
    else
    {
        ERROR_MSG ( "BCBase::robinCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

bool BCBase::isDataAVector() const
{
    return M_isStored_BcVector;
}

bool BCBase::finalized() const
{
    return M_finalized;
}

bool BCBase::isUDep() const
{
    return M_isStored_BcFunctionVectorDependent;
}


// ===================================================
// Private Methods
// ===================================================

void
BCBase::copyIdSetIntoIdVector()
{
    if ( ! M_idSet.empty() )
    {
        M_idVector.clear();
        M_idVector.reserve ( M_idSet.size() );
        std::copy ( M_idSet.begin(), M_idSet.end(), std::inserter ( M_idVector, M_idVector.end() ) );
        M_idSet.clear();
    }
    M_finalized = true;
}


}
