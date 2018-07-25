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
 *   @file
     @brief This file contains the definition of the IntegrateMatrixElement class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_MATRIX_ELEMENT_HPP
#define INTEGRATE_MATRIX_ELEMENT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETMatrixElemental.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{


//! The class to actually perform the loop over the elements to assemble a matrix
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a matrix and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (This convertion is done
  within a typedef).
 */
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
class IntegrateMatrixElement
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::S_fieldDim,
            SolutionSpaceType::S_fieldDim,
            MeshType::S_geoDimensions >::evaluation_Type  evaluation_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                            const QuadratureRule& quadrature,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                            const ExpressionType& expression);

    //! Copy constructor
    IntegrateMatrixElement (const IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>& integrator);

    //! Destructor
    ~IntegrateMatrixElement();

    //@}


    //! @name Operators
    //@{

    //! Operator wrapping the addTo method
    template <typename MatrixType>
    inline void operator>> (MatrixType& mat)
    {
        addTo (mat);
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename MatrixType>
    inline void operator>> (boost::shared_ptr<MatrixType> mat)
    {
        addTo (mat);
    }

    //@}


    //! @name Methods
    //@{

    //! Ouput method
    void check (std::ostream& out = std::cout);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.
     */
    template <typename MatrixType>
    void addTo (MatrixType& mat);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.

      Specialized for the case where the matrix is passed as a shared_ptr
     */
    template <typename MatrixType>
    inline void addTo (boost::shared_ptr<MatrixType> mat)
    {
        ASSERT (mat != 0, " Cannot assemble with an empty matrix");
        addTo (*mat);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    IntegrateMatrixElement();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature to be used
    QuadratureRule M_quadrature;

    // Shared pointer on the Spaces
    boost::shared_ptr<TestSpaceType> M_testSpace;
    boost::shared_ptr<SolutionSpaceType> M_solutionSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;


    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE;
    ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim>* M_testCFE;
    ETCurrentFE<SolutionSpaceType::S_spaceDim, SolutionSpaceType::S_fieldDim>* M_solutionCFE;

    ETMatrixElemental M_elementalMatrix;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>::
IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                        const QuadratureRule& quadrature,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                        const ExpressionType& expression)
    :   M_mesh (mesh),
        M_quadrature (quadrature),
        M_testSpace (testSpace),
        M_solutionSpace (solutionSpace),
        M_evaluation (expression),

        M_testCFE (new ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), quadrature) ),
        M_solutionCFE (new ETCurrentFE<SolutionSpaceType::S_spaceDim, SolutionSpaceType::S_fieldDim> (solutionSpace->refFE(), testSpace->geoMap(), quadrature) ),

        M_elementalMatrix (TestSpaceType::S_fieldDim * testSpace->refFE().nbDof(),
                           SolutionSpaceType::S_fieldDim * solutionSpace->refFE().nbDof() )
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case TRIANGLE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case QUAD:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case TETRA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case HEXA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
    M_evaluation.setTestCFE (M_testCFE);
    M_evaluation.setSolutionCFE (M_solutionCFE);
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>::
IntegrateMatrixElement (const IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_quadrature (integrator.M_quadrature),
        M_testSpace (integrator.M_testSpace),
        M_solutionSpace (integrator.M_solutionSpace),
        M_evaluation (integrator.M_evaluation),

        M_testCFE (new ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_quadrature) ),
        M_solutionCFE (new ETCurrentFE<SolutionSpaceType::S_spaceDim, SolutionSpaceType::S_fieldDim> (M_solutionSpace->refFE(), M_solutionSpace->geoMap(), M_quadrature) ),

        M_elementalMatrix (integrator.M_elementalMatrix)
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case TRIANGLE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case QUAD:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case TETRA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case HEXA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (M_quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
    M_evaluation.setTestCFE (M_testCFE);
    M_evaluation.setSolutionCFE (M_solutionCFE);
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>::
~IntegrateMatrixElement()
{
    delete M_globalCFE;
    delete M_testCFE;
    delete M_solutionCFE;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental matrix : " << std::endl;
    M_elementalMatrix.showMe (out);
    out << std::endl;
}


template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
template <typename MatrixType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>::
addTo (MatrixType& mat)
{
    UInt nbElements (M_mesh->numElements() );
    UInt nbQuadPt (M_quadrature.nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    UInt nbSolutionDof (M_solutionSpace->refFE().nbDof() );

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Zeros out the matrix
        M_elementalMatrix.zero();

        // Update the currentFEs
        M_globalCFE->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
        M_testCFE->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
        M_solutionCFE->update (M_mesh->element (iElement), evaluation_Type::S_solutionUpdateFlag);

        // Update the evaluation
        M_evaluation.update (iElement);

        // Loop on the blocks

        for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
        {
            for (UInt jblock (0); jblock < SolutionSpaceType::S_fieldDim; ++jblock)
            {

                // Set the row global indices in the local matrix
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalMatrix.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                // Set the column global indices in the local matrix
                for (UInt j (0); j < nbSolutionDof; ++j)
                {
                    M_elementalMatrix.setColumnIndex
                    (j + jblock * nbSolutionDof,
                     M_solutionSpace->dof().localToGlobalMap (iElement, j) + jblock * M_solutionSpace->dof().numTotalDof() );
                }

                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        for (UInt j (0); j < nbSolutionDof; ++j)
                        {
                            M_elementalMatrix.element (i + iblock * nbTestDof, j + jblock * nbSolutionDof) +=
                                M_evaluation.value_qij (iQuadPt, i + iblock * nbTestDof, j + jblock * nbSolutionDof)
                                * M_globalCFE->wDet (iQuadPt);

                        }
                    }
                }
            }
        }

        M_elementalMatrix.pushToGlobal (mat);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
