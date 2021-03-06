SET(fem_HEADERS
  fem/Assembly.hpp
  fem/AssemblyElemental.hpp
  fem/BCBase.hpp
  fem/BCDataInterpolator.hpp
  fem/BCFunction.hpp
  fem/BCHandler.hpp
  fem/BCIdentifier.hpp
  fem/BCManage.hpp
  fem/BCManageNormal.hpp
  fem/BCVector.hpp
  fem/CurrentFE.hpp
  fem/CurrentFEManifold.hpp
  fem/DOF.hpp
  fem/DOFInterface.hpp
  fem/DOFInterface3Dto2D.hpp
  fem/DOFInterface3Dto3D.hpp
  fem/DOFLocalPattern.hpp
  fem/FEField.hpp
  fem/FEFunction.hpp
  fem/FESpace.hpp
  fem/GeometricMap.hpp
  fem/HyperbolicFluxNumerical.hpp
  fem/PostProcessingBoundary.hpp
  fem/QuadraturePoint.hpp
  fem/QuadratureRule.hpp
  fem/QuadratureRuleProvider.hpp
  fem/ReferenceElement.hpp
  fem/ReferenceFE.hpp
  fem/ReferenceFEHdiv.hpp
  fem/ReferenceFEHybrid.hpp
  fem/ReferenceFEScalar.hpp
  fem/SobolevNorms.hpp
  fem/TimeAdvance.hpp
  fem/TimeAdvanceBDF.hpp
  fem/TimeAdvanceBDFVariableStep.hpp
  fem/TimeAdvanceData.hpp
  fem/TimeAdvanceNewmark.hpp
  fem/TimeData.hpp
CACHE INTERNAL "")

SET(fem_SOURCES
  fem/AssemblyElemental.cpp
  fem/BCBase.cpp
  fem/BCDataInterpolator.cpp
  fem/BCFunction.cpp
  fem/BCHandler.cpp
  fem/BCIdentifier.cpp
  fem/BCManage.cpp
  fem/BCVector.cpp
  fem/CurrentFE.cpp
  fem/CurrentFEManifold.cpp
  fem/DOF.cpp
  fem/DOFInterface.cpp
  fem/DOFInterface3Dto2D.cpp
  fem/DOFInterface3Dto3D.cpp
  fem/DOFLocalPattern.cpp
  fem/FEDefinitions.cpp
  fem/GeometricMap.cpp
  fem/QuadraturePoint.cpp
  fem/QuadratureRule.cpp
  fem/QuadratureRuleProvider.cpp
  fem/ReferenceElement.cpp
  fem/ReferenceFE.cpp
  fem/ReferenceFEHdiv.cpp
  fem/ReferenceFEHybrid.cpp
  fem/ReferenceFEScalar.cpp
  fem/TimeAdvanceData.cpp
  fem/TimeData.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fem_HEADERS})
  STRING(REPLACE "fem/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/fem"
  COMPONENT ${PACKAGE_NAME}
)
