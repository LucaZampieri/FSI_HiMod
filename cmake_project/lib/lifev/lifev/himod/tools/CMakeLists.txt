SET(tools_HEADERS
  tools/BCstructure.hpp
  tools/ReferenceMap.hpp
  tools/utilityFunctions.hpp
  tools/Poiseuille.hpp
  tools/RealGeoData.hpp
  tools/utilityFunctionsCone.hpp
  tools/utilityFunctionsPipe.hpp
  tools/utilityUnsteadyFunctions.hpp
  tools/utilityUnsteadyPressureFunctions.hpp
  tools/utilityUnsteadyPressureFunctionsCone.hpp
  tools/utilityNSFunctions.hpp
  tools/HiModExporterVtk.hpp
  tools/HiModPostProcessor.hpp
CACHE INTERNAL "")

SET(tools_SOURCES
   tools/ReferenceMap.cpp
   tools/HiModExporterVtk.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${tools_HEADERS})
  STRING(REPLACE "tools/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/himod/tools"
  COMPONENT ${PACKAGE_NAME}
)
