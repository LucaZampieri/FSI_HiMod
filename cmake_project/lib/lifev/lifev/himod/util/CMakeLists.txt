SET(util_HEADERS
  util/CheckModalBasis1D.hpp
  util/CheckModalBasis2D.hpp
  util/GeneralConvergenceTest.hpp
  util/ConvergenceErrorHandler.hpp
  util/CaseTest.hpp
  util/GeneralTest.hpp
CACHE INTERNAL "")

SET(util_SOURCES
   util/CheckModalBasis1D.cpp
   util/CheckModalBasis2D.cpp
   util/GeneralConvergenceTest.cpp
   util/ConvergenceErrorHandler.cpp
   util/CaseTest.cpp
   util/GeneralTest.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${util_HEADERS})
  STRING(REPLACE "util/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/himod/util"
  COMPONENT ${PACKAGE_NAME}
)
