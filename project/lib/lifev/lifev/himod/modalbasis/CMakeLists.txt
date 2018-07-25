SET(modalbasis_HEADERS
  modalbasis/ModalSpaceAbstract.hpp
  modalbasis/ModalSpaceCircular.hpp
  modalbasis/ModalSpaceRectangular.hpp
  modalbasis/HiModAssembler.hpp
  modalbasis/HiModAssemblerCircular.hpp
  modalbasis/HiModAssemblerRectangular.hpp
  modalbasis/Projector.hpp
  modalbasis/HMAaddNSAssembler.hpp
  modalbasis/HMAaddNSAssemblerPipe.hpp
  modalbasis/HMAaddNSAssemblerAxialDependent.hpp
  modalbasis/HMAbcHandler.hpp
  modalbasis/HMAbcHandlerPipe.hpp
  modalbasis/HMAinterpolation.hpp
  modalbasis/HMAinterpolationPipe.hpp
  modalbasis/HMAutilities.hpp
  modalbasis/HMAutilitiesPipe.hpp
  modalbasis/NSHiModAssembler.hpp
  modalbasis/NSHiModAssemblerPipe.hpp
  modalbasis/NSHiModAssemblerCircular.hpp
  modalbasis/NSModalSpaceAbstract.hpp
  modalbasis/NSModalSpaceCircular.hpp
  modalbasis/NSModalSpacePipe.hpp
CACHE INTERNAL "")

SET(modalbasis_SOURCES
  modalbasis/ModalSpaceCircular.cpp
  modalbasis/ModalSpaceRectangular.cpp
  modalbasis/NSModalSpaceCircular.cpp
  modalbasis/NSModalSpacePipe.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${modalbasis_HEADERS})
  STRING(REPLACE "modalbasis/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/himod/modalbasis"
  COMPONENT ${PACKAGE_NAME}
)
