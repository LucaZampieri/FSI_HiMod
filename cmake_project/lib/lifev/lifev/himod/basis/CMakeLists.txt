SET(basis_HEADERS
  basis/EducatedBasisFunctorAbstract.hpp
  basis/EducatedBasisFunctorImplemented.hpp
  basis/EducatedBasis2DFunctorImplemented.hpp
  basis/Basis1DAbstract.hpp
  basis/Basis2DAbstract.hpp
  basis/BasisAbstract.hpp
  basis/EducatedBasisAbstract.hpp
  basis/UneducatedBasisAbstract.hpp
  basis/EducatedBasisRR.hpp
  basis/EducatedBasisDD.hpp
  basis/EducatedBasisNN.hpp
  basis/EducatedBasisDR.hpp
  basis/EducatedBasisD.hpp
  basis/EducatedBasisN.hpp
  basis/EducatedBasisR.hpp
  basis/Robert.hpp
  basis/ChebyshevLinearShift.hpp
  basis/ChebyshevQuadraticShift.hpp
  basis/ChebyshevQuadShift.hpp
  basis/ChebyshevQuadraticShiftNatural.hpp
  basis/ChebyshevQuadShiftNatural.hpp
  basis/Bessel.hpp
  basis/Zernike.hpp
  basis/ZernikeNatural.hpp
  basis/FakeBasis.hpp
CACHE INTERNAL "")

SET(basis_SOURCES
  basis/EducatedBasisD.cpp
  basis/EducatedBasisN.cpp
  basis/EducatedBasisR.cpp
  basis/EducatedBasisRR.cpp
  basis/EducatedBasisDD.cpp
  basis/EducatedBasisDR.cpp
  basis/EducatedBasisNN.cpp
  basis/Robert.cpp
  basis/ChebyshevLinearShift.cpp
  basis/ChebyshevQuadraticShift.cpp
  basis/ChebyshevQuadShift.cpp
  basis/ChebyshevQuadraticShiftNatural.cpp
  basis/ChebyshevQuadShiftNatural.cpp
  basis/Bessel.cpp
  basis/Zernike.cpp
  basis/ZernikeNatural.cpp
  basis/FakeBasis.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${basis_HEADERS})
  STRING(REPLACE "basis/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/himod/basis"
  COMPONENT ${PACKAGE_NAME}
)
