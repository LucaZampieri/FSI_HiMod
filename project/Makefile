CXX = mpicxx
CXXFLAGS = -std=c++11# -Wall

MODULE_DIR = u/sw
PKGS_DIR = u/sw/pkgs/toolchains/gcc-glibc/5/pkgs

CPPFLAGS = -I. -I./lib/lifev \
 -I/$(PKGS_DIR)/muparser/2.2.5/include\
 -I/$(PKGS_DIR)/trilinos/12.6.3/include \
 -I/$(PKGS_DIR)/eigen/3.2.8/include/eigen3 \
 -I/$(PKGS_DIR)/mumps/5.0.1/include \
 -I/$(PKGS_DIR)/blacs/1.1/include \
 -I/$(PKGS_DIR)/hypre/2.11.0/include \
 -I/$(PKGS_DIR)/suitesparse/4.5.1/include \
 -I/$(PKGS_DIR)/hdf5/1.8.16/include \
 -I/$(LIB_DIR2)/pkgs/toolchains/gcc-glibc/5/base/include \
 -I/$(PKGS_DIR)/metis/5/include \
 -I/$(PKGS_DIR)/scotch/6.0.4/include \
 -I/$(PKGS_DIR)/boost/1.60.0/include \
 -I/$(PKGS_DIR)/openblas/0.2.17/include \
 -I/$(PKGS_DIR)/tbb/4.4/include \
 -I./lib/parmetis/include

LIBFLAGS = -L/$(PKGS_DIR)/trilinos/12.6.3/lib \
-L/$(PKGS_DIR)/scotch/6.0.4/lib \
-L./lib/lifev/static_lib \
-L/$(PKGS_DIR)/muparser/2.2.5/lib \
-llifevhimod -llifevnavierstokes -llifeveta -llifevcore -lmuparser \
-L/$(PKGS_DIR)/trilinos/12.6.3/lib \
-llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lrythmos -lteko -lifpack2-adapters -lifpack2 -lamesos2 \
-lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lshylu -lml -lgaleri-xpetra \
-lgaleri-epetra -lisorropia -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lifpack -lamesos -lanasazitpetra -lModeLaplace -lanasaziepetra \
-lanasazi -lthyraepetra -lthyracore -lrtop -lbelostpetra -lbelosepetra -lbelos -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels \
-ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -lkokkosalgorithms -lkokkoscontainers -laztecoo -lzoltan -lepetraext -ltriutils -lepetra \
-lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkoscore \
/$(PKGS_DIR)/glpk/4.60/lib/libglpk.so \
/$(PKGS_DIR)/boost/1.60.0/lib/libboost_program_options.so \
/$(PKGS_DIR)/boost/1.60.0/lib/libboost_system.so \
/$(PKGS_DIR)/matio/1.5.6/lib/libmatio.so \
/$(PKGS_DIR)/netcdf/4.4.0/lib/libnetcdf.so \
/$(PKGS_DIR)/suitesparse/4.5.1/lib/libcholmod.a \
/$(PKGS_DIR)/mumps/5.0.1/lib/libsmumps.so \
/$(PKGS_DIR)/mumps/5.0.1/lib/libdmumps.so \
/$(PKGS_DIR)/mumps/5.0.1/lib/libcmumps.so \
/$(PKGS_DIR)/mumps/5.0.1/lib/libzmumps.so \
/$(PKGS_DIR)/mumps/5.0.1/lib/libmumps_common.so \
/$(PKGS_DIR)/mumps/5.0.1/lib/libpord.so \
/$(PKGS_DIR)/scalapack/2.0.2/lib/libscalapack.so \
/$(PKGS_DIR)/blacs/1.1/lib/libblacs.so \
/$(PKGS_DIR)/blacs/1.1/lib/libblacsF77init.so \
/$(PKGS_DIR)/hypre/2.11.0/lib/libHYPRE.so \
/$(PKGS_DIR)/suitesparse/4.5.1/lib/libumfpack.so \
/$(PKGS_DIR)/suitesparse/4.5.1/lib/libamd.so \
/$(PKGS_DIR)/hdf5/1.8.16/lib/libhdf5.so \
/$(PKGS_DIR)/tbb/4.4/lib/libtbb.so -lz \
/$(PKGS_DIR)/metis/5/lib/libparmetis.so \
/$(MODULE_DIR)/pkgs/toolchains/gcc-glibc/5/prefix/lib/libgomp.so.1 \
/$(MODULE_DIR)/pkgs/toolchains/gcc-glibc/5/base/lib/libmpi.so.12 \
/$(PKGS_DIR)/metis/5/lib/libmetis.so \
-lptscotch -lptscotcherr -lscotch -lscotcherr \
/$(PKGS_DIR)/openblas/0.2.17/lib/libopenblas.so -lhwloc ./lib/parmetis/lib/libparmetis.a \
/$(PKGS_DIR)/lapack/3.6.0/lib64/liblapack.so \
/$(PKGS_DIR)/lapack/3.6.0/lib64/libblas.so -Wl,-rpath,/$(PKGS_DIR)/muparser/2.2.5/lib:$(PKGS_DIR)/trilinos/12.6.3/lib:./lib/lifev/shared_lib:/$(PKGS_DIR)/glpk/4.60/lib:/$(PKGS_DIR)/boost/1.60.0/lib:/$(PKGS_DIR)/matio/1.5.6/lib:/$(PKGS_DIR)/netcdf/4.4.0/lib:/$(PKGS_DIR)/mumps/5.0.1/lib:/$(PKGS_DIR)/scalapack/2.0.2/lib:/$(PKGS_DIR)/blacs/1.1/lib:/$(PKGS_DIR)/hypre/2.11.0/lib:/$(PKGS_DIR)/suitesparse/4.5.1/lib:/$(PKGS_DIR)/hdf5/1.8.16/lib:/$(PKGS_DIR)/tbb/4.4/lib:/$(PKGS_DIR)/metis/5/lib:/$(PKGS_DIR)/scotch/6.0.4/lib:/$(PKGS_DIR)/openblas/0.2.17/lib:/$(PKGS_DIR)/lapack/3.6.0/lib64


OBJS = build/FSIData.o build/FSISolver.o build/NSModalSpaceCircular.o build/ReferenceMap.o
HEADERS = include/FSIData.hpp include/FSISolver.hpp include/NSModalSpaceCircular.hpp include/NSHiModAssemblerCircular.hpp include/ReferenceMap.hpp
EXE = main

.PHONY : all clean distclean

all : $(EXE)

$(EXE) : $(OBJS) build/main.o build/FEDefinitions.o
	$(CXX) $(CXXFLAGS) $^ $(LIBFLAGS) -o $@

build/FEDefinitions.o : src/FEDefinitions.cpp include/QuadratureRule.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

build/main.o : src/main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJS) : build/%.o : src/%.cpp include/%.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean :
	$(RM) build/*

distclean : clean
	$(RM) $(EXE) output/*
