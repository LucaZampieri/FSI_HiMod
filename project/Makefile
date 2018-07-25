CXX = mpicxx
CXXFLAGS = -std=c++11# -Wall

CPPFLAGS = -I. -I./lib/lifev -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/eigen/3.2.8/include/eigen3 -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/blacs/1.1/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hypre/2.11.0/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/suitesparse/4.5.1/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hdf5/1.8.16/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/base/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/metis/5/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scotch/6.0.4/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.60.0/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/include -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/tbb/4.4/include -I./lib/parmetis/include

LIBFLAGS = -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/lib -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scotch/6.0.4/lib -L./lib/lifev/static_lib -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/lib -llifevhimod -llifevnavierstokes -llifeveta -llifevcore -lmuparser -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/lib -llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lrythmos -lteko -lifpack2-adapters -lifpack2 -lamesos2 -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lshylu -lml -lgaleri-xpetra -lgaleri-epetra -lisorropia -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lifpack -lamesos -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lthyraepetra -lthyracore -lrtop -lbelostpetra -lbelosepetra -lbelos -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -lkokkosalgorithms -lkokkoscontainers -laztecoo -lzoltan -lepetraext -ltriutils -lepetra -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkoscore /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/glpk/4.60/lib/libglpk.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.60.0/lib/libboost_program_options.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.60.0/lib/libboost_system.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/matio/1.5.6/lib/libmatio.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/netcdf/4.4.0/lib/libnetcdf.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/suitesparse/4.5.1/lib/libcholmod.a /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libsmumps.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libdmumps.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libcmumps.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libzmumps.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libmumps_common.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib/libpord.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scalapack/2.0.2/lib/libscalapack.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/blacs/1.1/lib/libblacs.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/blacs/1.1/lib/libblacsF77init.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hypre/2.11.0/lib/libHYPRE.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/suitesparse/4.5.1/lib/libumfpack.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/suitesparse/4.5.1/lib/libamd.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hdf5/1.8.16/lib/libhdf5.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/tbb/4.4/lib/libtbb.so -lz /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/metis/5/lib/libparmetis.so /u/sw/pkgs/toolchains/gcc-glibc/5/prefix/lib/libgomp.so.1 /u/sw/pkgs/toolchains/gcc-glibc/5/base/lib/libmpi.so.12 /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/metis/5/lib/libmetis.so -lptscotch -lptscotcherr -lscotch -lscotcherr /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/lib/libopenblas.so -lhwloc ./lib/parmetis/lib/libparmetis.a /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/liblapack.so /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/libblas.so -Wl,-rpath,/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/lib:./lib/lifev/shared_lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/glpk/4.60/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.60.0/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/matio/1.5.6/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/netcdf/4.4.0/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scalapack/2.0.2/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/blacs/1.1/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hypre/2.11.0/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/suitesparse/4.5.1/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hdf5/1.8.16/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/tbb/4.4/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/metis/5/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scotch/6.0.4/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/lib:/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64


OBJS = build/FSIData.o build/FSISolver.o build/NSModalSpaceCircular.o build/ReferenceMap.o
HEADERS = include/FSIData.hpp include/FSISolver.hpp include/NSModalSpaceCircular.hpp include/NSHiModAssemblerCircular.hpp include/ReferenceMap.hpp
EXE = bin/main

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
	$(RM) bin/main output/*
