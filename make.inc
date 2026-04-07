# some definitions controlling compilation
export F90=gfortran
# for HDF5 output, assuming MESASDK exists
export HDF5_HOME=$(MESASDK_ROOT)
export F90FLAGS=-g -O3 -Warray-bounds -fbounds-check -I$(HDF5_HOME)/include
export LDFLAGS=$(F90FLAGS)
#use the following line when compiling on unix
export LAPACKLIBS=-llapack  # -lg2c
# use the following line when compiling on mac
# export LAPACKLIBS=-framework Accelerate
export HDF5LIBS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -Wl,-rpath,$(HDF5_HOME)/lib
export LIBS+=$(LAPACKLIBS) $(HDF5LIBS)
