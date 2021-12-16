#LIBINT_ROOT = /global/cfs/cdirs/m1027/dbwy/mpqc4/li/install/haswell/release
#EIGEN_DIR   = /global/common/sw/cray/cnl7/haswell/eigen/3.3.7/gcc/8.2.0/2wwrykb/include/eigen3

CPU_CC = CC
CC = g++

#CONST_LIB = ../../../../../../include/
#LIBINT_ROOT = /home/dtpopovici/Executables/libint
#EIGEN_DIR   = /usr/local/include/eigen3

CONST_LIB = ../../../../../../include/
LIBINT_ROOT = /global/homes/t/thom13/Executables/libint_gnu/
EIGEN_DIR   = /global/homes/t/thom13/Executables/Eigen/include/eigen3/

compile:
	$(CC) test.cpp ../obara_saika.a $(LIBINT_ROOT)/lib/libint2.a -o test.x  -I$(CONST_LIB) -I$(LIBINT_ROOT)/include -I$(EIGEN_DIR) -I../include/ -std=c++17 -I/usr/common/software/sles15_cgpu/cuda/11.1.1/include/ -L/usr/common/software/sles15_cgpu/cuda/11.1.1/lib64/ -lcudart -lcudadevrt
