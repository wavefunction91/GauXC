AR = ar
ARFLGAS = -rc

CC = nvcc

SRC = ./src
INCLUDE = ./include
TEST = ./test

BOYS_FUNCTION = -I../../../../../../include/

CFLAGS = --ptxas-options=-v -arch=sm_70 -gencode=arch=compute_70,code=compute_70 -maxrregcount=128 -std=c++17 --expt-relaxed-constexpr -lm -I$(INCLUDE)

#DEBUG=
DEBUG = -DDEBUG

compile:
	$(CC) -c $(SRC)/chebyshev_boys_computation.cu -o $(SRC)/chebyshev_boys_computation.o $(CFLAGS) $(BOYS_FUNCTION) -I$(INCLUDE)

	$(CC) -c $(SRC)/integral_0.cu -o $(SRC)/integral_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_1.cu -o $(SRC)/integral_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2.cu -o $(SRC)/integral_2.o $(CFLAGS) $(BOYS_FUNCTION)

	$(CC) -c $(SRC)/integral_0_0.cu -o $(SRC)/integral_0_0.o $(CFLAGS) $(BOYS_FUNCTION) 
	$(CC) -c $(SRC)/integral_1_0.cu -o $(SRC)/integral_1_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_1_1.cu -o $(SRC)/integral_1_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_0.cu -o $(SRC)/integral_2_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_1.cu -o $(SRC)/integral_2_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_2.cu -o $(SRC)/integral_2_2.o $(CFLAGS) $(BOYS_FUNCTION)

	$(CC) -c $(SRC)/obara_saika_integrals.cu -o $(SRC)/obara_saika_integrals.o $(CFLAGS)

	$(AR) $(ARFLAGS) ./obara_saika.a $(SRC)/*.o

clean:
	rm -rf ./*.a $(SRC)/*.o
