AR = ar
ARFLGAS = -rc

CC = g++

SRC = ./src
INCLUDE = ./include
TEST = ./test

BOYS_FUNCTION = -DX86_AVX -I../../../../../include/

CFLAGS = -Wall -O2 -lm -mavx -mfma -std=c++1z -I$(INCLUDE)

#DEBUG=
DEBUG = -DDEBUG

compile:
	$(CC) -c $(SRC)/chebyshev_boys_computation.cxx -o $(SRC)/chebyshev_boys_computation.o $(CFLAGS) $(BOYS_FUNCTION) -I$(INCLUDE)

	$(CC) -c $(SRC)/integral_0.cxx -o $(SRC)/integral_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_1.cxx -o $(SRC)/integral_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2.cxx -o $(SRC)/integral_2.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_3.cxx -o $(SRC)/integral_3.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4.cxx -o $(SRC)/integral_4.o $(CFLAGS) $(BOYS_FUNCTION)

	$(CC) -c $(SRC)/integral_0_0.cxx -o $(SRC)/integral_0_0.o $(CFLAGS) $(BOYS_FUNCTION) 
	$(CC) -c $(SRC)/integral_1_0.cxx -o $(SRC)/integral_1_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_1_1.cxx -o $(SRC)/integral_1_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_0.cxx -o $(SRC)/integral_2_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_1.cxx -o $(SRC)/integral_2_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_2_2.cxx -o $(SRC)/integral_2_2.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_3_0.cxx -o $(SRC)/integral_3_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_3_1.cxx -o $(SRC)/integral_3_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_3_2.cxx -o $(SRC)/integral_3_2.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_3_3.cxx -o $(SRC)/integral_3_3.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4_0.cxx -o $(SRC)/integral_4_0.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4_1.cxx -o $(SRC)/integral_4_1.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4_2.cxx -o $(SRC)/integral_4_2.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4_3.cxx -o $(SRC)/integral_4_3.o $(CFLAGS) $(BOYS_FUNCTION)
	$(CC) -c $(SRC)/integral_4_4.cxx -o $(SRC)/integral_4_4.o $(CFLAGS) $(BOYS_FUNCTION)

	$(CC) -c $(SRC)/obara_saika_integrals.cxx -o $(SRC)/obara_saika_integrals.o $(CFLAGS)

	$(AR) $(ARFLAGS) ./obara_saika.a $(SRC)/*.o

clean:
	rm -rf ./*.a $(SRC)/*.o
