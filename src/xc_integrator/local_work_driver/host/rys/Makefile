AR = ar
ARFLAGS = -rc

CC=gcc
CFLAGS=-Wall -O2 -lm -mavx

SRC=./src
INCLUDE=./include
TEST=./test

#DEBUG=
DEBUG=-DDEBUG

compile:
	$(CC) -c $(SRC)/boys_table.c -o $(SRC)/boys_table.o $(CFLAGS)
	$(CC) -c $(SRC)/jacobi_table.c -o $(SRC)/jacobi_table.o $(CFLAGS)

	$(CC) -c $(SRC)/rys_1rw.c -o $(SRC)/rys_1rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_2rw.c -o $(SRC)/rys_2rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_3rw.c -o $(SRC)/rys_3rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_4rw.c -o $(SRC)/rys_4rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_5rw.c -o $(SRC)/rys_5rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_xrw.c -o $(SRC)/rys_xrw.o $(CFLAGS)

	$(CC) -c $(SRC)/rys_rw.c -o $(SRC)/rys_rw.o $(CFLAGS)
	$(CC) -c $(SRC)/rys_integral.c -o $(SRC)/rys_integral.o $(CFLAGS) -I$(INCLUDE)

	$(AR) $(ARFLAGS) ./rys_integral.a $(SRC)/*.o

test: compile
	$(CC) -o $(TEST)/test_int_v0.x $(TEST)/test_int_v0.c ./rys_integral.a $(CFLAGS) -I$(INCLUDE) $(DEBUG)


clean:
	rm -rf ./*.a $(SRC)/*.o $(TEST)/*.x
