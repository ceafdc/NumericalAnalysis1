CC=cc
CFLAGS= -std=c11 -Wall


all: main.app

main.app: main.o matrix.o
	$(CC) $(CFLAGS) *.o -o main.app

matrix.o: matrix/matrix.c
	$(CC) $(CFLAGS) -c matrix/matrix.c

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -rf *.o
	rm -rf *.app
