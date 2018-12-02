CC=gcc
CFLAGS=-I.
DEPS = matrix.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

qr: qr.o matrix.o
	$(CC) -o qr qr.o matrix.o
