CC=gcc
CFLAGS=-I.
DEPS = helpers.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

qr: qr.o functions.o
	$(CC) -o qr qr.o functions.o
