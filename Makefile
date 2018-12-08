CC=gcc

IDIR=include
SDIR=src
ODIR=$(SDIR)/obj
LDIR=lib
BDIR=bin

CFLAGS=-I$(IDIR)


LIBS=""

_DEPS = mem.h matrix.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = matrix.o mem.o qr.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/QR: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

run:
	./$(BDIR)/QR

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
