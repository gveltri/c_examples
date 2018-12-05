IDIR=include
CC=gcc
CFLAGS=-I$(IDIR)

SDIR=src
ODIR=$(SDIR)/obj
LDIR=lib
BDIR=bin

LIBS=""

_DEPS = matrix.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = qr.o matrix.o
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
