CC=gcc

IDIR=include
SDIR=src
ODIR=$(SDIR)/obj
LDIR=lib
BDIR=bin

CFLAGS=-I$(IDIR)

LIBS=""

_DEPS = mem.h matrix.h factorization.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = matrix.o mem.o factorization.o linalg.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/linalg: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

test:
	./$(BDIR)/linalg ge -v


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
