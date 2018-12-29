CC=gcc

IDIR=include
SDIR=src
ODIR=$(SDIR)/obj
LDIR=lib
BDIR=bin

CFLAGS=-I$(IDIR) -g -Wall -Wextra

LIBS="-lm"

_DEPS = mem.h matrix.h factorization.h estimation.h precision.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =  mem.o matrix.o factorization.o estimation.o precision.o linalg.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/linalg: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

gdb:
	gdb $(BDIR)/linalg
