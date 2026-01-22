# project name (generate executable with this name)
TARGET   = chemins.exe
PYMOL = pymol

EXTRA_CFLAGS ?=
CFLAGS=-Wall -g -O2 $(EXTRA_CFLAGS) #-O3 #-march=native -mtune=native -flto -DNDEBUG
LDFLAGS=-lm
INCDIR=-Iinclude
OBJDIR=obj
SRCDIR=src
BINDIR=bin

CC=gcc
EXEC=clean dir $(BINDIR)/$(TARGET)
SRC:=$(wildcard $(SRCDIR)/*.c)
OBJ:=$(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: $(EXEC)

$(BINDIR)/$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJ) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(INCPATH) $(INCDIR) -o $@ -c $< $(CFLAGS)

dir:
	mkdir $(OBJDIR)
	mkdir $(BINDIR)

demo: all
	./$(BINDIR)/$(TARGET) -i demos_realSubstrat/ABINOS/ -n 0 -s 10

# Target to open PyMOL with the generated script
open_pymol:
	$(PYMOL) open_pymol.pml

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	ls

mrproper: clean
	rm -rf results