# Warning: dependences on .h files are not given. If you change
# these, type "make clean; make"

OBJS=check2xsf.o check_read.o xsf_write.o molecule_fix.o cube_write.o \
     cell_read.o periodic_table.o basis.o super.o xplor_write.o \
     pdb_write.o cell_write.o dx_write.o vasp_write.o pdb_read.o \
     xyz_write.o cml_write.o fdf_write.o gpfa.o rotate.o

# Don't define QSORT unless your libc provides qsort()

# Linux / IA32 / gcc
CFLAGS=-Wall -O -D_FILE_OFFSET_BITS=64 -DQSORT
# Tru64
#CFLAGS=-O -DQSORT
# Solaris
#CFLAGS=-O -xarch=native64

CC=cc

check2xsf: c2xsf.h $(OBJS)
	$(CC) $(CFLAGS) -o  check2xsf $(OBJS) -lm

clean:
	rm check2xsf $(OBJS)

