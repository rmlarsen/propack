#
#     (C) Rasmus Munk Larsen, Stanford University, 2004
#
include ../make.inc
IPATH = -I.

all: lib

### Object files.

OBJ = iblasext.o printstat.o

.f.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.f
.F.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.F
.c.o:
	$(CC)  $(CFLAGS) $(IPATH) -c -o $*.o $*.c


lib: libcommon_$(PLAT).a

libcommon_$(PLAT).a : $(OBJ)
	rm -f libcommon_$(PLAT).a
	$(MAKELIB) libcommon_$(PLAT).a $(OBJ) 
	$(RANLIB)  libcommon_$(PLAT).a

clean:	
	rm -f  *.o *.il libcommon_$(PLAT).a

cleanall:	
	rm -f  *.o  *.a *~ *.il
	rm -rf rii_files


