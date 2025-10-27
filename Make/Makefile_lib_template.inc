#
#     (C) Rasmus Munk Larsen, Stanford University, 2004
#
include ../make.inc
IPATH = -I.

all: lib exampleprogs

### Object files.

REAL_OBJ_PREFIX = ../$(REAL_DIR)/$(REAL_TYPE_PREFIX)
REAL_OBJS = $(REAL_OBJ_PREFIX)blasext.o $(REAL_OBJ_PREFIX)ortho_rec.o \
            $(REAL_OBJ_PREFIX)bsvd.o  $(REAL_OBJ_PREFIX)gemm_ovwr.o
OBJ = $(TYPE_PREFIX)lanbpro.o $(TYPE_PREFIX)reorth.o $(TYPE_PREFIX)getu0.o $(TYPE_PREFIX)safescal.o \
       $(TYPE_PREFIX)blasext.o  $(TYPE_PREFIX)lansvd.o $(TYPE_PREFIX)gemm_ovwr.o \
       $(TYPE_PREFIX)lansvd_irl.o $(TYPE_PREFIX)ritzvec.o $(REAL_OBJS)

LIBFILE = lib$(TYPE_PREFIX)propack_$(PLAT).a 

.f.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.f
.F.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.F
.c.o:
	$(CC)  $(CFLAGS) $(IPATH) -c -o $*.o $*.c


lib: $(LIBFILE)
	@( cd ../common; \
	$(MAKE) lib; \
	cd  ../$(TYPE_DIR) )


$(LIBFILE) : $(OBJ)  $(REAL_OBJS)
	rm -f $(LIBFILE)
	$(MAKELIB) $(LIBFILE) $(OBJ) 
	$(RANLIB)  $(LIBFILE)

exampleprogs: lib
	@( cd Examples; \
	$(MAKE) all; \
	cd .. )

test: exampleprogs
	@( cd Examples; \
	$(MAKE) test; \
	cd .. )

verify: test
	@( cd Examples; \
	$(MAKE) verify; \
	cd .. )

clean:	
	rm -f  *.o *.il lib*_$(PLAT).a
	rm -rf rii_files
	@( cd Examples; \
	$(MAKE) clean; \
	cd .. )

cleanall:	
	rm -f  *.o  *.a *~ *.il
	rm -rf rii_files
	@( cd Examples; \
	$(MAKE) cleanall; \
	cd .. )

