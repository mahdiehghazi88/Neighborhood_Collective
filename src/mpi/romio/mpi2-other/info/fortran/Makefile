CC          = gcc
AR          = ar
LIBNAME     = /home/mahdieh/project/mahdieh/mvapich_nbrhood11/mvapich2-2.2/lib/libmpi.a
RANLIB      = ranlib
srcdir      = /home/mahdieh/project/mahdieh/mvapich_nbrhood11/mvapich2-2.2/src/mpi/romio/mpi2-other/info/fortran

INCLUDE_DIR = -I. -I../../../include -I${srcdir}/../../../adio/include -I${srcdir}/../../../mpi-io -I../../../adio/include
CFLAGS      =    -DNDEBUG -DNVALGRIND -O2  -DFORTRANUNDERSCORE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_ROMIOCONF_H -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 $(MPIOPROFILE) $(INCLUDE_DIR)

VPATH=.:${srcdir}

MPIO_OBJECTS = info_createf.o info_dupf.o info_getf.o info_getnthf.o \
      info_setf.o info_deletef.o info_freef.o info_getnksf.o info_getvlnf.o

MPIO_TMP_POBJECTS = info_createf.p info_dupf.p info_getf.p info_getnthf.p \
      info_setf.p info_deletef.p info_freef.p info_getnksf.p info_getvlnf.p

MPIO_REAL_POBJECTS = _info_createf.o _info_dupf.o _info_getf.o _info_getnthf.o\
    _info_setf.o _info_deletef.o _info_freef.o _info_getnksf.o _info_getvlnf.o

all: $(LIBNAME)

.SUFFIXES: $(SUFFIXES) .p

.c.o:
	$(CC) $(CFLAGS) -c $<

.c.p:
	@cp $(srcdir)/$*.c _$*.c
	$(CC) $(CFLAGS) -c _$*.c
	@rm -f _$*.c

profile: $(MPIO_REAL_POBJECTS)
	$(AR) $(LIBNAME) $(MPIO_REAL_POBJECTS)
	$(RANLIB) $(LIBNAME) 
	@rm -f _*.o

$(LIBNAME): $(MPIO_OBJECTS)
	$(AR) $(LIBNAME) $(MPIO_OBJECTS)
	$(RANLIB) $(LIBNAME) 

clean:
	@rm -f *.o *.lo *.gcno *.gcda *.bb *.bbg
	@rm -f ${srcdir}/*.gcno ${srcdir}/*.gcda 
	@rm -f ${srcdir}/*.bb ${srcdir}/*.bbg

# Rules for the profiling objects
_info_createf.o: info_createf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_createf.c -o _info_createf.o
_info_dupf.o: info_dupf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_dupf.c -o _info_dupf.o
_info_getf.o: info_getf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_getf.c -o _info_getf.o
_info_getnthf.o: info_getnthf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_getnthf.c -o _info_getnthf.o
_info_setf.o: info_setf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_setf.c -o _info_setf.o
_info_deletef.o: info_deletef.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_deletef.c -o _info_deletef.o
_info_freef.o: info_freef.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_freef.c -o _info_freef.o
_info_getnksf.o: info_getnksf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_getnksf.c -o _info_getnksf.o
_info_getvlnf.o: info_getvlnf.c
	$(CC) $(CFLAGS) -c $(srcdir)/info_getvlnf.c -o _info_getvlnf.o
