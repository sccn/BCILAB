LFSRC=liblfev.c liblocf.c libmut.c libtube.c
LIBDIR=../lib
INSTALLDIR=/usr/local
LDFLAGS=-L$(LIBDIR) -Wl,-rpath $(INSTALLDIR)/lib -llfev -ltube -llocf -lmut
CFLAGS=-I. -I../include
all: mexlf mexpp
mexlf: mexlf.c mlfut.c
	mex $(MXFLAGS) $(CFLAGS) mexlf.c mlfut.c $(LDFLAGS)
mexpp: mexpp.c mlfut.c
	mex $(MXFLAGS) $(CFLAGS) mexpp.c mlfut.c $(LDFLAGS)
nodlls: mexlf.c mlfut.c mexpp.c $(LFSRC)
	mex $(MXFLAGS) mexlf.c mlfut.c $(LFSRC)
	mex $(MXFLAGS) mexpp.c mlfut.c $(LFSRC)
very-clean: clean
	rm -f mexlf.mexglx mexpp.mexglx
clean:
	rm -f *.o
FORCE:
	
