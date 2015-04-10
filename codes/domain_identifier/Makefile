CC = mpicc

GSL_INCL =  -I/home/sebastian/Software/gsl-1.16/include -I/home/sebastian/Software/gsl-1.16 -I/home/sebastian/Software/gsl-1.16/include/gsl
GSL_LIBS =  -L/home/sebastian/Software/gsl-1.16/.libs -L/home/sebastian/Software/gsl-1.16/cblas/.libs

CFLAGS = -c -O3 -I/export/$(USER)/local/include/ -I/usr/include/ $(GSL_INCL)
LFLAGS = -lm -L/export/$(USER)/local/lib -Wl,"-R /export/$(USER)/local/lib" 

#CFLAGS = $(OPTIONS) -c -O3 -I/lrz/sys/libraries/gsl/1.10/include
#LFLAGS = -L/lrz/sys/libraries/gsl/1.10/lib -Wl,"-R /export/$(USER)/local/lib"

PROGRAM = Domain_identifier

Domain_identifier:
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) $(GSL_LIBS) -lgsl $(GSL_LIBS) -lgslcblas -lm -o $@
	mv Domain_identifier Domain_identifier.x

clean:
	rm -rf $(PROGRAM)
	rm -rf *~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a	
	rm -rf *.so
	rm Domain_identifier.x
