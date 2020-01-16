DBG = 
ODIR = ../bin_summit/
CC = mpicc
FORT = mpif90
LIB = -L/sw/summit/xl/16.1.1-1/xlf/16.1.1/lib

all:
	make param
	make awpdata
	make xapiir
	make peak

param:
	$(CC) -c $(DBG) fd3dparam.c

awpdata:
	$(CC) -c $(DBG) awp_data.c -O3

peak:
	$(CC) $(LIB) $(INC) -o $(ODIR)/pgv2 pgv2.c fd3dparam.o awp_data.o xapiir.o -O3 -lxlf90

extract:
	$(CC) $(DBG) -lm -DUSE_IN3D -o $(ODIR)extrts extrts.c awp_data.o  fd3dparam.o

xapiir:
	$(FORT) -c xapiir.f -O1
clean:
	rm -f *.o 
