FC=gfortran -ffixed-line-length-none
CC=gcc
FFLAGS=-C -O -g 
CFLAGS=-O -m64
#if SAC library has been installed, uncomment the next two lines
#CFLAGS=$(optimize) -DSAC_LIB
#SACLIB=-L/usr/local/sac/lib -lsac -lsacio


PROG= ./multitel3
SRCS=sub.bodyw3.f \
	sub.readray.f \
	sub.interpolation.f \
	sub.cfft.f \
	sub.clear.f \
	sub.cnvr.f \
	sub.cnvrsh.f \
	sub.instg.f \
	sub.prod.f \
	sub.qf.f \
	sub.radp3.f \
	sub.refl.f \
	sub.reflsh.f 
OBJS=$(SRCS:%.f=%.o) sacio.o

SUBS = fft.o Complex.o sacio.o

all: multisyn multitel3

multisyn: multisyn.o $(SUBS) radiats.o futterman.o
	$(LINK.f) -o $@ $^ $(LIBS) -lm

multitel3: $(OBJS) teleseis3.f
	$(FC) $(FFLAGS) -o $(PROG) teleseis3.f $(OBJS) $(LIBS)
	$(CC) -o rayinformation main.c -lm 

clean:
	rm -f *.o $(OBJS) *~ core
