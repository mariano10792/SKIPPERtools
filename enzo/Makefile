FITSIOROOT = /home/dario/CCD/Datos_Fe55/cfitsio/
CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS +=
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++11
else
  CFLAGS += -std=c++0x
endif
OBJECTS =  enzo.o 
#HEADERS = globalConstants.h

ALL : enzo.exe
	@echo "Listo!"

enzo.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o enzo.exe $(LIBS) $(GLIBS) $(CFLAGS)

enzo.o : enzo.C $(HEADERS)
	$(CPP) -c enzo.C -o enzo.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


