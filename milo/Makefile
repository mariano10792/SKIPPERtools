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
OBJECTS =  milo.o 
#HEADERS = globalConstants.h

ALL : milo.exe
	@echo "Listo!"

milo.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o milo.exe $(LIBS) $(GLIBS) $(CFLAGS)

milo.o : milo.C $(HEADERS)
	$(CPP) -c milo.C -o milo.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


