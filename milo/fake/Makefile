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
OBJECTS =  fakemilo.o 
#HEADERS = globalConstants.h

ALL : fakemilo.exe
	@echo "Listo!"

fakemilo.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o fakemilo.exe $(LIBS) $(GLIBS) $(CFLAGS)

fakemilo.o : fakemilo.C $(HEADERS)
	$(CPP) -c fakemilo.C -o fakemilo.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


