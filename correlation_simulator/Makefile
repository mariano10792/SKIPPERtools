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
OBJECTS =  CCD_simulation_correlation.o 
#HEADERS = globalConstants.h

ALL : CCD_simulation_correlation.exe
	@echo "Listo!"

CCD_simulation_correlation.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o CCD_simulation_correlation.exe $(LIBS) $(GLIBS) $(CFLAGS)

CCD_simulation_correlation.o : CCD_simulation_correlation.C $(HEADERS)
	$(CPP) -c CCD_simulation_correlation.C -o CCD_simulation_correlation.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


