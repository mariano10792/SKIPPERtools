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
OBJECTS =  poisson_test_1erate.o 
#HEADERS = globalConstants.h

ALL : poisson_test_1erate.exe
	@echo "Listo!"

poisson_test_1erate.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o poisson_test_1erate.exe $(LIBS) $(GLIBS) $(CFLAGS)

poisson_test_1erate.o : poisson_test_1erate.C $(HEADERS)
	$(CPP) -c poisson_test_1erate.C -o poisson_test_1erate.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


