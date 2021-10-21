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
  CFLAGS += -std=c++14
else
  CFLAGS += -std=c++0x
endif
OBJECTS =  KS_SENSEI_test.o 
#HEADERS = globalConstants.h

ALL : KS_SENSEI_test.exe
	@echo "Listo!"

KS_SENSEI_test.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o KS_SENSEI_test.exe $(LIBS) $(GLIBS) $(CFLAGS)

KS_SENSEI_test.o : KS_SENSEI_test.C $(HEADERS)
	$(CPP) -c KS_SENSEI_test.C -o KS_SENSEI_test.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


