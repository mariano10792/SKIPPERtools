ifdef FITSIOROOT
    CFITSIO = $(FITSIOROOT)
endif
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags) -O3
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++14
else
  CFLAGS += -std=c++11
  CFLAGS += -std=c++0x
endif
LIBS = -lcfitsio -lm $(shell root-config --libs)  
ifdef CFITSIO
    LIBS += -L$(CFITSIO)/lib
endif
GLIBS = 
GLIBS += 
OBJECTS = DM.o tinyxml2.o gConfig.o gCal.o globalConstants.o 
HEADERS = globalConstants.h tinyxml2.h gConfig.h gCal.h

ALL : DM.exe
	@echo "Listo!"

DM.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o DM.exe $(LIBS) $(GLIBS) $(CFLAGS)

DM.o : DM.cc $(HEADERS)
	$(CPP) -c DM.cc -o DM.o $(CFLAGS)

gConfig.o : gConfig.cc $(HEADERS)
	$(CPP) -c gConfig.cc -o gConfig.o $(CFLAGS)

gCal.o : gCal.cc $(HEADERS)
	$(CPP) -c gCal.cc -o gCal.o $(CFLAGS)

tinyxml2.o : tinyxml2.cpp $(HEADERS)
	$(CPP) -c tinyxml2.cpp -o tinyxml2.o $(CFLAGS)

globalConstants.o : globalConstants.cc $(HEADERS)
	$(CPP) -c globalConstants.cc -o globalConstants.o $(CFLAGS)
	
clean:
	rm -f *~ *.o *.exe
