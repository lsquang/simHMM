

CFLAGS= -O3 -c -Wall 
gzFLAGS= -lstdc++ -lgzstream -lz

CXX      = g++   # for Linux RedHat 6.1, g++ version 2.95.2
AR       = ar cr
CPPFLAGS = -I. -O


all: hmm

hmm: libgzstream.a main.o CRec.o  CMCMC.o  CData.o COneSite.o COption.o CHMM.o
	g++ ${CPPFLAGS} -lpthread -O3 -Wall gzstream/gzstream.o build/CRec.o build/CMCMC.o build/CHMM.o build/CData.o build/main.o build/COneSite.o build/COption.o -lstdc++  -lz -lpthread -o hmm



main.o: main.cpp
	g++ ${CPPFLAGS}  -c -Wall main.cpp -lstdc++ -lgzstream -lz -lpthread -o build/main.o


CRec.o: CRec.cpp
	g++ -I. -c -Wall CRec.cpp -lstdc++ -lgzstream -lz -lpthread -o build/CRec.o



CMCMC.o: CMCMC.cpp
	g++ -I. -c -Wall CMCMC.cpp -lstdc++ -lgzstream -lz -lpthread -o build/CMCMC.o

CHMM.o: CHMM.cpp
	g++ -I. -c -Wall CHMM.cpp -lstdc++ -lgzstream -lz -lpthread -o build/CHMM.o

COneSite.o: COneSite.cpp
	g++ -I. -c -Wall COneSite.cpp -lstdc++ -lgzstream -lz -lpthread -o build/COneSite.o

COption.o: COption.cpp
	g++ -I. -c -Wall COption.cpp -lstdc++ -lgzstream -lz -lpthread -o build/COption.o


CData.o: CData.cpp
	g++ ${CPPFLAGS}  -c -Wall CData.cpp -lstdc++ -lgzstream -lz -lpthread -o build/CData.o

gzstream.o : gzstream/gzstream.C 
	${CXX} ${CPPFLAGS} -c -o gzstream/gzstream.o gzstream/gzstream.C

libgzstream.a : gzstream.o
	${AR} gzstream/libgzstream.a gzstream/gzstream.o

clean:
	rm -rf *.o */*.o hmm
