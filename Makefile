
CXX = g++
#CXXFLAGS = -O0 -g -std=c++11
CXXFLAGS = -O2 -Wall -march=native -std=c++11
INCLUDE = 
TLIB = -lm

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

.SUFFIXES: .cpp

.o :
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) 

#-----File Dependencies----------------------

SRC = $(SRC1) $(SRC2) $(SRC3)

SRC1 = chromosome.cpp dsmga2.cpp fastcounting.cpp global.cpp main.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp

SRC2 = chromosome.cpp dsmga2.cpp fastcounting.cpp global.cpp sweep.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp

SRC3 = genZobrist.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

OBJ1 = $(addsuffix .o, $(basename $(SRC1)))
OBJ2 = $(addsuffix .o, $(basename $(SRC2)))
OBJ3 = $(addsuffix .o, $(basename $(SRC3)))

all: DSMGA2 sweep genZobrist


DSMGA2: $(OBJ1)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ1)

sweep: $(OBJ2) DSMGA2
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ2)

genZobrist: $(OBJ3)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ3)

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ)

# DO NOT DELETE

chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h
main.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
main.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h zkey.h
main.o: sat.h trimatrix.h fastcounting.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h
nk-wa.o: nk-wa.h
sat.o: sat.h
chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h
sweep.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
sweep.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h
sweep.o: zkey.h sat.h trimatrix.h fastcounting.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h
nk-wa.o: nk-wa.h
sat.o: sat.h
