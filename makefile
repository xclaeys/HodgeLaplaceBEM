GCC  = g++ -O3
INCLUDE = -I./bemtool/

all: Formulation-N clean

######################################################

Formulation-N: Formulation-N.o
	$(GCC) Formulation-N.o -o Formulation-N

Formulation-N.o: Formulation-N.cpp
	$(GCC) $(INCLUDE) -c Formulation-N.cpp -o Formulation-N.o

######################################################

Formulation-T: Formulation-T.o
	$(GCC) Formulation-T.o -o Formulation-T

Formulation-T.o: Formulation-T.cpp
	$(GCC) $(INCLUDE) -c Formulation-T.cpp -o Formulation-T.o

######################################################

clean:
	rm *.o
