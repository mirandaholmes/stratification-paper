#SYSTEM = WormLoop
SYSTEM = Polymer
#SYSTEM = Trimer
PROGRAM = testBD$(SYSTEM)

CC=g++ -std=gnu++11 
OPTFLAGS = -O3 -I /home/holmes/Programs/
EXECUTABLE=./$(PROGRAM) 

all: $(EXECUTABLE) run 


$(PROGRAM): $(PROGRAM).o Particles.o BDSampler.o System.o
	$(CC) -o $(PROGRAM) $^

$(PROGRAM).o: $(PROGRAM).cpp Particles.hpp BDSampler.hpp  makefile
	$(CC) -c $(OPTFLAGS) $(PROGRAM).cpp 

BDSampler.o: BDSampler.cpp BDSampler.hpp Particles.hpp 
	$(CC) -c $(OPTFLAGS) BDSampler.cpp

Particles.o: Particles.hpp Particles.cpp makefile
	$(CC) -c $(OPTFLAGS) Particles.cpp 

System.o: Particles$(SYSTEM).cpp Particles.hpp makefile
	$(CC) -c $(OPTFLAGS) Particles$(SYSTEM).cpp -o System.o

run: $(EXECUTABLE)
	$(EXECUTABLE) $(ARGS) 

clean:
	rm -f *.o 


