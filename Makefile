all: executable

executable: meshgen.o meshestimator.o pathestimator.o
	g++ meshgen.o meshestimator.o pathestimator.o -o executable

meshgen.o: meshgen.cpp
	g++ -c meshgen.cpp

meshestimator.o: meshestimator.cpp
	g++ -c meshestimator.cpp

pathestimator.o: pathestimator.cpp
	g++ -c pathestimator.cpp

clean:
	rm -rf *o mesh
