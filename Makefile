all: pointSetCrossMin
	make clean;

pointSetCrossMin: main.o JSONParser.o 
	g++ -o pointSetCrossMin main.o JSONParser.o

main.o:
	g++ -c -std=c++20 main.cpp

JSONParser.o:
	g++ -c -std=c++20 JSONParser.cpp 
		
clean:
	rm *.o;
