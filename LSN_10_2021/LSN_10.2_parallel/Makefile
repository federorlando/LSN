CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o functions.o TSP_position.o
	$(CC) random.o main.o functions.o TSP_position.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
TSP_position.o : TSP_position.cpp TSP_position.h
	$(CC) -c TSP_position.cpp -o TSP_position.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
