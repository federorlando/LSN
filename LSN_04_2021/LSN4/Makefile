CC = g++
CFLAGS = -Wall -O3 --std=c++11

#main.exe : main.o random.o functions.o MolDyn_NVE.o
	#$(CC) random.o main.o functions.o MolDyn_NVE.o -o main.exe
#main.o : main.cpp
	#$(CC) -c main.cpp -o main.o $(CFLAGS)
MolDyn_NVE.exe : MolDyn_NVE.o random.o
	$(CC) random.o MolDyn_NVE.o -o MolDyn_NVE.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
