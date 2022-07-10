CC = g++
CFLAGS = -Wall -O3 --std=c++11

MD.exe : MD.o random.o
	$(CC) random.o MD.o -o MD.exe
MD.o : MD.cpp
	$(CC) -c MD.cpp -o MD.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o MD.exe
