CC=g++
CFLAGS=-c -Wall

all: MOPSI

MOPSI: main.o pi.o simplex.o matrice.o
	$(CC) -lX11 -lpthread main.o pi.o simplex.o matrice.o -o MOPSI

main.o: ../source/main.cpp
	$(CC) $(CFLAGS) ../source/main.cpp

pi.o: ../source/pi.cpp
	$(CC) $(CFLAGS) ../source/pi.cpp

simplex.o: ../source/simplex.cpp
	$(CC) $(CFLAGS) ../source/simplex.cpp

matrice.o: ../source/matrice.cpp
	$(CC) $(CFLAGS) ../source/matrice.cpp
