CC=g++

CFLAGS = -c -O3 -I/usr/include/eigen3 -fopenmp

all: hello

hello: TITAN.o Cell.o Geo_param.o Phys_param.o Gran.o Help_funk.o Luch.o Setka.o Setka_physics.o Setka_move.o Surfaces.o Yzel.o
	$(CC) TITAN.o Cell.o Geo_param.o Phys_param.o Gran.o Help_funk.o Luch.o Setka.o Setka_physics.o Setka_move.o Surfaces.o Yzel.o -o hello  -fopenmp

TITAN.o: TITAN.cpp
	$(CC) $(CFLAGS) TITAN.cpp

Cell.o: Cell.cpp
	$(CC) $(CFLAGS) Cell.cpp

Geo_param.o: Geo_param.cpp
	$(CC) $(CFLAGS) Geo_param.cpp

Gran.o: Gran.cpp
	$(CC) $(CFLAGS) Gran.cpp

Help_funk.o: Help_funk.cpp
	$(CC) $(CFLAGS) Help_funk.cpp

Luch.o: Luch.cpp
	$(CC) $(CFLAGS) Luch.cpp

Setka.o: Setka.cpp
	$(CC) $(CFLAGS) Setka.cpp

Setka_physics.o: Setka_physics.cpp
	$(CC) $(CFLAGS) Setka_physics.cpp

Setka_move.o: Setka_move.cpp
	$(CC) $(CFLAGS) Setka_move.cpp


Surfaces.o: Surfaces.cpp
	$(CC) $(CFLAGS) Surfaces.cpp

Yzel.o: Yzel.cpp
	$(CC) $(CFLAGS) Yzel.cpp

Phys_param.o: Phys_param.cpp
	$(CC) $(CFLAGS) Phys_param.cpp


clean:
	rm -rf *.o hello