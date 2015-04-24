OBJS = main2.o
CC = gcc
DEBUG = -g

proter.out : $(OBJS)
	g++  -g -O -lGL -lglut $(OBJS) -o proter.out
	rm *.o

main2.o: main2.cpp prot_gl.h
	g++ -g -c main2.cpp -o main2.o
	
clean:
	rm *.o proter.out
