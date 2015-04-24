proter.out : main2.o
	g++ main2.o -o proter.out
	rm *.o

main1.o: main2.cpp
	g++ -c main2.cpp -o main2.o
	
clean:
	&nbps rm *.o proter.out
