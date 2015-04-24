proter.out : main1.o
	g++ main1.o -o proter.out
	rm *.o

main1.o: main1.cpp
	g++ -c main1.cpp -o main1.o
	
clean:
	&nbps rm *.o proter.out
