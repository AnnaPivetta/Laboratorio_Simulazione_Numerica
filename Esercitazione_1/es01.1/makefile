
es011: es011.cpp random.o
	g++ -Wall -o $@ $^ 
%.o: %.C %.h 
	g++ -Wall -c -o $@ $<
clean: 
	rm *.o
	rm *~
