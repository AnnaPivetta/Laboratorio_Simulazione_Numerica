
es031: es031.cpp random.o
	g++ -Wall -o $@ $^ 
%.o: %.C %.h 
	g++ -Wall -c -o $@ $<
clean: 
	rm *.o
	rm *~
