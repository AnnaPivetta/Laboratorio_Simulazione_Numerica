
es051: es051.cpp random.o
	g++ -Wall --std=c++11 -o $@ $^ 
%.o: %.cpp %.h 
	g++ -Wall --std=c++11 -c -o $@ $< 
clean: 
	rm *.o
	rm *~
