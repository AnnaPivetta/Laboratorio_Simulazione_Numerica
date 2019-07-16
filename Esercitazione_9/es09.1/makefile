
es091: es091.cpp random.o classi.o
	g++ -Wall --std=c++11 -o $@  $^ 
%.o: %.C %.h 
	g++ -Wall --std=c++11 -c -o $@  $<
clean: 
	rm *.o
	rm *~
