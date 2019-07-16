
es10: es10.cpp random.o classi.o
	mpicxx -Wall --std=c++11 -o $@  $^ 
%.o: %.C %.h 
	mpicxx -Wall --std=c++11 -c -o $@  $<
clean: 
	rm *.o
	rm *~
