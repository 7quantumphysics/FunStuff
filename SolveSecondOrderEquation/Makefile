Executable: RungeKutta.o TemplateCode2.o
	g++ --std=c++11 RungeKutta.o TemplateCode2.o -o TemplateCode2.x

RungeKutta.o: RungeKutta.cc RungeKutta.h
	g++ --std=c++11 -c RungeKutta.cc

ChaosTest.o: ChaosTest.cpp
	g++ --std=c++11 -c TemplateCode2.cpp

clean_obj:
	rm *.o

clean_exe:
	rm *.x
clean_data:
	rm *.dat
