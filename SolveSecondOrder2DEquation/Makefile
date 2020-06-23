Executable: RungeKutta.o TemplateCode.o
	g++ --std=c++11 RungeKutta.o TemplateCode.o -o TemplateCode.x

RungeKutta.o: RungeKutta.cc RungeKutta.h
	g++ --std=c++11 -c RungeKutta.cc

ChaosTest.o: ChaosTest.cpp
	g++ --std=c++11 -c TemplateCode.cpp

clean_obj:
	rm *.o

clean_exe:
	rm *.x
