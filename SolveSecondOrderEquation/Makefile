Executable: RungeKutta.o ChaosTest.o
	g++ --std=c++11 RungeKutta.o ChaosTest.o -o ChaosTest.x

RungeKutta.o: RungeKutta.cc RungeKutta.h
	g++ --std=c++11 -c RungeKutta.cc

ChaosTest.o: ChaosTest.cpp
	g++ --std=c++11 -c ChaosTest.cpp

clean_obj:
	rm *.o

clean_exe:
	rm ChaosTest.x
