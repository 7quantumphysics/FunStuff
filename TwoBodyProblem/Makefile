Executable: RungeKutta.o EarthSunOrbit.o
	g++ --std=c++11 RungeKutta.o EarthSunOrbit.o -o EarthSunOrbit.x

RungeKutta.o: RungeKutta.cc RungeKutta.h
	g++ --std=c++11 -c RungeKutta.cc

ChaosTest.o: ChaosTest.cpp
	g++ --std=c++11 -c EarthSunOrbit.cpp

clean_obj:
	rm *.o

clean_exe:
	rm *.x

clean_data:
	rm *.dat	
