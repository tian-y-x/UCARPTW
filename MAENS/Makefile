GPP=g++
CXX_STD = -std=gnu++11
CPP = $(GPP) $(CXX_STD)

program: main.cpp functions.h functions.cpp heutistic.h heuristic.cpp arrayoperations.cpp ins30.cpp ins30.h readmap.cpp MAENS/MAENS.cpp MAENS/MAENS.h MAENS/operators.cpp
	@mkdir -p bin/
	$(CPP) main.cpp functions.cpp heuristic.cpp arrayoperations.cpp ins30.cpp readmap.cpp MAENS/MAENS.cpp MAENS/operators.cpp -o bin/program
	@chmod +x bin/program
clean:
	@rm -rf bin/
.PHONY: clean