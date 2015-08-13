CXX = clang++ -std=c++11 -Xclang -fcolor-diagnostics -O2
LIBS = -lglpk -lboost_program_options

all: mwis

mwis: mwis.cpp
	$(CXX) -o mwis mwis.cpp $(LIBS)
