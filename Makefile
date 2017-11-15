.PHONY: all

libs = -larmadillo 
args = -O2 -g -std=c++11 
src = src/main.cpp

lsesolver: $(src) 
	g++ $(src) $(args) $(libs) -o lsesolver

all: lsesolver

