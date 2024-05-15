#!/bin/bash

g++ -std=c++17 -Wall -O setup.cc bezier.cc new.cc -I/usr/include/eigen3 -lSDL2 -lm -o new.out && ./new.out
