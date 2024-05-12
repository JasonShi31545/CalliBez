#!/bin/bash

g++ -std=c++17 -Wall -O setup.cc bezier.cc new.cc -lSDL2 -lm -o new
./new
