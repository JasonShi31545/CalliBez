##
# CalliBez
#
# @file
# @version 0.1

CC=g++
CPPFLAGS=-std=c++17 -Wall -Wpedantic -g -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -O0 -Wno-unused-result
LIBS=-lSDL2 -lm

%: %.cc
	$(CC) $(CPPFLAGS) $(LIBS) $< -o $@.out

# end
