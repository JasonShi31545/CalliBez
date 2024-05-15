##
# CalliBez
#
# @file
# @version 0.1

CC=g++
CPPFLAGS=-std=c++17 -Wall -Wpedantic -g -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -O0 -Wno-unused-result -fstack-protector
INCLUDES=-I/usr/include/eigen3
LIBS=-lSDL2 -lm

PROD_CPPFLAGS=-std=c++17 -Wall -Wpedantic -O2 -Wno-unused-result -fstack-protector -s

# %: %.cc
# 	$(CC) $(CPPFLAGS) $(LIBS) $< -o $@.out

.DEFAULT_GOAL := debug

debug:
	$(CC) *.cc $(INCLUDES) $(CPPFLAGS) -c
	$(CC) *.o $(LIBS) $(CPPFLAGS) -o db.out

prod:
	$(CC) *.cc $(INCLUDES) $(PROD_CPPFLAGS) -c
	$(CC) *.o $(LIBS) $(PROD_CPPFLAGS) -o prod.out

clean:
	rm -f *.o
	rm -f *.out


# end
