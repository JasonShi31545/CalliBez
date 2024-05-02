##
# CalliBez
#
# @file
# @version 0.1

CC=g++
CPPFLAGS=-std=c++17 -Wall -Wpedantic -g -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -O0 -Wno-unused-result -fstack-protector
LIBS=-lSDL2 -lm

PROD_CPPFLAGS=-std=c++17 -Wall -Wpedantic -O2 -Wno-unused-result -fstack-protector -s
PROD_LIBS=-lSDL2 -lm

# %: %.cc
# 	$(CC) $(CPPFLAGS) $(LIBS) $< -o $@.out

.DEFAULT_GOAL := debug

debug:
	$(CC) *.cc $(CPPFLAGS) -c
	$(CC) *.o $(LIBS) $(CPPFLAGS) -o db.out

prod:
	$(CC) *.cc $(PROD_CPPFLAGS) -c
	$(CC) *.o $(PROD_LIBS) $(PROD_CPPFLAGS) -o prod.out

clean:
	rm -f *.o
	rm -f *.out


# end
