# Simple Makefile with variables (compiler+flags+include)
# and explicit paths on targets + dependencies
#
#
# Author: 

NAME := main

CC := g++-7
CCFLAGS := -std=c++14
ROOT := `root-config --cflags --glibs`
GTK1 := `pkg-config --cflags gtk+-3.0`
GTK2 := `pkg-config --libs gtk+-3.0`

INCLUDE := src

SOURCES1 := $(shell ls src | grep cpp & ls memes | grep cpp)
SOURCES2 := $(shell ls src | grep cc & ls memes | grep cc)
OBJECTS1 := $(SOURCES1:.cpp=.o)
OBJECTS2 := $(SOURCES2:.cc=.o)
OBJECTS := $(addprefix obj/, $(OBJECTS1)) $(addprefix obj/, $(OBJECTS2))

#####################################################

all: bin/$(NAME)

bin/$(NAME): $(OBJECTS)
	@echo " yo compiling $^  [$@]"
	@$(CC) $(GTK1) -o $@ $^ -I $(INCLUDE) $(GTK2) $(ROOT)

obj/%.o: memes/%.cpp
	@echo " ye compiling $^  [$@]"
	@$(CC) $(GTK1) $(CCFLAGS) -c $^ -I $(INCLUDE) $(ROOT) $(GTK2) -o $@

obj/%.o: src/%.cpp
	@echo " ye compiling $^  [$@]"
	@$(CC) $(GTK1) $(CCFLAGS) -c $^ $(ROOT) $(GTK2) -o $@

clean:
	@find obj/ -type f -name '*.o' -delete
	@find bin/ -type f -name '*' -delete

run:
	./bin/$(NAME)

meme:
	make clean
	make 
	make run
