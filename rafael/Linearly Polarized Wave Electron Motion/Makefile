NAME := Main.exe

CC := g++-7
CCFLAGS := -std=c++14

INCLUDE := src

SOURCES := $(shell ls src | grep cpp & ls | grep cpp)
OBJECTS := $(SOURCES:.cpp=.o)

#####################################################

all: bin/$(NAME)

bin/$(NAME): $(OBJECTS)
	$(CC) -o $@ $(addprefix bin/, $^) -I $(INCLUDE)  `root-config --cflags --glibs`

%.o: %.cpp
	$(CC) $(CCFLAGS) -o bin/$@ -c $^ -I $(INCLUDE) `root-config --cflags --glibs`

%.o: src/%.cpp	
	$(CC) $(CCFLAGS) -o bin/$@ -c $^ `root-config --cflags --glibs`

clean:
	rm bin/*

run:
	./bin/$(NAME)

val: 
	valgrind ./bin/$(NAME)

