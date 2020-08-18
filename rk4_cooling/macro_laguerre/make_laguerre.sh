g++ -o2 macro_laguerre.cpp -lm -o macro_laguerre `root-config --cflags --glibs`

# Make 3 equal folders and run each block in one
# Block 1
./macro_laguerre 0 0
./macro_laguerre 0 1
./macro_laguerre 0 2

# Block 2
./macro_laguerre 1 0
./macro_laguerre 1 1
./macro_laguerre 1 2

# Block 3
./macro_laguerre 2 0
./macro_laguerre 2 1
./macro_laguerre 2 2
