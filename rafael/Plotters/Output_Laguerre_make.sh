
g++ -o2 Output_Laguerre_v$1.cpp -lm -o Output_Laguerre_v$1 `root-config --cflags --glibs`
./Output_Laguerre_v$1
rm Output_Laguerre_v$1
xdg-open Output_Laguerre_v$1_data/$2

# v1 - Laguerre Gaussian without Bx, px=0
# v2 - Simplified Laguerre Gaussian
# v3 - Simplified Laguerre Gaussian, kdamped
# v4 - Simplified Laguerre Gaussian with Ex
# v5 - Simplified Laguerre Gaussian with Ex and A0 = 30
