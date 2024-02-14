#! /bin/bash
g++ -c -std=c++0x \
  -Wextra -Wno-unused-function -Wno-unused-value -Wno-parentheses \
  -DLINUX *.cpp  \
  -static-libgcc -m32 -msse -pthread \
  -Wl,-Bstatic -lstdc++ -Wl,-Bstatic -lz -Wl,-Bstatic -lm -Wl,-Bstatic -ldl

ar rcs libnavmesh.a *.o