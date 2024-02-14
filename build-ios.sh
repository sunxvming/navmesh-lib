#! /bin/bash
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang  -c -std=c++0x \
  -Wextra -Wno-unused-function -Wno-unused-value -Wno-parentheses \
   *.cpp  \
  -arch arm64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS11.1.sdk \
  -Wl,-Bstatic -lstdc++ -Wl,-Bstatic -lz -Wl,-Bstatic -lm -Wl,-Bstatic -ldl

ar rcs libnavmesh.a *.o
