APP_PLATFORM = android-9  
APP_ABI := armeabi armeabi-v7a  
APP_STL := gnustl_static  
APP_OPTIM := release
APP_BUILD_SCRIPT := Android.mk
APP_CPPFLAGS +=-std=c++11 #允许使用c++11的函数等功能