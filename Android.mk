#编译动态库    
LOCAL_PATH := $(call my-dir)    
include $(CLEAR_VARS)    
LOCAL_MODULE = NavPath    
LOCAL_CFLAGS = $(L_CFLAGS)

#声明一个变量MY_CPP_PATH表示源码目录
MY_CPP_PATH := $(LOCAL_PATH)

My_All_Files := $(shell find $(MY_CPP_PATH)/.)
My_All_Files := $(My_All_Files:$(MY_CPP_PATH)/./%=$(MY_CPP_PATH)%)

#从My_All_Files中再次提取所有的cpp文件,这里也可以使用filter函数
MY_CPP_LIST := $(foreach c_file,$(My_All_Files), $(wildcard $(c_file)/*.cpp) )
MY_CPP_LIST := $(MY_CPP_LIST:$(LOCAL_PATH)/%=%)

$(warning $(MY_CPP_LIST))
LOCAL_SRC_FILES := $(MY_CPP_LIST)
APP_STL := gnustl_static
LOCAL_CPPFLAGS :=--std=c++11
LOCAL_LDLIBS    := -lm -llog
include $(BUILD_SHARED_LIBRARY) 