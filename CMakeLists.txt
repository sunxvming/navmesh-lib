cmake_minimum_required(VERSION 3.0)

project(navmesh)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")

# 定义工程根目录; CMAKE_SOURCE_DIR为内建变量，表示工程根目录的CMakeLists.txt文件路径
SET(ROOT_DIR ${CMAKE_SOURCE_DIR})

# 选项：控制库的类型
option(BUILD_SHARED_LIBS "Build shared library" ON)

# 编译那些源码
file(GLOB SOURCES "src/*.cpp")

# 添加库
if(BUILD_SHARED_LIBS)
    message("build share lib")
    add_library(${PROJECT_NAME} SHARED ${SOURCES})
else()
    message("build static lib")
    add_library(${PROJECT_NAME} STATIC ${SOURCES})
endif()


set_target_properties(${PROJECT_NAME}  PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
)


add_subdirectory(example)