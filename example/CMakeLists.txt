cmake_minimum_required(VERSION 3.0)

project(navmesh-test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")

# 定义工程根目录; CMAKE_SOURCE_DIR为内建变量，表示工程根目录的CMakeLists.txt文件路径
SET(ROOT_DIR ${CMAKE_SOURCE_DIR})

file(GLOB SOURCES "*.cpp")


add_executable(${PROJECT_NAME} ${SOURCES})

set_target_properties(${PROJECT_NAME}  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)