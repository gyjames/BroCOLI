#CMake版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(BroCOLI LANGUAGES C CXX)

# 设置 C++ 标准为 C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# 设置线程支持
find_package(Threads REQUIRED)

# 添加 Eigen 头文件路径
# include_directories(${CMAKE_SOURCE_DIR}/src)

# 创建可执行文件
add_executable(BroCOLI_sc src/BroCOLI_sc.cpp src/edlib.cpp)
target_include_directories(BroCOLI_sc PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(BroCOLI_sc Threads::Threads)

# 创建可执行文件
add_executable(BroCOLI_bulk src/BroCOLI_bulk.cpp)
target_include_directories(BroCOLI_bulk PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(BroCOLI_bulk Threads::Threads)


set(CMAKE_C_STANDARD 99)
add_subdirectory(WFA2-lib)
if(TARGET wfa2cpp_static)
  target_compile_options(wfa2cpp_static PRIVATE -std=gnu++11)
  set_property(TARGET wfa2cpp_static PROPERTY CXX_STANDARD 11)
  set_property(TARGET wfa2cpp_static PROPERTY CXX_STANDARD_REQUIRED ON)
endif()

target_compile_options(wfa2cpp PRIVATE -std=c++11)
set_property(TARGET wfa2cpp PROPERTY CXX_STANDARD 11)
set_property(TARGET wfa2cpp PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(preBroCOLI src/preBroCOLI.cpp src/edlib.cpp)
target_include_directories(preBroCOLI PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(preBroCOLI
  Threads::Threads
  wfa2cpp
  wfa2
)











