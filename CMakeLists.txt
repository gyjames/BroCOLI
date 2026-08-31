cmake_minimum_required(VERSION 3.16)
project(brocoli CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

option(BROCOLI_NATIVE "Build with -march=native (not portable across nodes)" OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib>=1.10)
find_package(Threads REQUIRED)

add_executable(brocoli src/main.cpp)
target_include_directories(brocoli PRIVATE src)
target_link_libraries(brocoli PRIVATE PkgConfig::HTSLIB Threads::Threads)
target_compile_options(brocoli PRIVATE -Wall -Wextra -Wno-unused-parameter)
if(BROCOLI_NATIVE)
  target_compile_options(brocoli PRIVATE -march=native)
endif()

enable_testing()
add_executable(brocoli_test tests/test_pipeline.cpp)
target_include_directories(brocoli_test PRIVATE src)
add_test(NAME pipeline COMMAND brocoli_test)

# preBroCOLI: built only when edlib is vendored and WFA2-lib is available.
set(EDLIB_DIR ${CMAKE_SOURCE_DIR}/third_party/edlib)
if(EXISTS ${EDLIB_DIR}/edlib.cpp AND DEFINED ENV{CONDA_PREFIX}
   AND EXISTS $ENV{CONDA_PREFIX}/include/wfa2lib/bindings/cpp/WFAligner.hpp)
  find_package(ZLIB REQUIRED)
  add_executable(prebrocoli pre/main.cpp ${EDLIB_DIR}/edlib.cpp)
  target_include_directories(prebrocoli PRIVATE src pre ${EDLIB_DIR}
                             $ENV{CONDA_PREFIX}/include/wfa2lib
                             $ENV{CONDA_PREFIX}/include)
  target_link_directories(prebrocoli PRIVATE $ENV{CONDA_PREFIX}/lib)
  target_link_libraries(prebrocoli PRIVATE wfa2cpp wfa2 ZLIB::ZLIB Threads::Threads)
  target_compile_options(prebrocoli PRIVATE -Wall -Wextra -Wno-unused-parameter)
  add_executable(prebrocoli_seqio_test tests/test_seqio.cpp)
  target_include_directories(prebrocoli_seqio_test PRIVATE pre)
  target_link_libraries(prebrocoli_seqio_test PRIVATE ZLIB::ZLIB)
  add_test(NAME seqio COMMAND prebrocoli_seqio_test)
  install(TARGETS prebrocoli RUNTIME DESTINATION bin)
else()
  message(STATUS "preBroCOLI skipped (need third_party/edlib and conda wfa2-lib)")
endif()

install(TARGETS brocoli RUNTIME DESTINATION bin)
