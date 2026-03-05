mkdir build
cd build
cmake .. -DENABLE_AVX512=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build .
cd ..
rm -rf align_benchmark
