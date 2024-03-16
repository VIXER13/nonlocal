#!/bin/bash
export FOLDER=/mnt/c/repos/nonlocal
cd ${FOLDER}

rm -r build
mkdir build
cd build

conan config install ${FOLDER}/conan_settings/conan_default_settings/settings.yml 
conan profile update settings.compiler.libcxx=libstdc++11 default
conan profile update settings.build_type=Release default
conan install ./..

export FOLDER_BIN=/usr/bin
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j8



