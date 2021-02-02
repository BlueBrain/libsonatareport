![libSonata Logo](logo/libSonataLogo.jpg)
C++ / Python reader for SONATA reports files:
https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md

[![Coverage Status](https://coveralls.io/repos/github/BlueBrain/libsonatareport/badge.svg)](https://coveralls.io/github/BlueBrain/libsonatareport)

# Installation

## Building the C++ library

```shell
git clone git@github.com:BlueBrain/libsonatareport.git --recursive
cd libsonatareport
mkdir build && cd build
cmake  -DCMAKE_BUILD_TYPE=Release  -DEXTLIB_FROM_SUBMODULES=ON ..
make -j
```
