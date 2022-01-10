![libSonata Logo](logo/libSonataLogo.jpg)
C++ / Python reader for SONATA reports files:
https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md

[![Coverage Status](https://coveralls.io/repos/github/BlueBrain/libsonatareport/badge.svg)](https://coveralls.io/github/BlueBrain/libsonatareport)
![clang-format](https://github.com/BlueBrain/libsonatareport/workflows/clang-format-check/badge.svg)
![unit tests](https://github.com/BlueBrain/libsonatareport/workflows/run-test/badge.svg)

# Installation

## Building the C++ library

```shell
git clone git@github.com:BlueBrain/libsonatareport.git --recursive
cd libsonatareport
mkdir build && cd build
cmake  -DCMAKE_BUILD_TYPE=Release  -DSONATA_REPORT_ENABLE_SUBMODULES=ON ..
make -j
```


# Funding & Acknowledgment
 
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright © 2019-2021 Blue Brain Project/EPFL
