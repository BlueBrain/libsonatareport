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
cmake -B build -DCMAKE_BUILD_TYPE=Release -DSONATA_REPORT_ENABLE_SUBMODULES=ON -GNinja
cmake --build build
```

# Usage

This section provides guidance on how to use the tools provided by libsonatareport for converting report and spike files into the SONATA format.

## Prerequisites

Before using the conversion tools, ensure the library is compiled with the converter enabled by including the cmake flag `-DSONATA_REPORT_ENABLE_CONVERTER=ON` during the installation process described in the "Installation" section.

## Tools

### Report Converter

Converts `.bbp` format files, typically containing voltage or current reports from simulations, into the SONATA format.

```shell
reports_converter <report_filename> <report_type> [population_name]

Example:
    reports_converter soma.bbp --soma PopulationA
```
This command will convert the soma report `soma.bbp` for population 'PopulationA' into the SONATA format as `soma.bbp.h5`

### Spikes Converter

The spikes converter tool is designed to convert old `.dat` spike files into the SONATA format.

```shell
spikes_converter <spike_filename> [population_name]

Example:
    spikes_converter out.dat PopulationB
```

This command will convert the spike file `out.dat` for population 'PopulationB' into the SONATA format as `out.dat.h5`

## Additional Information

For more detailed information on the usage of these converters, including additional options and examples, you can run them with the `--help` flag to display the usage:

```shell
reports_converter --help
spikes_converter --help
```

# Funding & Acknowledgment
 
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright (c) 2021-2022 Blue Brain Project/EPFL
