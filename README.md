## Introduction

The Genetic Algorithm based calibration module introduced to EPANET 3 enables the calibration of hydraulic models of water distribution networks. In this context, 4 new files named GeneticAlgorithm.cpp, GeneticAlgorithm.h, Chromosome.cpp and Chromosome.h were generated in src. Through these files, the Observed.dat file in the Network folder is read and the content of the Observed.dat file is compared with the content of the Simulated file produced by the model. The goal of the program is to converge the model outputs given in the Simulated file to those in the Observed file. Thus, the calibration process can be performed. 

The Genetic Algorithm parameters such as the number of polulations, number of generations, mutation rate and crossover rate are in the main.cpp file and can be modified from there.

I hope this module will be useful to all researchers and engineers interested in the subject.

## Building

The source code can be compiled as both a shared library and a command-line executable. Any C++ compiler that supports the C++11 language standard can be used.

To build using CMake on Linux/Mac:

mkdir build
cd build
cmake .. 
make

The shared library (libepanet3.so) will be found in the /lib sub-directory and the command-line executable (run-epanet3) will be in the /bin sub-directory.

To build, use CMake on Windows with Visual Studio:
```
mkdir build
cd build
cmake -G "Visual Studio n yyyy" ..
cmake --build . --config Release
```
## Contributors
```
Mehmet Melih Koşucu,		Istanbul Technical University
