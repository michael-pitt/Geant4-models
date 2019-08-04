# Geant4-models
This repository contains few Geant4 models:

* [ATLAS-simplified](ATLAS-simplified): A simplified version of ATLAS detector calorimeter.

## Download and install
To clone the examples from the repository just run:
```bash
git clone https://USERNAME:PASSWORD@github.com/mpitt82/Geant4-models.git
```
The examples are build with CMake cross-platform, and compiled with GEANT4 instalation. 
To build a specific model, one need to setup corresponding ENV:
```bash
mkdir -pv Geant4-models/EXAMPLE/build; cd Geant4-models/EXAMPLE/build
cmake ../
make -j
#lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"
#lsetup "sft releases/LCG_95/Geant4/10.05"
#lsetup "sft releases/LCG_95/clhep/2.4.1.0"
```

now one can build the model:
```bash
mkdir -pv Geant4-models/EXAMPLE/build; cd Geant4-models/EXAMPLE/build
cmake ../
```
where 

## ATLAS-simplified

A simplified version of the ATLAS detector calorimeter. This setup incorporates the concept of a sampling calorimeter.
The detector granularity matches the ATLAS geometry at &eta;=0 (barrel). The detector composited from 6 layers:
3 form the Electro-Magnetic Calorimeter (ECAL) and 3 forms the Hadronic Calorimeter (HCAL).

