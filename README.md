# Geant4-models
This repository contains few Geant4 models:

* [ATLAS-simplified](ATLAS-simplified): A simplified version of ATLAS detector calorimeter.

## Download and install
To clone the examples from a git just run:
```bash
git clone https://USERNAME:PASSWORD@github.com/mpitt82/Geant4-models.git
```
The examples are build with CMake cross-platform, to build a specific model, run:
```bash
mkdir -pv Geant4-models/EXAMPLE/build; cd Geant4-models/EXAMPLE/build
cmake ../
```

## ATLAS-simplified

A simplified version of the ATLAS detector calorimeter. This setup incorporates the concept of a sampling calorimeter.
The detector granularity matches the ATLAS geometry at &eta;=0 (barrel). The detector composited from 6 layers:
3 form the Electro-Magnetic Calorimeter (ECAL) and 3 forms the Hadronic Calorimeter (HCAL).

