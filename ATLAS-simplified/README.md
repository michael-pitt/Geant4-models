# Simplified ATLAS calorimeter geometry

A (very) simplified version of [ATLAS detector](https://iopscience.iop.org/article/10.1088/1748-0221/3/08/S08003) callorimeter incorporating the concept of a sampling callorimeter.

### Instalation and execution
```bash
export Geant4_DIR=PATH_TO_GEANT4_INSTALL
#if using CERN/CVMFS with CENTOS7
#source /cvmfs/sft.cern.ch/lcg/releases/LCG_95/Geant4/10.05/x86_64-centos7-gcc8-opt/Geant4-env.sh
mkdir build; cd build
cmake ../; make -j
```

## Geometry discription:

Full technical details on the electro-magnetic calorimeter (ECAL) design can be found in [ECAL ref](https://cds.cern.ch/record/331061/files/CERN-LHCC-96-41.pdf)
and  about hadronic calorimeter (HCAL) design can be found in [HCAL ref](https://cds.cern.ch/record/2004868/files/ATL-TILECAL-PROC-2015-002.pdf).

The detector granularity matches ATLAS geometry at &eta;=0 (barrel). The detector conposited from 6 layers:
3 form the Electro-Magnetic Callorimeter and 3 form the Hadronid Callorimeter.
A side view of the detector setup is shown in Fig. 1.

![Fig 1: Scheme of detector layers](images/calorimeter_layers.png)

### ECAL geometry:

Electro-Magnetic calorimeter composed from 3 layers. It consist of pasive layers,
1.5 mm each of lead (Pb) and a 4.5 mm gap of liquid argon (LAr) used as a sensitive layer. The passibe lead layers 
are placed in accordion geometry shape and the gaps were filled with LAr as shown in Fig 2.

![Fig 2: Accordion geometry of ECAL](images/ECAL_L1.png)

The energy resolution of ECAL geometry was simulated to be &sigma;/E=6.5%/&sqrt;E, which is controled by electro-magnetic energy
 sampling fraction of about 25%.
 
Cell granularity was defined similar to ATLAS ECAL granularity at &eta;=0:

.... (to be filled)



