//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file Geant4-models/ATLAS-simplified/include/Constants.hh
/// \brief Definition of detector constants used in ATLAS-simplified project.

#ifndef CaloRConstants_h
#define CaloRConstants_h 1

#include "G4SystemOfUnits.hh"

#include "globals.hh" // used for G4int, G4double

enum class CaloIdx { WORLD, ECAL1, ECAL2, ECAL3, HCAL1, HCAL2, HCAL3 };

// Cell energy threshold
constexpr G4double CELL_THRESHOLD = 1e-5;

// Geometry parameters cell size of dR=0.1 in mm
constexpr G4double dR_01 = 144 * mm;

constexpr G4int fnLayers = 6;

//Material constants:
constexpr G4double denPb  = 11.35*g/cm3, aPb = 207.19*g/mole, X0Pb = 5.612 * mm, LintPb = 18.247 * cm;
constexpr G4double denLAr = 1.39*g/cm3, aLAr = 39.95*g/mole, X0LAr = 14.065 * cm, LintLAr = 86.078 * cm;
constexpr G4double denFe  = 7.87*g/cm3, aFe = 55.85*g/mole, X0Fe = 1.759 * cm, LintFe = 16.999 * cm;
constexpr G4double denSci  = 1.032*g/cm3, X0Sci = 42.544 * cm, LintSci = 69.969 * cm;


//ECAL (use Pb-LAr mix with 1.5:4.5 ratio)
constexpr G4double  dECAL_abs = 1.5*mm, dECAL_gap = 4.5*mm;
constexpr G4double  X0abs=X0Pb, X0gap=X0LAr, LintECALabs=LintPb, LintECALgap=LintLAr;
// In the following line we control ECAL X0, the ratio between Pb/LAr set to 1.1:4.5
constexpr G4double  X0ECAL=(1.1+4.5)/(1.1/X0abs+4.5/X0gap);// in mm
constexpr G4double  LintECAL=(1.1+4.5)/(1.1/LintECALabs+4.5/LintECALgap);// in mm
constexpr G4double  ECAL_nLayerPerX0 = 1;
constexpr G4double  ECAL1_X0 = 6, ECAL2_X0 = 16, ECAL3_X0 = 3;
constexpr G4double  ECAL1_dx = dR_01/32, ECAL1_dy = dR_01, ECAL1_dz = ECAL1_X0*X0ECAL;
constexpr G4double  ECAL2_dx = dR_01/4, ECAL2_dy = dR_01/4, ECAL2_dz = ECAL2_X0*X0ECAL;
constexpr G4double  ECAL3_dx = dR_01/2, ECAL3_dy = dR_01/4, ECAL3_dz = ECAL3_X0*X0ECAL;

//HCAL (use boxes of Fe+Scint with 5:1 ratio)
constexpr G4double  dHCAL_abs = 15.0*mm, dHCAL_gap = 3.0*mm, Lintabs=LintFe, Lintgap=LintSci;
constexpr G4double  HCAL1_Lint = 1.5, HCAL2_Lint = 4.1, HCAL3_Lint = 1.8;
constexpr G4double  HCAL_nLayerPerLint = 2; // 1.36
constexpr G4double  LintHCAL=1./(0.5*mm/(Lintabs)+0.5*mm/(Lintgap)); // in mm
//constexpr G4double  LintHCAL=1./(dHCAL_abs/(dHCAL_gap+dHCAL_abs)*mm/(Lintabs)+dHCAL_gap/(dHCAL_gap+dHCAL_abs)*mm/(Lintgap)); // in mm
constexpr G4double  HCAL1_dx = dR_01, HCAL1_dy = dR_01, HCAL1_dz = HCAL1_Lint*LintHCAL;
constexpr G4double  HCAL2_dx = dR_01, HCAL2_dy = dR_01, HCAL2_dz = HCAL2_Lint*LintHCAL;
constexpr G4double  HCAL3_dx = dR_01*2, HCAL3_dy = dR_01*2, HCAL3_dz = HCAL3_Lint*LintHCAL;

//Calibration paramters:
constexpr G4double ECAL_cell_calib = 3.69048;
constexpr G4double HCAL_cell_calib = 41.9115;

// Special constanst related to geometry settings:
constexpr G4int nSubCellDivECAL = 2*dECAL_gap/dECAL_abs;// (since dECAL_abs<dECAL_gap);
constexpr G4int nSubCellDivHCAL = 2;
constexpr G4int nSubCellDivECAL2 = (dECAL_gap+dECAL_abs)/dECAL_abs;
constexpr G4int nSubCellDivHCAL2 = (dHCAL_gap+dHCAL_abs)/dHCAL_gap;

// Subcell dZ width:
constexpr G4double dZHCAL1_subcell = HCAL1_Lint*LintHCAL/(G4int(HCAL1_Lint)*HCAL_nLayerPerLint);
constexpr G4double dZHCAL2_subcell = HCAL2_Lint*LintHCAL/(G4int(HCAL2_Lint)*HCAL_nLayerPerLint);
constexpr G4double dZHCAL3_subcell = HCAL3_Lint*LintHCAL/(G4int(HCAL3_Lint)*HCAL_nLayerPerLint);
  

// Detector geometry
constexpr G4double ECALThickness = ECAL1_dz + ECAL2_dz + ECAL3_dz;
constexpr G4double HCALThickness = HCAL1_dz + HCAL2_dz + HCAL3_dz;

constexpr G4double calorSizeXY = dR_01*16; // to cover dR  = 0.8 should be 16
constexpr G4double calorThickness = ECALThickness + HCALThickness;

// Calorimeter cell geometry
constexpr G4int kNofEm1Columns = (calorSizeXY/ECAL1_dx);
constexpr G4int kNofEm1Rows    = (calorSizeXY/ECAL1_dy);
constexpr G4int kNofEm1Cells   = kNofEm1Columns * kNofEm1Rows;
constexpr G4int kNofEm1ZRows   = ECAL_nLayerPerX0 * ECAL1_X0;
constexpr G4int kNofEm1SubCells= (ECAL1_dy / (dECAL_abs + dECAL_gap)) * kNofEm1ZRows;
constexpr G4int kNofEm2Columns = (calorSizeXY/ECAL2_dx);
constexpr G4int kNofEm2Rows    = (calorSizeXY/ECAL2_dy);
constexpr G4int kNofEm2Cells   = kNofEm2Columns * kNofEm2Rows;
constexpr G4int kNofEm2ZRows   = ECAL_nLayerPerX0 * ECAL2_X0;
constexpr G4int kNofEm2SubCells= (ECAL2_dy / (dECAL_abs + dECAL_gap)) * kNofEm2ZRows;
constexpr G4int kNofEm3Columns = (calorSizeXY/ECAL3_dx);
constexpr G4int kNofEm3Rows    = (calorSizeXY/ECAL3_dy);
constexpr G4int kNofEm3Cells   = kNofEm3Columns * kNofEm3Rows;
constexpr G4int kNofEm3ZRows   = ECAL_nLayerPerX0 * ECAL3_X0;
constexpr G4int kNofEm3SubCells= (ECAL3_dy / (dECAL_abs + dECAL_gap)) * kNofEm3ZRows;

constexpr G4int kNofHad1Columns = (calorSizeXY/HCAL1_dx);
constexpr G4int kNofHad1Rows    = (calorSizeXY/HCAL1_dy);
constexpr G4int kNofHad1Cells   = kNofHad1Columns * kNofHad1Rows;
constexpr G4int kNofHad1ZRows   = HCAL_nLayerPerLint * G4int(HCAL1_Lint);
constexpr G4int kNofHad1SubCells= (HCAL1_dy / (dHCAL_abs + dHCAL_gap)) * kNofHad1ZRows;
constexpr G4int kNofHad2Columns = (calorSizeXY/HCAL2_dx);
constexpr G4int kNofHad2Rows    = (calorSizeXY/HCAL2_dy);
constexpr G4int kNofHad2Cells   = kNofHad2Columns * kNofHad2Rows;
constexpr G4int kNofHad2ZRows   = HCAL_nLayerPerLint * G4int(HCAL2_Lint);
constexpr G4int kNofHad2SubCells= (HCAL2_dy / (dHCAL_abs + dHCAL_gap)) * kNofHad2ZRows;
constexpr G4int kNofHad3Columns = (calorSizeXY/HCAL3_dx);
constexpr G4int kNofHad3Rows    = (calorSizeXY/HCAL3_dy);
constexpr G4int kNofHad3Cells   = kNofHad3Columns * kNofHad3Rows;
constexpr G4int kNofHad3ZRows   = HCAL_nLayerPerLint * G4int(HCAL3_Lint);
constexpr G4int kNofHad3SubCells= (HCAL3_dy / (dHCAL_abs + dHCAL_gap)) * kNofHad3ZRows;

// Scintilator responce, return LY_electron / LY_proton for BC400 from 
// https://www.crystals.saint-gobain.com/sites/imdf.crystals.com/files/documents/sgc-organics-plastic-scintillators_0.pdf
constexpr G4double HCAL_HE_RATIO(G4double x){
if (x<0.1000) return 0.068401;
else if (x<0.1300) return 0.062265 + 0.061365*x;
else if (x<0.1700) return 0.054570 + 0.120551*x;
else if (x<0.2000) return 0.057777 + 0.101690*x;
else if (x<0.2370) return 0.072470 + 0.028222*x;
else if (x<0.3000) return 0.054462 + 0.104206*x;
else if (x<0.3330) return 0.054052 + 0.105574*x;
else if (x<0.4000) return 0.065021 + 0.072634*x;
else if (x<0.4700) return 0.025832 + 0.170605*x;
else if (x<0.6000) return 0.057359 + 0.103526*x;
else if (x<0.7100) return 0.046434 + 0.121735*x;
else if (x<0.8200) return 0.061743 + 0.100174*x;
else if (x<1.0000) return 0.040794 + 0.125721*x;
else if (x<1.2500) return 0.071927 + 0.094588*x;
else if (x<1.6677) return 0.117920 + 0.057793*x;
else if (x<2.0000) return 0.125106 + 0.053485*x;
else if (x<2.3500) return 0.043749 + 0.094163*x;
else if (x<3.0000) return 0.143393 + 0.051761*x;
else if (x<3.3600) return 0.231689 + 0.022329*x;
else if (x<4.0000) return 0.149860 + 0.046683*x;
else if (x<4.7000) return 0.204554 + 0.033010*x;
else if (x<6.0000) return 0.270425 + 0.018995*x;
else if (x<7.1000) return 0.327963 + 0.009405*x;
else if (x<8.3000) return 0.167233 + 0.032043*x;
else if (x<10.0000) return 0.347234 + 0.010356*x;
else if (x<12.8300) return 0.248581 + 0.020221*x;
else if (x<17.0000) return 0.444498 + 0.004951*x;
else if (x<20.0000) return 0.406916 + 0.007162*x;
else if (x<23.6400) return 0.385253 + 0.008245*x;
else if (x<30.0000) return 0.505827 + 0.003145*x;
else if (x<33.8000) return 0.443631 + 0.005218*x;
else if (x<40.0000) return 0.529026 + 0.002691*x;
else return 0.636681;
}

constexpr G4double HCAL_CE_RATIO(G4double x){
if (x<0.1000) return 0.010518;
else if (x<0.1300) return 0.011436 - 0.00918831*x;
else if (x<0.1700) return 0.012791 - 0.01960885*x;
else if (x<0.2000) return 0.010862 - 0.00826226*x;
else if (x<0.2370) return 0.010514 - 0.00652359*x;
else if (x<0.3000) return 0.010286 - 0.00555951*x;
else if (x<0.3330) return 0.008618;
else if (x<0.4000) return 0.011370 - 0.00826292*x;
else if (x<0.4700) return 0.008672 - 0.00151972*x;
else if (x<0.6000) return 0.009806 - 0.00393244*x;
else if (x<0.7100) return 0.007983 - 0.00089303*x;
else if (x<0.8200) return 0.009201 - 0.00260903*x;
else if (x<1.0000) return 0.007486 - 0.00051751*x;
else if (x<1.2500) return 0.009103 - 0.00213466*x;
else if (x<1.6677) return 0.007187 - 0.00060164*x;
else if (x<2.0000) return 0.007395 - 0.00072673*x;
else if (x<2.3500) return 0.006832 - 0.00044495*x;
else if (x<3.0000) return 0.005786 - 0.00000000*x;
else if (x<3.3600) return 0.006422 - 0.00021202*x;
else if (x<4.0000) return 0.006496 - 0.00023383*x;
else if (x<6.0000) return 0.005560;
else if (x<7.1000) return 0.005960 - 0.00006668*x;
else if (x<8.3000) return 0.005915 - 0.00006032*x;
else if (x<10.0000) return 0.005061 + 0.00004258*x;
else if (x<12.8300) return 0.005487 - 0.00000000*x;
else if (x<17.0000) return 0.003842 + 0.00012817*x;
else if (x<20.0000) return 0.006021;
else if (x<23.6400) return 0.003278 + 0.00013719*x;
else if (x<30.0000) return 0.005536 + 0.00004167*x;
else if (x<33.8000) return 0.002343 + 0.00014810*x;
else if (x<40.0000) return 0.004598 + 0.00008137*x;
else return 0.007853;
}

#endif
