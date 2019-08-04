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
/// \file Geant4-models/ATLAS-simplified/src/CellParameterisation.cc
/// \brief Implementation of the CellParameterisation class

#include "CellParameterisation.hh"
#include "Constants.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

// Primary cells
Em1CellParameterisation::Em1CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm1Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::ECAL1);
    fXEm1Cell[copyNo] = XY[0];
    fYEm1Cell[copyNo] = XY[1];
  }
}
Em1CellParameterisation::~Em1CellParameterisation() {}
void Em1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXEm1Cell[copyNo],fYEm1Cell[copyNo],0.));}


Em2CellParameterisation::Em2CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm2Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::ECAL2);
    fXEm2Cell[copyNo] = XY[0];
    fYEm2Cell[copyNo] = XY[1];
  }
}
Em2CellParameterisation::~Em2CellParameterisation() {}
void Em2CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXEm2Cell[copyNo],fYEm2Cell[copyNo],0.));}

Em3CellParameterisation::Em3CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm3Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::ECAL3);
    fXEm3Cell[copyNo] = XY[0];
    fYEm3Cell[copyNo] = XY[1];
  }
}
Em3CellParameterisation::~Em3CellParameterisation() {}
void Em3CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXEm3Cell[copyNo],fYEm3Cell[copyNo],0.));}

Had1CellParameterisation::Had1CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad1Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::HCAL1);
    fXHad1Cell[copyNo] = XY[0];
    fYHad1Cell[copyNo] = XY[1];
  }
}
Had1CellParameterisation::~Had1CellParameterisation() {}
void Had1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXHad1Cell[copyNo],fYHad1Cell[copyNo],0.));}


Had2CellParameterisation::Had2CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad2Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::HCAL2);
    fXHad2Cell[copyNo] = XY[0];
    fYHad2Cell[copyNo] = XY[1];
  }
}
Had2CellParameterisation::~Had2CellParameterisation() {}
void Had2CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXHad2Cell[copyNo],fYHad2Cell[copyNo],0.));}

Had3CellParameterisation::Had3CellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad3Cells; copyNo++) {
	G4double *XY = CopyToXY(copyNo,CaloIdx::HCAL3);
    fXHad3Cell[copyNo] = XY[0];
    fYHad3Cell[copyNo] = XY[1];
  }
}
Had3CellParameterisation::~Had3CellParameterisation() {}
void Had3CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(fXHad3Cell[copyNo],fYHad3Cell[copyNo],0.));}


//sub-cells (for absorber/scintilator)
Em1SubCellParameterisation::Em1SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm1SubCells; copyNo++) {
    auto column = copyNo / kNofEm1ZRows ;
    auto row = copyNo % kNofEm1ZRows;
    fYEm1SubCell[copyNo] = (column+0.5)*(dECAL_abs+dECAL_gap) - ECAL1_dy/2;
    fZEm1SubCell[copyNo] = (row+0.5)*X0ECAL/ECAL_nLayerPerX0 - ECAL1_dz/2;
  }
}
Em1SubCellParameterisation::~Em1SubCellParameterisation() {}
void Em1SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYEm1SubCell[copyNo],fZEm1SubCell[copyNo]));}

Em2SubCellParameterisation::Em2SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm2SubCells; copyNo++) {
    auto column = copyNo / kNofEm2ZRows ;
    auto row = copyNo % kNofEm2ZRows;
    fYEm2SubCell[copyNo] = (column+0.5)*(dECAL_abs+dECAL_gap) - ECAL2_dy/2;
    fZEm2SubCell[copyNo] = (row+0.5)*X0ECAL/ECAL_nLayerPerX0 - ECAL2_dz/2;
  }
}
Em2SubCellParameterisation::~Em2SubCellParameterisation() {}
void Em2SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYEm2SubCell[copyNo],fZEm2SubCell[copyNo]));}

Em3SubCellParameterisation::Em3SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofEm3SubCells; copyNo++) {
    auto column = copyNo / kNofEm3ZRows ;
    auto row = copyNo % kNofEm3ZRows;
    fYEm3SubCell[copyNo] = (column+0.5)*(dECAL_abs+dECAL_gap) - ECAL3_dy/2;
    fZEm3SubCell[copyNo] = (row+0.5)*X0ECAL/ECAL_nLayerPerX0 - ECAL3_dz/2;
  }
}
Em3SubCellParameterisation::~Em3SubCellParameterisation() {}
void Em3SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYEm3SubCell[copyNo],fZEm3SubCell[copyNo]));}


Had1SubCellParameterisation::Had1SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad1SubCells; copyNo++) {
    auto column = copyNo / kNofHad1ZRows ;
    auto row = copyNo % kNofHad1ZRows;
    fYHad1SubCell[copyNo] = (column+0.5)*(dHCAL_abs+dHCAL_gap) - HCAL1_dy/2;
    fZHad1SubCell[copyNo] = (row+0.5)*dZHCAL1_subcell - HCAL1_dz/2;
  }
}
Had1SubCellParameterisation::~Had1SubCellParameterisation() {}
void Had1SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYHad1SubCell[copyNo],fZHad1SubCell[copyNo]));}

Had2SubCellParameterisation::Had2SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad2SubCells; copyNo++) {
    auto column = copyNo / kNofHad2ZRows ;
    auto row = copyNo % kNofHad2ZRows;
    fYHad2SubCell[copyNo] = (column+0.5)*(dHCAL_abs+dHCAL_gap) - HCAL2_dy/2;
    fZHad2SubCell[copyNo] = (row+0.5)*dZHCAL2_subcell - HCAL2_dz/2;
  }
}
Had2SubCellParameterisation::~Had2SubCellParameterisation() {}
void Had2SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYHad2SubCell[copyNo],fZHad2SubCell[copyNo]));}

Had3SubCellParameterisation::Had3SubCellParameterisation() : G4VPVParameterisation() {
  for (auto copyNo=0; copyNo<kNofHad3SubCells; copyNo++) {
    auto column = copyNo / kNofHad3ZRows ;
    auto row = copyNo % kNofHad3ZRows;
    fYHad3SubCell[copyNo] = (column+0.5)*(dHCAL_abs+dHCAL_gap) - HCAL3_dy/2;
    fZHad3SubCell[copyNo] = (row+0.5)*dZHCAL3_subcell - HCAL3_dz/2;
  }
}
Had3SubCellParameterisation::~Had3SubCellParameterisation() {}
void Had3SubCellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  physVol->SetTranslation(G4ThreeVector(0.,fYHad3SubCell[copyNo],fZHad3SubCell[copyNo]));}

