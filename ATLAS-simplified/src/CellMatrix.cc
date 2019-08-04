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
/// \file Geant4-models/ATLAS-simplified/src/CellMatrix.cc
/// \brief Implementation of the CellMatrix class

#include "CellMatrix.hh"

// Primary cells
CaloRCellMatrix::CaloRCellMatrix(){}
CaloRCellMatrix::~CaloRCellMatrix() {}

// Converters:

G4double* CaloRCellMatrix::CopyToXY(G4int copyNo, CaloIdx layer){
  G4double X(0), Y(0);
  if(layer==CaloIdx::ECAL1){
	auto column = copyNo / kNofEm1Rows ;
	auto row = copyNo % kNofEm1Rows;
	X = (column+0.5)*ECAL1_dx - calorSizeXY/2;
	Y = (row+0.5)*ECAL1_dy - calorSizeXY/2;
  }
  else if(layer==CaloIdx::ECAL2){
	auto column = copyNo / kNofEm2Rows ;
	auto row = copyNo % kNofEm2Rows;
	X = (column+0.5)*ECAL2_dx - calorSizeXY/2;
	Y = (row+0.5)*ECAL2_dy - calorSizeXY/2;
  }   
  else if(layer==CaloIdx::ECAL3){
	auto column = copyNo / kNofEm3Rows ;
	auto row = copyNo % kNofEm3Rows;
	X = (column+0.5)*ECAL3_dx - calorSizeXY/2;
	Y = (row+0.5)*ECAL3_dy - calorSizeXY/2;
  }  
  else if(layer==CaloIdx::HCAL1){
	auto column = copyNo / kNofHad1Rows ;
	auto row = copyNo % kNofHad1Rows;
	X = (column+0.5)*HCAL1_dx - calorSizeXY/2;
	Y = (row+0.5)*HCAL1_dy - calorSizeXY/2;
  }
  else if(layer==CaloIdx::HCAL2){
	auto column = copyNo / kNofHad2Rows ;
	auto row = copyNo % kNofHad2Rows;
	X = (column+0.5)*HCAL2_dx - calorSizeXY/2;
	Y = (row+0.5)*HCAL2_dy - calorSizeXY/2;
  }   
  else if(layer==CaloIdx::HCAL3){
	auto column = copyNo / kNofHad3Rows ;
	auto row = copyNo % kNofHad3Rows;
	X = (column+0.5)*HCAL3_dx - calorSizeXY/2;
	Y = (row+0.5)*HCAL3_dy - calorSizeXY/2;
  }
  static G4double _XY[2] = {-1};
  _XY[0] = X; _XY[1]=Y;
  return _XY;
}

G4int CaloRCellMatrix::XYToCopy(G4double * XY, CaloIdx layer){
  G4int copyNo = -1;
  if(layer==CaloIdx::ECAL1){
	G4int column = G4int((XY[0]+calorSizeXY/2)/ECAL1_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/ECAL1_dy);
	copyNo = row + column*kNofEm1Rows;
  }
  else if(layer==CaloIdx::ECAL2){
	G4int column = G4int((XY[0]+calorSizeXY/2)/ECAL2_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/ECAL2_dy);
	copyNo = row + column*kNofEm2Rows;
  }   
  else if(layer==CaloIdx::ECAL3){
	G4int column = G4int((XY[0]+calorSizeXY/2)/ECAL3_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/ECAL3_dy);
	copyNo = row + column*kNofEm3Rows;
  }  
  else if(layer==CaloIdx::HCAL1){
	G4int column = G4int((XY[0]+calorSizeXY/2)/HCAL1_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/HCAL1_dy);
	copyNo = row + column*kNofHad1Rows;
  }
  else if(layer==CaloIdx::HCAL2){
	G4int column = G4int((XY[0]+calorSizeXY/2)/HCAL2_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/HCAL2_dy);
	copyNo = row + column*kNofHad2Rows;
  }   
  else if(layer==CaloIdx::HCAL3){
	G4int column = G4int((XY[0]+calorSizeXY/2)/HCAL3_dx);
	G4int row = G4int((XY[1]+calorSizeXY/2)/HCAL3_dy);
	copyNo = row + column*kNofHad3Rows;
  }  
  return copyNo;
}

