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
/// \file Geant4-models/ATLAS-simplified/src/DetectorConstruction.cc
/// \brief Implementation of the CaloRDetectorConstruction class
//

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "DetectorConstruction.hh"
#include "CellParameterisation.hh"

// Sensitive Detector managers
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

G4ThreadLocal
G4GlobalMagFieldMessenger* CaloRDetectorConstruction::fMagFieldMessenger = 0; 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRDetectorConstruction::CaloRDetectorConstruction()
 : G4VUserDetectorConstruction() ,
 fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRDetectorConstruction::~CaloRDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* CaloRDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloRDetectorConstruction::DefineMaterials()
{
  // ==============================================================
  // Materials
  // ==============================================================

  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_lAr");
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); //PolyVinylToluene, C_9H_10

  // Vacuum
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density;
  G4int ncomponents, natoms; 
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= CLHEP::universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
  
  // Basic elements:
  G4Element* elH  = nistManager->FindOrBuildElement(1);
  G4Element* elC  = nistManager->FindOrBuildElement(6);
  
  //Compounds
  G4Material* Sci = new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(elC, natoms=9);
  Sci->AddElement(elH, natoms=10);
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);  
   				  
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* CaloRDetectorConstruction::DefineVolumes()
{

  // Geometry parameters are defined in CaloRConstants.hh
  // G4double  calorSizeXY  = ...
  // G4double  calorThickness  = ...

  auto  worldSizeXY = 1.5 * calorSizeXY;
  auto  worldSizeZ  = 1.2 * (calorThickness + GunDinsance); 
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto ECALabs = G4Material::GetMaterial("G4_Pb");
  auto ECALgap = G4Material::GetMaterial("G4_lAr");
  auto HCALabs = G4Material::GetMaterial("G4_Fe");
  auto HCALgap = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  if ( ! defaultMaterial || ! HCALabs || ! HCALgap || ! ECALabs || ! ECALgap ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("CaloRDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }
  
  G4int counter=0;
  G4int SkipRowECAL[nSubCellDivECAL]; G4int SkipRowHCAL[nSubCellDivHCAL];
  G4double runECAL_dY[nSubCellDivECAL]; G4double runHCAL_dY[nSubCellDivHCAL];
  for(int i=0;i<(nSubCellDivECAL/2);++i){
    SkipRowECAL[i]=i;SkipRowECAL[nSubCellDivECAL-i-1]=i+1;
    runECAL_dY[i]=dECAL_abs*i;runECAL_dY[nSubCellDivECAL-i-1]=dECAL_abs*(i+1); 
  }
  SkipRowHCAL[0] = 1; SkipRowHCAL[1] = 4;
  for(int i=0;i<(nSubCellDivHCAL);++i){
    runHCAL_dY[i]=dHCAL_gap*SkipRowHCAL[i]; 
  }

  // ==============================================================
  // Experimental Hall (world)
  // ==============================================================
  auto worldS=
    new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2 );

  auto worldLV=
    new G4LogicalVolume(worldS, defaultMaterial, "World");

  auto worldPV= new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  //                               
  // Calorimeters
  //

  auto calorimeterS=
    new G4Box("Calorimeter", calorSizeXY/2, calorSizeXY/2, calorThickness/2 );

  auto calorLV=
    new G4LogicalVolume(calorimeterS, defaultMaterial, "Calorimeter");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,0,(GunDinsance-calorThickness)/2),  // its position
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);   // checking overlaps 
	
  //                               
  // ECAL
  //  
  auto ECALS
    = new G4Box("ECAL",     // its name
                 calorSizeXY/2, calorSizeXY/2, ECALThickness/2); // its size
                         
  auto ECALLV
    = new G4LogicalVolume(
                 ECALS,    // its solid
                 defaultMaterial, // its material
                 "ECALLV");  // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -HCALThickness/2),  // its position
                 ECALLV,           // its logical volume                         
                 "ECAL",           // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
				 
  //                                 
  // ECAL Layers
  //
  auto ECAL1S = new G4Box("ECALLayer1", calorSizeXY/2, calorSizeXY/2, ECAL1_dz/2);               
  auto ECAL1LV = new G4LogicalVolume(ECAL1S,defaultMaterial,"ECALLayer1LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(ECAL2_dz+ECAL3_dz)/2),  // its position
                 ECAL1LV,          // its logical volume                         
                 "ECALLayer1",     // its name
                 ECALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
				 
  auto ECAL2S = new G4Box("ECALLayer2", calorSizeXY/2, calorSizeXY/2, ECAL2_dz/2);               
  auto ECAL2LV = new G4LogicalVolume(ECAL2S,defaultMaterial,"ECALLayer2LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(ECAL3_dz-ECAL1_dz)/2),  // its position
                 ECAL2LV,          // its logical volume                         
                 "ECALLayer2",     // its name
                 ECALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
				 
  auto ECAL3S = new G4Box("ECALLayer3", calorSizeXY/2, calorSizeXY/2, ECAL3_dz/2);               
  auto ECAL3LV = new G4LogicalVolume(ECAL3S,defaultMaterial,"ECALLayer3LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (ECAL1_dz+ECAL2_dz)/2),  // its position
                 ECAL3LV,          // its logical volume                         
                 "ECALLayer3",     // its name
                 ECALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
				 
  //                               
  // ECAL Cells
  //	 
  auto ECAL1cellSolid = new G4Box("ECAL1cell",ECAL1_dx/2,ECAL1_dy/2,ECAL1_dz/2);
  auto ECAL1cellLV = new G4LogicalVolume(ECAL1cellSolid,defaultMaterial,"ECAL1cellLogical");
  new G4PVParameterised("ECAL1cellPhysical",ECAL1cellLV,ECAL1LV,
                        kXAxis,kNofEm1Cells,new Em1CellParameterisation()); 
  
  auto dZECAL_layer = X0ECAL/ECAL_nLayerPerX0;
  auto ECAL1subcellSolid = new G4Box("ECAL1subcell",ECAL1_dx/2,(dECAL_abs+dECAL_gap)/2,dZECAL_layer/2);
  auto ECAL1subcellLV = new G4LogicalVolume(ECAL1subcellSolid,defaultMaterial,"ECAL1subcellLogical");
  new G4PVParameterised("ECAL1subcellPhysical",ECAL1subcellLV,ECAL1cellLV,
                        kYAxis,kNofEm1SubCells,new Em1SubCellParameterisation());
  auto ECAL1subcellAbs = new G4Box("ECAL1subcell"+absWord,ECAL1_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL1subcellAbsLV = new G4LogicalVolume(ECAL1subcellAbs,ECALabs,"ECAL1subcell"+absWord);
  for(int i=0;i<nSubCellDivECAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dECAL_gap/2+runECAL_dY[i], -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
                 ECAL1subcellAbsLV,"ECAL1subcell"+absWord,ECAL1subcellLV,false,i,false);}
  auto ECAL1subcellGap = new G4Box("ECAL1subcell"+gapWord,ECAL1_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL1subcellGapLV = new G4LogicalVolume(ECAL1subcellGap,ECALgap,"ECAL1subcell"+gapWord);
  counter=0;
  for(int i=0;i<nSubCellDivECAL;i++){for(int row=0;row<nSubCellDivECAL2;row++){if(SkipRowECAL[i]==row) continue; counter++;
	            new G4PVPlacement(0,G4ThreeVector(0., row*dECAL_abs - dECAL_gap/2, -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
				ECAL1subcellGapLV,"ECAL1subcell"+gapWord,ECAL1subcellLV,false,counter,false);}}

  auto ECAL2cellSolid = new G4Box("ECAL2cell",ECAL2_dx/2,ECAL2_dy/2,ECAL2_dz/2);
  auto ECAL2cellLV = new G4LogicalVolume(ECAL2cellSolid,defaultMaterial,"ECAL2cellLogical");
  new G4PVParameterised("ECAL2cellPhysical",ECAL2cellLV,ECAL2LV,
                        kXAxis,kNofEm2Cells,new Em2CellParameterisation());  
  						
  auto ECAL2subcellSolid = new G4Box("ECAL2subcell",ECAL2_dx/2,(dECAL_abs+dECAL_gap)/2,dZECAL_layer/2);
  auto ECAL2subcellLV = new G4LogicalVolume(ECAL2subcellSolid,defaultMaterial,"ECAL2subcellLogical");
  new G4PVParameterised("ECAL2subcellPhysical",ECAL2subcellLV,ECAL2cellLV,
                        kYAxis,kNofEm2SubCells,new Em2SubCellParameterisation());
  auto ECAL2subcellAbs = new G4Box("ECAL2subcell"+absWord,ECAL2_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL2subcellAbsLV = new G4LogicalVolume(ECAL2subcellAbs,ECALabs,"ECAL2subcell"+absWord);
  for(int i=0;i<nSubCellDivECAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dECAL_gap/2+runECAL_dY[i], -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
                 ECAL2subcellAbsLV,"ECAL2subcell"+absWord,ECAL2subcellLV,false,i,false);}
  auto ECAL2subcellGap = new G4Box("ECAL2subcell"+gapWord,ECAL2_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL2subcellGapLV = new G4LogicalVolume(ECAL2subcellGap,ECALgap,"ECAL2subcell"+gapWord);
  counter=0;
  for(int i=0;i<nSubCellDivECAL;i++){for(int row=0;row<nSubCellDivECAL2;row++){if(SkipRowECAL[i]==row) continue; counter++;
	            new G4PVPlacement(0,G4ThreeVector(0., row*dECAL_abs - dECAL_gap/2, -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
				ECAL2subcellGapLV,"ECAL2subcell"+gapWord,ECAL2subcellLV,false,counter,false);}}
  
  auto ECAL3cellSolid = new G4Box("ECAL3cell",ECAL3_dx/2,ECAL3_dy/2,ECAL3_dz/2);
  auto ECAL3cellLV = new G4LogicalVolume(ECAL3cellSolid,defaultMaterial,"ECAL3cellLogical");
  new G4PVParameterised("ECAL3cellPhysical",ECAL3cellLV,ECAL3LV,
                        kXAxis,kNofEm3Cells,new Em3CellParameterisation());  
  auto ECAL3subcellSolid = new G4Box("ECAL3subcell",ECAL3_dx/2,(dECAL_abs+dECAL_gap)/2,dZECAL_layer/2);
  auto ECAL3subcellLV = new G4LogicalVolume(ECAL3subcellSolid,defaultMaterial,"ECAL3subcellLogical");
  new G4PVParameterised("ECAL3subcellPhysical",ECAL3subcellLV,ECAL3cellLV,
                        kYAxis,kNofEm3SubCells,new Em3SubCellParameterisation());
  auto ECAL3subcellAbs = new G4Box("ECAL3subcell"+absWord,ECAL3_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL3subcellAbsLV = new G4LogicalVolume(ECAL3subcellAbs,ECALabs,"ECAL3subcell"+absWord);
  for(int i=0;i<nSubCellDivECAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dECAL_gap/2+runECAL_dY[i], -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
                 ECAL3subcellAbsLV,"ECAL3subcell"+absWord,ECAL3subcellLV,false,i,false);}
  auto ECAL3subcellGap = new G4Box("ECAL3subcell"+gapWord,ECAL3_dx/2,dECAL_abs/2,dZECAL_layer/2/nSubCellDivECAL);
  auto ECAL3subcellGapLV = new G4LogicalVolume(ECAL3subcellGap,ECALgap,"ECAL3subcell"+gapWord);
  counter=0;
  for(int i=0;i<nSubCellDivECAL;i++){for(int row=0;row<nSubCellDivECAL2;row++){if(SkipRowECAL[i]==row) continue; counter++;
	            new G4PVPlacement(0,G4ThreeVector(0., row*dECAL_abs - dECAL_gap/2, -(nSubCellDivECAL-1-2*i)*dZECAL_layer/2/nSubCellDivECAL),
				ECAL3subcellGapLV,"ECAL3subcell"+gapWord,ECAL3subcellLV,false,counter,false);}}
  
  //                               
  // HCAL
  //  
  auto HCALS
    = new G4Box("HCAL",     // its name
                 calorSizeXY/2, calorSizeXY/2, HCALThickness/2); // its size
                         
  auto HCALLV
    = new G4LogicalVolume(
                 HCALS,    // its solid
                 defaultMaterial, // its material
                 "HCALLV");  // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., ECALThickness/2),  // its position
                 HCALLV,           // its logical volume                         
                 "HCAL",           // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // HCAL Layers
  //
  auto HCAL1S = new G4Box("HCALLayer1", calorSizeXY/2, calorSizeXY/2, HCAL1_dz/2);                     
  auto HCAL1LV = new G4LogicalVolume(HCAL1S,defaultMaterial,"HCALLayer1LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(HCAL2_dz+HCAL3_dz)/2),  // its position
                 HCAL1LV,          // its logical volume                         
                 "HCALLayer1",     // its name
                 HCALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  auto HCAL2S = new G4Box("HCALLayer2", calorSizeXY/2, calorSizeXY/2, HCAL2_dz/2);                     
  auto HCAL2LV = new G4LogicalVolume(HCAL2S,defaultMaterial,"HCALLayer2LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(HCAL3_dz-HCAL1_dz)/2),  // its position
                 HCAL2LV,          // its logical volume                         
                 "HCALLayer2",     // its name
                 HCALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
				 
  auto HCAL3S = new G4Box("HCALLayer3", calorSizeXY/2, calorSizeXY/2, HCAL3_dz/2);                     
  auto HCAL3LV = new G4LogicalVolume(HCAL3S,defaultMaterial,"HCALLayer3LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (HCAL1_dz+HCAL2_dz)/2),  // its position
                 HCAL3LV,          // its logical volume                         
                 "HCALLayer3",     // its name
                 HCALLV,           // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  

  //                               
  // HCAL Cells
  //
  auto HCAL1cellSolid = new G4Box("HCAL1cell",HCAL1_dx/2,HCAL1_dy/2,HCAL1_dz/2);
  auto HCAL1cellLV = new G4LogicalVolume(HCAL1cellSolid,defaultMaterial,"HCAL1cellLogical");
  new G4PVParameterised("HCAL1cellPhysical",HCAL1cellLV,HCAL1LV,
                        kXAxis,kNofHad1Cells,new Had1CellParameterisation());  
  auto HCAL1subcellSolid = new G4Box("HCAL1subcell",HCAL1_dx/2,(dHCAL_abs+dHCAL_gap)/2,dZHCAL1_subcell/2);
  auto HCAL1subcellLV = new G4LogicalVolume(HCAL1subcellSolid,defaultMaterial,"HCAL1subcellLogical");
  new G4PVParameterised("HCAL1subcellPhysical",HCAL1subcellLV,HCAL1cellLV,
                        kYAxis,kNofHad1SubCells,new Had1SubCellParameterisation());					
  auto HCAL1subcellGap = new G4Box("HCAL1subcell"+gapWord,HCAL1_dx/2,dHCAL_gap/2,dZHCAL1_subcell/2/nSubCellDivHCAL);
  auto HCAL1subcellGapLV = new G4LogicalVolume(HCAL1subcellGap,HCALgap,"HCAL1subcell"+gapWord);
  for(int i=0;i<nSubCellDivHCAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2+runHCAL_dY[i], -(nSubCellDivHCAL-1-2*i)*dZHCAL1_subcell/2/nSubCellDivHCAL),
				HCAL1subcellGapLV,"HCAL1subcell"+gapWord,HCAL1subcellLV,false,i,false);}
  counter=0;
  auto HCAL1subcellAbs = new G4Box("HCAL1subcell"+absWord,HCAL1_dx/2,dHCAL_gap/2,dZHCAL1_subcell/2/nSubCellDivHCAL);
  auto HCAL1subcellAbsLV = new G4LogicalVolume(HCAL1subcellAbs,HCALabs,"HCAL1subcell"+absWord);
  for(int i=0;i<nSubCellDivHCAL;i++){for(int row=0;row<nSubCellDivHCAL2;row++){if(SkipRowHCAL[i]==row) continue; counter++;
				new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2 + row*dHCAL_gap, -(nSubCellDivHCAL-1-2*i)*dZHCAL1_subcell/2/nSubCellDivHCAL),
                HCAL1subcellAbsLV,"HCAL1subcell"+absWord,HCAL1subcellLV,false,i,false);
  }}

  auto HCAL2cellSolid = new G4Box("HCAL2cell",HCAL2_dx/2,HCAL2_dy/2,HCAL2_dz/2);
  auto HCAL2cellLV = new G4LogicalVolume(HCAL2cellSolid,defaultMaterial,"HCAL2cellLogical");
  new G4PVParameterised("HCAL2cellPhysical",HCAL2cellLV,HCAL2LV,
                        kXAxis,kNofHad2Cells,new Had2CellParameterisation());  
  auto HCAL2subcellSolid = new G4Box("HCAL2subcell",HCAL2_dx/2,(dHCAL_abs+dHCAL_gap)/2,dZHCAL2_subcell/2);
  auto HCAL2subcellLV = new G4LogicalVolume(HCAL2subcellSolid,defaultMaterial,"HCAL2subcellLogical");
  new G4PVParameterised("HCAL2subcellPhysical",HCAL2subcellLV,HCAL2cellLV,
                        kYAxis,kNofHad2SubCells,new Had2SubCellParameterisation());					
  auto HCAL2subcellGap = new G4Box("HCAL2subcell"+gapWord,HCAL2_dx/2,dHCAL_gap/2,dZHCAL2_subcell/2/nSubCellDivHCAL);
  auto HCAL2subcellGapLV = new G4LogicalVolume(HCAL2subcellGap,HCALgap,"HCAL2subcell"+gapWord);
  for(int i=0;i<nSubCellDivHCAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2+runHCAL_dY[i], -(nSubCellDivHCAL-1-2*i)*dZHCAL2_subcell/2/nSubCellDivHCAL),
				HCAL2subcellGapLV,"HCAL2subcell"+gapWord,HCAL2subcellLV,false,i,false);}
  counter=0;
  auto HCAL2subcellAbs = new G4Box("HCAL2subcell"+absWord,HCAL2_dx/2,dHCAL_gap/2,dZHCAL2_subcell/2/nSubCellDivHCAL);
  auto HCAL2subcellAbsLV = new G4LogicalVolume(HCAL2subcellAbs,HCALabs,"HCAL2subcell"+absWord);
  for(int i=0;i<nSubCellDivHCAL;i++){for(int row=0;row<nSubCellDivHCAL2;row++){if(SkipRowHCAL[i]==row) continue; counter++;
				new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2 + row*dHCAL_gap, -(nSubCellDivHCAL-1-2*i)*dZHCAL2_subcell/2/nSubCellDivHCAL),
                HCAL2subcellAbsLV,"HCAL2subcell"+absWord,HCAL2subcellLV,false,i,false);
  }}

  auto HCAL3cellSolid = new G4Box("HCAL3cell",HCAL3_dx/2,HCAL3_dy/2,HCAL3_dz/2);
  auto HCAL3cellLV = new G4LogicalVolume(HCAL3cellSolid,defaultMaterial,"HCAL3cellLogical");
  new G4PVParameterised("HCAL3cellPhysical",HCAL3cellLV,HCAL3LV,
                        kXAxis,kNofHad3Cells,new Had3CellParameterisation());  
  auto HCAL3subcellSolid = new G4Box("HCAL3subcell",HCAL3_dx/2,(dHCAL_abs+dHCAL_gap)/2,dZHCAL3_subcell/2);
  auto HCAL3subcellLV = new G4LogicalVolume(HCAL3subcellSolid,defaultMaterial,"HCAL3subcellLogical");
  new G4PVParameterised("HCAL3subcellPhysical",HCAL3subcellLV,HCAL3cellLV,
                        kYAxis,kNofHad3SubCells,new Had3SubCellParameterisation());	
  auto HCAL3subcellGap = new G4Box("HCAL3subcell"+gapWord,HCAL3_dx/2,dHCAL_gap/2,dZHCAL3_subcell/2/nSubCellDivHCAL);
  auto HCAL3subcellGapLV = new G4LogicalVolume(HCAL3subcellGap,HCALgap,"HCAL3subcell"+gapWord);
  for(int i=0;i<nSubCellDivHCAL;i++){new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2+runHCAL_dY[i], -(nSubCellDivHCAL-1-2*i)*dZHCAL3_subcell/2/nSubCellDivHCAL),
				HCAL3subcellGapLV,"HCAL3subcell"+gapWord,HCAL3subcellLV,false,i,false);}
  counter=0;
  auto HCAL3subcellAbs = new G4Box("HCAL3subcell"+absWord,HCAL3_dx/2,dHCAL_gap/2,dZHCAL3_subcell/2/nSubCellDivHCAL);
  auto HCAL3subcellAbsLV = new G4LogicalVolume(HCAL3subcellAbs,HCALabs,"HCAL3subcell"+absWord);
  for(int i=0;i<nSubCellDivHCAL;i++){for(int row=0;row<nSubCellDivHCAL2;row++){if(SkipRowHCAL[i]==row) continue; counter++;
				new G4PVPlacement(0,G4ThreeVector(0., -dHCAL_abs/2 + row*dHCAL_gap, -(nSubCellDivHCAL-1-2*i)*dZHCAL3_subcell/2/nSubCellDivHCAL),
                HCAL3subcellAbsLV,"HCAL3subcell"+absWord,HCAL3subcellLV,false,i,false);
  }}
  
  //
  // print parameters
  //
  
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> Detector size : "
    << calorSizeXY/cm << " x " << calorSizeXY/cm << " x " << calorThickness/cm << " [ cm^3 ] " << G4endl
    << "------------------------------------------------------------" << G4endl 
    << "---> ECAL depth : "
    << ECALThickness/cm << " cm, " << ECALThickness/X0ECAL
	<< " X0 " << ECALThickness/LintECAL << " Lint " << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> HCAL depth : "
    << HCALThickness/cm << " cm, " << HCALThickness/LintHCAL << " Lint " << G4endl
    << "------------------------------------------------------------" << G4endl;  
	
G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> ECAL layers : " << G4endl
    << "L1 = " << ECAL_nLayerPerX0*ECAL1_X0 << " sub-layers, total=" << ECAL1_dz/X0ECAL << " X0 " << ECAL1_dz/LintECAL << " Lint " << G4endl
    << "L2 = " << ECAL_nLayerPerX0*ECAL2_X0 << " sub-layers, total=" << ECAL2_dz/X0ECAL << " X0 " << ECAL2_dz/LintECAL << " Lint " << G4endl
    << "L3 = " << ECAL_nLayerPerX0*ECAL3_X0 << " sub-layers, total=" << ECAL3_dz/X0ECAL << " X0 " << ECAL3_dz/LintECAL << " Lint " << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> HCAL layers : " << G4endl
    << "L1 = " <<  HCAL_nLayerPerLint*G4int(HCAL1_Lint) << " sub-layers, total=" << HCAL1_dz/LintHCAL << " Lint " << G4endl
    << "L2 = " <<  HCAL_nLayerPerLint*G4int(HCAL2_Lint) << " sub-layers, total=" << HCAL2_dz/LintHCAL << " Lint " << G4endl
    << "L3 = " <<  HCAL_nLayerPerLint*G4int(HCAL3_Lint) << " sub-layers, total=" << HCAL3_dz/LintHCAL << " Lint " << G4endl
    << "------------------------------------------------------------" << G4endl;  	

  
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  calorLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  ECALLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  HCALLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  ECAL1cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  ECAL2cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  ECAL3cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());   
  ECAL1subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  ECAL2subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  ECAL3subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  

  ECAL1subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  ECAL1subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  ECAL2subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  ECAL2subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  ECAL3subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  ECAL3subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 

  HCAL1cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  HCAL2cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  HCAL3cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());   
  HCAL1subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  HCAL2subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  
  HCAL3subcellLV->SetVisAttributes(G4VisAttributes::GetInvisible());  

  HCAL1subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  HCAL1subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  HCAL2subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  HCAL2subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  HCAL3subcellAbsLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  HCAL3subcellGapLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 

  //auto redHollow= new G4VisAttributes(true, G4Colour(1.0, 0.0, 0.0));
  //auto greenHollow= new G4VisAttributes(true, G4Colour(0.0, 1.0, 0.0));
  //auto blueSolid= new G4VisAttributes(true, G4Colour(0.33, 1.0, 1.0)); blueSolid->SetForceSolid();
  //auto mustardSolid= new G4VisAttributes(true, G4Colour(0.66, 0.66, 0.0)); mustardSolid->SetForceSolid();
  //ECAL1cellLV->SetVisAttributes(redHollow);
  //ECAL2cellLV->SetVisAttributes(redHollow);
  //ECAL3cellLV->SetVisAttributes(redHollow);
  //ECAL1subcellLV->SetVisAttributes(greenHollow);
  //ECAL2subcellLV->SetVisAttributes(greenHollow);
  //ECAL3subcellLV->SetVisAttributes(greenHollow);
  //ECAL1subcellGapLV->SetVisAttributes(blueSolid);
  //ECAL1subcellAbsLV->SetVisAttributes(mustardSolid);
  //ECAL1subcellGapLV->SetVisAttributes(blueSolid);
  //ECAL3subcellGapLV->SetVisAttributes(blueSolid);
  //HCAL2cellLV->SetVisAttributes(redHollow);
  //HCAL3cellLV->SetVisAttributes(redHollow);
  
  auto ECAL_VisAtt= new G4VisAttributes(true, G4Colour(0.0, 1.0, 1.0)); 
  ECAL1LV-> SetVisAttributes(ECAL_VisAtt);
  ECAL2LV-> SetVisAttributes(ECAL_VisAtt);
  ECAL3LV-> SetVisAttributes(ECAL_VisAtt);
  
  auto HCAL_VisAtt= new G4VisAttributes(true, G4Colour(1.0, 0.0, 1.0)); 
  HCAL1LV-> SetVisAttributes(HCAL_VisAtt);
  HCAL2LV-> SetVisAttributes(HCAL_VisAtt);
  HCAL3LV-> SetVisAttributes(HCAL_VisAtt);
  

  //
  // Always return the physical World
  //
  return worldPV;
}

void CaloRDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  // 
  // Scorers
  //
  // declare ECAL/HCAL as a MultiFunctionalDetector scorer
  //  
  /* 
  auto ECALDetector1 = new G4MultiFunctionalDetector("Calorimeter");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector1);
  ECALDetector1->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("Calorimeter",ECALDetector1);

  auto ECALDetector2 = new G4MultiFunctionalDetector("ECAL2");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector2);
  ECALDetector2->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("ECAL2subcell"+gapWord,ECALDetector2);
  auto ECALDetector3 = new G4MultiFunctionalDetector("ECAL3");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector3);
  ECALDetector3->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("ECAL3subcell"+gapWord,ECALDetector3);
  
  auto HCALDetector1 = new G4MultiFunctionalDetector("HCAL1");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector1);
  HCALDetector1->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL1subcell"+gapWord,HCALDetector1);
  auto HCALDetector2 = new G4MultiFunctionalDetector("HCAL2");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector2);
  HCALDetector2->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL2subcell"+gapWord,HCALDetector2);
  auto HCALDetector3 = new G4MultiFunctionalDetector("HCAL3");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector3);
  HCALDetector3->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL3subcell"+gapWord,HCALDetector3);
  
  // Collect energy from the abrorbers
  auto ECALDetector1abs = new G4MultiFunctionalDetector("ECAL1abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector1abs);
  ECALDetector1abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("ECAL1subcell"+absWord,ECALDetector1abs);
  auto ECALDetector2abs = new G4MultiFunctionalDetector("ECAL2abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector2abs);
  ECALDetector2abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("ECAL2subcell"+absWord,ECALDetector2abs);
  auto ECALDetector3abs = new G4MultiFunctionalDetector("ECAL3abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(ECALDetector3abs);
  ECALDetector3abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("ECAL3subcell"+absWord,ECALDetector3abs);
  
  auto HCALDetector1abs = new G4MultiFunctionalDetector("HCAL1abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector1abs);
  HCALDetector1abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL1subcell"+absWord,HCALDetector1abs);
  auto HCALDetector2abs = new G4MultiFunctionalDetector("HCAL2abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector2abs);
  HCALDetector2abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL2subcell"+absWord,HCALDetector2abs);
  auto HCALDetector3abs = new G4MultiFunctionalDetector("HCAL3abs");
  G4SDManager::GetSDMpointer()->AddNewDetector(HCALDetector3abs);
  HCALDetector3abs->RegisterPrimitive(new G4PSEnergyDeposit("Edep"));
  SetSensitiveDetector("HCAL3subcell"+absWord,HCALDetector3abs);  
  */
  
  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

