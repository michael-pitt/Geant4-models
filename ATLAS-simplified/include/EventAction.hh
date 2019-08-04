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
/// \file Geant4-models/ATLAS-simplified/include/EventAction.hh
/// \brief Definition of the CaloREventAction class
//
//
#ifndef CALOR_EVENT_ACTION_H
#define CALOR_EVENT_ACTION_H

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "G4THitsMap.hh"
#include "Constants.hh"
#include "CellMatrix.hh"


class CaloREventAction : public G4UserEventAction, public CaloRCellMatrix {
public:
  CaloREventAction();
  ~CaloREventAction();

  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);
  std::vector<G4double>& GetCellEnergyVec() { return cell_e; };
  std::vector<G4double>& GetCellXVec() { return cell_x; };
  std::vector<G4double>& GetCelldXVec() { return cell_dx; };
  std::vector<G4double>& GetCellYVec() { return cell_y; };
  std::vector<G4double>& GetCelldYVec() { return cell_dy; };
  std::vector<G4int>& GetCellLayerVec() { return cell_l; };
  std::vector<G4double>& GetCellChFractionVec() { return cell_Chfrac; };
  std::vector<G4double>& GetParEnergyVec() { return particle_e; };
  std::vector<G4double>& GetParPXVec() { return particle_px; };
  std::vector<G4double>& GetParPYVec() { return particle_py; };
  std::vector<G4double>& GetParPZVec() { return particle_pz; };
  std::vector<G4double>& GetParXVec() { return particle_x; };
  std::vector<G4double>& GetParYVec() { return particle_y; };
  std::vector<G4double>& GetParZVec() { return particle_z; };
  std::vector<G4int>& GetParPdgIdVec() { return particle_pdgId; };
  
  void AddHit(G4double de, G4ThreeVector pos, G4String c_name, G4double fCharge, G4int EMflag);
  void AddEnergyTotal(G4double de){ Caltotal_e += de; };
 
private:
  // methods
  G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
                                          const G4Event* event) const;
  G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
  void SetCellInfo(G4double e_ch, G4double e_nu, G4double *XY, CaloIdx l);
  void SetParInfo(G4ThreeVector x3, G4double e, G4ThreeVector p3, G4int pdgId);
  G4double CalibrateCellEnergy(G4double de,G4String c_name,G4int EMflag);

  // data members                   
  G4double  Caltotal_e;
  std::vector<G4double> cell_e;
  std::vector<G4double> cell_x;
  std::vector<G4double> cell_dx;
  std::vector<G4double> cell_y;
  std::vector<G4double> cell_dy;
  std::vector<G4int> cell_l;
  std::vector<G4double> cell_Chfrac;
  std::vector<G4double> particle_e;
  std::vector<G4double> particle_x;
  std::vector<G4double> particle_y;
  std::vector<G4double> particle_z;
  std::vector<G4double> particle_px;
  std::vector<G4double> particle_py;
  std::vector<G4double> particle_pz;
  std::vector<G4int> particle_pdgId;
  
  // Cell Matrix:
  G4double fChECAL1_CellE[kNofEm1Cells];
  G4double fChECAL2_CellE[kNofEm2Cells];
  G4double fChECAL3_CellE[kNofEm3Cells];
  G4double fChHCAL1_CellE[kNofHad1Cells];
  G4double fChHCAL2_CellE[kNofHad2Cells];
  G4double fChHCAL3_CellE[kNofHad3Cells]; 
  G4double fNuECAL1_CellE[kNofEm1Cells];
  G4double fNuECAL2_CellE[kNofEm2Cells];
  G4double fNuECAL3_CellE[kNofEm3Cells];
  G4double fNuHCAL1_CellE[kNofHad1Cells];
  G4double fNuHCAL2_CellE[kNofHad2Cells];
  G4double fNuHCAL3_CellE[kNofHad3Cells]; 
  
	  
};

#endif
