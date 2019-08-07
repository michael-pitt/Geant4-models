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
/// \file Geant4-models/ATLAS-simplified/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloRSteppingAction::CaloRSteppingAction(
                      const CaloRDetectorConstruction* detectorConstruction,
                      CaloREventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloRSteppingAction::~CaloRSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloRSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  edep -= step->GetNonIonizingEnergyDeposit(); // make difference of 0.5%
  if(edep<=0) return; // if no energy deposited skip calculations
   
  // get initial step point (the post step point if fall on boundary will be defined in the next volume)
  auto preStep = step->GetPreStepPoint();
  
  // Cell to step association
  auto touchable = preStep->GetTouchable(); //GetTouchableHandle();
  if(touchable->GetHistoryDepth() < fDetConstruction->SensitiveToCell) return;

  // fill total energy deposit variable (used for calibration)
  fEventAction->AddEnergyTotal(edep);

  // check if step was generated in the gap (reject otherwise):
  if(touchable->GetVolume()->GetName().contains(fDetConstruction->absWord)) return;
  
  // Trace back to the initial particle:
  G4Track * track( step->GetTrack() );
  CaloRTrackInformation *  trackInfo( static_cast< CaloRTrackInformation * >(track->GetUserInformation() ) );
  G4double fCharge = trackInfo->GetOriginCharge();
  
  // For HCAL check the particel type (different light yields for e/p/carbons
  // assume LY(e) = LY(gamma)
  G4int DynParPDG = track->GetDynamicParticle()->GetPDGcode();
  G4int DynParflag = 1;	// protons
  if((DynParPDG==22) || (DynParPDG==11) || (DynParPDG==-11)) DynParflag = 0; // electrons
  else if((DynParPDG>10e8)) DynParflag = 2; // heavy nuclei
	  
  // Store hit in cell Matrix
  auto cell = touchable->GetVolume(fDetConstruction->SensitiveToCell);
  auto position = cell->GetObjectTranslation(); // or GetTranslation? looks like the same
  //G4cout << " energy deposit = " << edep << " particle type = " << DynParPDG << "("<<DynParflag<<") hit in " << cell->GetName() << G4endl;
  fEventAction->AddHit(edep,position,cell->GetName(),fCharge, DynParflag);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
