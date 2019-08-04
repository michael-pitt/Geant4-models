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
/// \file Geant4-models/ATLAS-simplified/src/RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "g4root.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloRRunAction::CaloRRunAction()
 : G4UserRunAction() 
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  G4String fileName = "CaloResponce";
  analysisManager->SetFileName(fileName);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloRRunAction::~CaloRRunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloRRunAction::BeginOfRunAction(const G4Run* myRun)
{ 

  //get event action (to be able to link vectors to Ntuple)
  const CaloREventAction* constEventAction 
		= static_cast<const CaloREventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  CaloREventAction* eventAction = const_cast<CaloREventAction*>(constEventAction);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  // Open an output file
  //

  analysisManager->OpenFile();
  
  // protection to avoid creating ntuple more times
  if ( myRun->GetRunID() == 0  ) {
	  analysisManager->CreateNtuple("physics", "Edep");
	  analysisManager->CreateNtupleDColumn("ECAL1_e");
	  analysisManager->CreateNtupleDColumn("ECAL2_e");
	  analysisManager->CreateNtupleDColumn("ECAL3_e");
	  analysisManager->CreateNtupleDColumn("HCAL1_e");
	  analysisManager->CreateNtupleDColumn("HCAL2_e");
	  analysisManager->CreateNtupleDColumn("HCAL3_e");
	  analysisManager->CreateNtupleDColumn("Cal_e");
	  analysisManager->CreateNtupleDColumn("Caltotal_e");
	  analysisManager->CreateNtupleDColumn("cell_e",eventAction->GetCellEnergyVec());  
	  analysisManager->CreateNtupleDColumn("cell_x",eventAction->GetCellXVec());  
	  analysisManager->CreateNtupleDColumn("cell_y",eventAction->GetCellYVec());  
	  analysisManager->CreateNtupleDColumn("cell_dx",eventAction->GetCelldXVec());  
	  analysisManager->CreateNtupleDColumn("cell_dy",eventAction->GetCelldYVec());  
	  analysisManager->CreateNtupleIColumn("cell_l",eventAction->GetCellLayerVec());  
	  analysisManager->CreateNtupleDColumn("cell_chfrac",eventAction->GetCellChFractionVec());  
	  analysisManager->CreateNtupleDColumn("particle_e",eventAction->GetParEnergyVec());  
	  analysisManager->CreateNtupleDColumn("particle_x",eventAction->GetParXVec());  
	  analysisManager->CreateNtupleDColumn("particle_y",eventAction->GetParYVec());  
	  analysisManager->CreateNtupleDColumn("particle_z",eventAction->GetParZVec());  
	  analysisManager->CreateNtupleDColumn("particle_px",eventAction->GetParPXVec());  
	  analysisManager->CreateNtupleDColumn("particle_py",eventAction->GetParPYVec());  
	  analysisManager->CreateNtupleDColumn("particle_pz",eventAction->GetParPZVec());  
	  analysisManager->CreateNtupleIColumn("particle_pdgId",eventAction->GetParPdgIdVec());  
	  analysisManager->FinishNtuple();
  }

}

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloRRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
