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
/// \file Geant4-models/ATLAS-simplified/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//

#include "EventAction.hh"
#include "RunAction.hh"
#include "g4root.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloREventAction::CaloREventAction()
 : G4UserEventAction(), Caltotal_e(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloREventAction::~CaloREventAction()
{
}

G4THitsMap<G4double>* 
CaloREventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("CaloREventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double CaloREventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;  
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CaloREventAction::SetCellInfo(G4double e_ch, G4double e_nu, G4double * XY, CaloIdx l)
{
  // Fill the corresponding vectors with the cell info
  cell_e.push_back(e_ch+e_nu);
  cell_x.push_back(XY[0]);
  cell_y.push_back(XY[1]);
  if(e_ch+e_nu) cell_Chfrac.push_back(e_ch/(e_ch+e_nu));
  else cell_Chfrac.push_back(-1);
  cell_l.push_back(G4int(l));
  if(l==CaloIdx::ECAL1){
    cell_dx.push_back(ECAL1_dx);
	cell_dy.push_back(ECAL1_dy);
  }
  else if(l==CaloIdx::ECAL2){
    cell_dx.push_back(ECAL2_dx);
	cell_dy.push_back(ECAL2_dy);
  }   
  else if(l==CaloIdx::ECAL3){
    cell_dx.push_back(ECAL3_dx);
	cell_dy.push_back(ECAL3_dy);
  }  
  else if(l==CaloIdx::HCAL1){
    cell_dx.push_back(HCAL1_dx);
	cell_dy.push_back(HCAL1_dy);
  }
  else if(l==CaloIdx::HCAL2){
    cell_dx.push_back(HCAL2_dx);
	cell_dy.push_back(HCAL2_dy);
  }   
  else if(l==CaloIdx::HCAL3){
    cell_dx.push_back(HCAL3_dx);
	cell_dy.push_back(HCAL3_dy);
  }    
}

void CaloREventAction::SetParInfo(G4ThreeVector x3, G4double e, G4ThreeVector p3, G4int pdgId)
{
  // Fill the corresponding vectors with the cell info
  particle_x.push_back(x3.x()/mm);
  particle_y.push_back(x3.y()/mm);
  particle_z.push_back(x3.z()/mm);
  particle_e.push_back(e/GeV);
  particle_px.push_back(p3.x()/GeV);
  particle_py.push_back(p3.y()/GeV);
  particle_pz.push_back(p3.z()/GeV);
  particle_pdgId.push_back(pdgId);
  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloREventAction::AddHit(G4double uncalib_de, G4ThreeVector pos, 
							  G4String c_name, G4double fCharge, G4int EMflag)
{
  // Fill the cell
  G4double XY[] = {pos.x(),pos.y()}; 
  G4double de = CalibrateCellEnergy(uncalib_de,c_name,EMflag);

  if(fCharge){
	if(c_name.contains("ECAL1"))	  
	  fChECAL1_CellE[XYToCopy(XY,CaloIdx::ECAL1)]+=de;
	if(c_name.contains("ECAL2"))
	  fChECAL2_CellE[XYToCopy(XY,CaloIdx::ECAL2)]+=de;
	if(c_name.contains("ECAL3"))
	  fChECAL3_CellE[XYToCopy(XY,CaloIdx::ECAL3)]+=de;
	if(c_name.contains("HCAL1"))
	  fChHCAL1_CellE[XYToCopy(XY,CaloIdx::HCAL1)]+=de;
	if(c_name.contains("HCAL2"))
	  fChHCAL2_CellE[XYToCopy(XY,CaloIdx::HCAL2)]+=de;
	if(c_name.contains("HCAL3"))
	  fChHCAL3_CellE[XYToCopy(XY,CaloIdx::HCAL3)]+=de;
  }
  else{
	if(c_name.contains("ECAL1"))
	  fNuECAL1_CellE[XYToCopy(XY,CaloIdx::ECAL1)]+=de;
	if(c_name.contains("ECAL2"))
	  fNuECAL2_CellE[XYToCopy(XY,CaloIdx::ECAL2)]+=de;
	if(c_name.contains("ECAL3"))
	  fNuECAL3_CellE[XYToCopy(XY,CaloIdx::ECAL3)]+=de;
	if(c_name.contains("HCAL1"))
	  fNuHCAL1_CellE[XYToCopy(XY,CaloIdx::HCAL1)]+=de;
	if(c_name.contains("HCAL2"))
	  fNuHCAL2_CellE[XYToCopy(XY,CaloIdx::HCAL2)]+=de;
	if(c_name.contains("HCAL3"))
	  fNuHCAL3_CellE[XYToCopy(XY,CaloIdx::HCAL3)]+=de;
  }
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CaloREventAction::BeginOfEventAction(const G4Event* event)
{
  // ResetVariables
  cell_e.clear();
  cell_x.clear();
  cell_y.clear();
  cell_dx.clear();
  cell_dy.clear();
  cell_l.clear();
  cell_Chfrac.clear();
  particle_x.clear();
  particle_y.clear();
  particle_z.clear();
  particle_e.clear();
  particle_px.clear();
  particle_py.clear();
  particle_pz.clear();
  particle_pdgId.clear();
  
  Caltotal_e = 0;

  for(G4int icopy = 0; icopy < kNofEm1Cells; icopy++){
      fChECAL1_CellE[icopy] = 0. ;
      fNuECAL1_CellE[icopy] = 0. ;
  }
  for(G4int icopy = 0; icopy < kNofEm2Cells; icopy++){
      fChECAL2_CellE[icopy] = 0. ;
      fNuECAL2_CellE[icopy] = 0. ;
  }    
  for(G4int icopy = 0; icopy < kNofEm3Cells; icopy++){
      fChECAL3_CellE[icopy] = 0. ;
      fNuECAL3_CellE[icopy] = 0. ;
  }
  for(G4int icopy = 0; icopy < kNofHad1Cells; icopy++){
      fChHCAL1_CellE[icopy] = 0. ;
      fNuHCAL1_CellE[icopy] = 0. ;
  }
  for(G4int icopy = 0; icopy < kNofHad2Cells; icopy++){
      fChHCAL2_CellE[icopy] = 0. ;
      fNuHCAL2_CellE[icopy] = 0. ;
  }    
  for(G4int icopy = 0; icopy < kNofHad3Cells; icopy++){
      fChHCAL3_CellE[icopy] = 0. ;
      fNuHCAL3_CellE[icopy] = 0. ;
  }  
  
  // Get information about the primary particles (might change during the run?)
  G4int nVtx= event-> GetNumberOfPrimaryVertex();
  for(auto i=0; i< nVtx; i++) {
	  const G4PrimaryVertex* primaryVertex= event-> GetPrimaryVertex(i);
	  auto nPar = primaryVertex->GetNumberOfParticle ();
	  for(auto j=0; j< nPar; j++) {
		auto particle = primaryVertex->GetPrimary(j);
		SetParInfo(primaryVertex->GetPosition(), particle->GetTotalEnergy(), particle->GetMomentum(), particle->GetPDGcode());
	  }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CaloREventAction::EndOfEventAction(const G4Event* /*event*/)
{
  
  //Store cells only above threshold
  G4double ECAL1Edep=0, ECAL2Edep=0, ECAL3Edep=0;
  G4double HCAL1Edep=0, HCAL2Edep=0, HCAL3Edep=0;
  
  for(G4int icopy = 0; icopy < kNofEm1Cells; icopy++){
	  if((fChECAL1_CellE[icopy]+fNuECAL1_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChECAL1_CellE[icopy],fNuECAL1_CellE[icopy],CopyToXY(icopy,CaloIdx::ECAL1),CaloIdx::ECAL1);
	  ECAL1Edep+=fChECAL1_CellE[icopy]+fNuECAL1_CellE[icopy];
  }}
  for(G4int icopy = 0; icopy < kNofEm2Cells; icopy++){
	  if((fChECAL2_CellE[icopy]+fNuECAL2_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChECAL2_CellE[icopy],fNuECAL2_CellE[icopy],CopyToXY(icopy,CaloIdx::ECAL2),CaloIdx::ECAL2);
	  ECAL2Edep+=fChECAL2_CellE[icopy]+fNuECAL2_CellE[icopy];
  }}
  for(G4int icopy = 0; icopy < kNofEm3Cells; icopy++){
	  if((fChECAL3_CellE[icopy]+fNuECAL3_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChECAL3_CellE[icopy],fNuECAL3_CellE[icopy],CopyToXY(icopy,CaloIdx::ECAL3),CaloIdx::ECAL3);
	  ECAL3Edep+=fChECAL3_CellE[icopy]+fNuECAL3_CellE[icopy];
  }}   
  for(G4int icopy = 0; icopy < kNofHad1Cells; icopy++){
	  if((fChHCAL1_CellE[icopy]+fNuHCAL1_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChHCAL1_CellE[icopy],fNuHCAL1_CellE[icopy],CopyToXY(icopy,CaloIdx::HCAL1),CaloIdx::HCAL1);
	  HCAL1Edep+=fChHCAL1_CellE[icopy]+fNuHCAL1_CellE[icopy];
  }}
  for(G4int icopy = 0; icopy < kNofHad2Cells; icopy++){
	  if((fChHCAL2_CellE[icopy]+fNuHCAL2_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChHCAL2_CellE[icopy],fNuHCAL2_CellE[icopy],CopyToXY(icopy,CaloIdx::HCAL2),CaloIdx::HCAL2);
	  HCAL2Edep+=fChHCAL2_CellE[icopy]+fNuHCAL2_CellE[icopy];
  }}
  for(G4int icopy = 0; icopy < kNofHad3Cells; icopy++){
	  if((fChHCAL3_CellE[icopy]+fNuHCAL3_CellE[icopy])>CELL_THRESHOLD){
	  SetCellInfo(fChHCAL3_CellE[icopy],fNuHCAL3_CellE[icopy],CopyToXY(icopy,CaloIdx::HCAL3),CaloIdx::HCAL3);
	  HCAL3Edep+=fChHCAL3_CellE[icopy]+fNuHCAL3_CellE[icopy];
  }}
  
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  // fill ntuple
  analysisManager->FillNtupleDColumn( 0, ECAL1Edep);
  analysisManager->FillNtupleDColumn( 1, ECAL2Edep);
  analysisManager->FillNtupleDColumn( 2, ECAL3Edep);
  analysisManager->FillNtupleDColumn( 3, HCAL1Edep);
  analysisManager->FillNtupleDColumn( 4, HCAL2Edep);
  analysisManager->FillNtupleDColumn( 5, HCAL3Edep);
  analysisManager->FillNtupleDColumn( 6, ECAL1Edep+ECAL2Edep+ECAL3Edep+HCAL1Edep+HCAL2Edep+HCAL3Edep);
  analysisManager->FillNtupleDColumn( 7, Caltotal_e);
  analysisManager->AddNtupleRow();

  //print per event (modulo n)
  //
  //auto eventID = event->GetEventID();
  //auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
  //  G4cout << "---> End of event: " << eventID << G4endl;    
  //}  

}


G4double CaloREventAction::CalibrateCellEnergy(G4double de, G4String c_name, G4int DynParflag)
{
 // Cell energy calibration includes Tile scintilator energy correction to hadrons (protons)
 // Responce to alphas and carbons ignored but the fractions coded in data/BC400_LYR.csv
 
  // cell EM calibration:
  if(c_name.contains("ECAL")){
	  de*=ECAL_cell_calib;
  }
  else{
	  de*=HCAL_cell_calib;
	  if(DynParflag==1) de*=HCAL_HE_RATIO(de);
	  else if(DynParflag==2) de*=HCAL_CE_RATIO(de);
  }
  return de;
}  

