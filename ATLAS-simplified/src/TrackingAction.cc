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
/// \file Geant4-models/ATLAS-simplified/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

#include "TrackingAction.hh"
#include "Trajectory.hh"
#include "TrackInformation.hh"
#include "EventAction.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
CaloRTrackingAction::CaloRTrackingAction()
:G4UserTrackingAction()
{;}

void CaloRTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory for all tracks (is it necessary?)
  // /tracking/storeTrajectory _bool_ does the same
  fpTrackingManager->SetStoreTrajectory(true);

  // Add Origin particle info to the trajectory
  G4int ParentID = aTrack->GetParentID() ;
  auto trajectory = new CaloRTrajectory(aTrack);
  OriginCharge =  trajectory->GetOriginCharge();
  OriginPDGEncoding = trajectory->GetOriginPDGEncoding();
  
  if(0==ParentID){ // primary particle
    OriginCharge = trajectory->GetCharge();
    OriginPDGEncoding = trajectory->GetPDGEncoding();
  }
  else{
    auto event
      = static_cast<const G4Event*>(
          G4RunManager::GetRunManager()->GetCurrentEvent());
	G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();  
	G4int n_trajectories = 0;
	if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

	// Store Origin PdgID and charge from it's parrent:
	for(int iTraj = n_trajectories - 1; iTraj >-1 ; iTraj--) {
		auto trj = (CaloRTrajectory*)((*(event->GetTrajectoryContainer()))[iTraj]);
		if(ParentID==trj->GetTrackID()){
			OriginCharge = trj->GetOriginCharge();
			OriginPDGEncoding = trj->GetOriginPDGEncoding();
		}
	}
  }
  
  // Store new trajectory
  trajectory->SetOriginCharge(OriginCharge);
  trajectory->SetOriginPDGEncoding(OriginPDGEncoding);
  fpTrackingManager->SetTrajectory(trajectory);

  // update TrackInfo
  CaloRTrackInformation * trackInfo( static_cast< CaloRTrackInformation * >(
                                                aTrack->GetUserInformation() ) );

  if ( trackInfo ){
	  trackInfo->SetOriginCharge(OriginCharge);
	  trackInfo->SetOriginPDG(OriginPDGEncoding);
	  aTrack->SetUserInformation(trackInfo);
  } else {
	G4Track *  theTrack( const_cast< G4Track * >( aTrack ) );
	trackInfo = new CaloRTrackInformation();
	trackInfo->SetOriginCharge(OriginCharge);
	trackInfo->SetOriginPDG(OriginPDGEncoding);
	theTrack->SetUserInformation( trackInfo );  
  }

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void CaloRTrackingAction::PostUserTrackingAction(const G4Track* /*aTrack*/)
{}


