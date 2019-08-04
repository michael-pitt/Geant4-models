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
/// \file Geant4-models/ATLAS-simplified/src/TrackInformation.cc
/// \brief Implementation of the TrackInformation class
//
//
//

#include "TrackInformation.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"    

G4ThreadLocal G4Allocator<CaloRTrackInformation> *
                                   aTrackInformationAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRTrackInformation::CaloRTrackInformation()
  : G4VUserTrackInformation()
{
    fOriginCharge = -999;
    fOriginPDG = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRTrackInformation::CaloRTrackInformation(const G4Track* aTrack)
  : G4VUserTrackInformation()
{
    fOriginCharge = aTrack->GetDynamicParticle()->GetCharge();
    fOriginPDG = aTrack->GetDynamicParticle()->GetPDGcode();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRTrackInformation
::CaloRTrackInformation(const CaloRTrackInformation* aTrackInfo)
  : G4VUserTrackInformation()
{
    fOriginCharge = aTrackInfo->fOriginCharge;
    fOriginPDG = aTrackInfo->fOriginPDG;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRTrackInformation::~CaloRTrackInformation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRTrackInformation& CaloRTrackInformation
::operator =(const CaloRTrackInformation& aTrackInfo)
{
    fOriginCharge = aTrackInfo.fOriginCharge;
    fOriginPDG = aTrackInfo.fOriginPDG;
    return *this;
}


