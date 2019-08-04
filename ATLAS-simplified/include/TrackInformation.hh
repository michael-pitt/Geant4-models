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
/// \file Geant4-models/ATLAS-simplified/include/TrackInformation.hh
/// \brief Definition of the CaloRTrackInformation class
//
//
//

#ifndef CaloRTrackInformation_h
#define CaloRTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class CaloRTrackInformation : public G4VUserTrackInformation 
{
public:
  CaloRTrackInformation();
  CaloRTrackInformation(const G4Track* aTrack);
  CaloRTrackInformation(const CaloRTrackInformation* aTrackInfo);
  virtual ~CaloRTrackInformation();
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);

  CaloRTrackInformation& operator =(const CaloRTrackInformation& right);
  

public:
//  inline G4int GetTrackingStatus() const {return fTrackingStatus;}
//  inline void  SetTrackingStatus(G4int i) {fTrackingStatus = i;}
//  inline G4int GetSourceTrackID() const {return fSourceTrackID;}
//  inline void  SetSuspendedStepID(G4int i) {fSuspendedStepID = i;}
  inline G4double GetOriginCharge() const {return fOriginCharge;}
  inline G4int GetOriginPDG() const {return fOriginPDG;}
  inline void SetOriginCharge(G4double charge) {fOriginCharge = charge;}
  inline void SetOriginPDG(G4int pdgId) {fOriginPDG = pdgId;}

private:
  // Information of the primary track
  G4double              fOriginCharge;  // charge of primary particle
  G4int                 fOriginPDG;  // PdgID of primary particle
};

extern G4ThreadLocal
 G4Allocator<CaloRTrackInformation> * aTrackInformationAllocator;

inline void* CaloRTrackInformation::operator new(size_t)
{
  if(!aTrackInformationAllocator)
    aTrackInformationAllocator = new G4Allocator<CaloRTrackInformation>;
  return (void*)aTrackInformationAllocator->MallocSingle();
}

inline void CaloRTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator->FreeSingle((CaloRTrackInformation*)aTrackInfo);}

#endif

