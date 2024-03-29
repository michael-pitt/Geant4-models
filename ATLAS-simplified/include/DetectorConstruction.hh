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
/// \file Geant4-models/ATLAS-simplified/include/DetectorConstruction.hh
/// \brief Definition of the CaloDetectorConstruction class
//
//

#ifndef CALOR_DETECTOR_CONSTRUCTION_H
#define CALOR_DETECTOR_CONSTRUCTION_H

#include "G4Types.hh"
#include "G4VUserDetectorConstruction.hh"
#include "Constants.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of G4MultiFunctionalDetector 
/// type with primitive scorers are created and associated with the ECAL 
/// and HCAL volumes. In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.


class CaloRDetectorConstruction : public G4VUserDetectorConstruction 
{
public:
  CaloRDetectorConstruction();
  ~CaloRDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  G4String absWord = "Abs";
  G4String gapWord = "Gap";
  G4String layerWord[fnLayers] = {"ECAL1","ECAL2","ECAL3","HCAL1","HCAL2","HCAL3"};
  G4int SensitiveToCell = 2; // depth between cell and sensitive colume (abs or gap)
  
private:
	void DefineMaterials();
	G4VPhysicalVolume* DefineVolumes();
	static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
	
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

};

#endif
