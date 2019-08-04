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
/// \file Geant4-models/ATLAS-simplified/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRPrimaryGeneratorAction::CaloRPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr)
{
	
  fParticleGun = new G4ParticleGun(1);
  // default particle kinematic (unless set in /gun/particle ...)
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("e+");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy((2 + 98 * G4UniformRand() ) * GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., Zinit));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CaloRPrimaryGeneratorAction::~CaloRPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CaloRPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  // If you whish to code the generated particles
  // comment this line and uncomment the next block
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  
  
  // TwoPions
  G4String ListParticles[13] = {"mu+", "mu-", "pi+", "pi-",
							 "e+", "e-", "proton", "anti_proton",
							 "gamma", "neutron", "anti_neutron","anti_kaon0","kaon0L"};
					
  G4double energy, dx0, dy0, dx_calo, dy_calo;
  int i_particle = int(12.999 * G4UniformRand() );
  i_particle = 0;
  
  // generate pi0
  fParticleGun = nullptr ;
  fParticleGun = new G4ParticleGun(1);
  dx_calo = -dR_01 + 2*(dR_01)*G4UniformRand(); dy_calo = -dR_01 + 2*(dR_01)*G4UniformRand();
  dx0 = dx_calo - 4*dR_01 + 8*(dR_01)*G4UniformRand();
  dy0 = dy_calo - 4*dR_01 + 8*(dR_01)*G4UniformRand();
  energy = (10 + 20 * G4UniformRand() ) * GeV;

  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ListParticles[i_particle]));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx0-dx_calo, dy0-dy_calo, GunDinsance));
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(G4ThreeVector(-dx0, -dy0, Zinit));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // generate pi+
  /*
  fParticleGun = nullptr ;
  fParticleGun = new G4ParticleGun(1);
  theta = CLHEP::pi/20. + CLHEP::pi/10. *G4UniformRand();
  phi = phi - CLHEP::pi/20. + CLHEP::pi/20. *G4UniformRand();
  energy = (15 + 15 * G4UniformRand() ) * GeV;

  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi+"));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector( sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)));
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., Zinit));
  fParticleGun->GeneratePrimaryVertex(anEvent);
  */
  
}
