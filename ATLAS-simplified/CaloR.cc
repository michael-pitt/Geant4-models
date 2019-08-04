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
// Basic tutorial: https://www.ge.infn.it/geant4/training/ornl_2008/retrievinginformationfromkernel.pdf
/// \file CaloR.cc
/// \brief Main program of the CaloR example
//
//
#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"

#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
//#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronicProcessStore.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

// package includes
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"


#include "TROOT.h"

int main(int argc, char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  G4UIExecutive* ui = nullptr;
  G4long seed(0);
  if ( argc == 1 ) {
	// to get rid of the nasty "LLVM SYMBOLS ARE EXPOSED TO CLING" error
	// https://root-forum.cern.ch/t/error-llvm-symbols-exposed-to-cling/23597
	gROOT->GetInterpreter(); 
    ui = new G4UIExecutive(argc, argv);
  }
  else if( argc == 2 ){
	  G4cout << "Using \"time(NULL)\" random seed" << G4endl;
	  seed = (long) time(NULL);
  }
  else if ( argc == 3 ) {
	  seed = (long) atoi(argv[2]);
	  G4cout << "Using user random seed = " <<  seed  << G4endl;
  }
  else {
    G4cout << "ERROR: wrong program usage" << G4endl;
    return 0;
  }

	  

  //choose the Random engine and set a random seed
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << "Seed: " << CLHEP::HepRandom::getTheSeed() << G4endl;

  #ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4);
  #else
  G4RunManager * runManager = new G4RunManager;
  #endif
  
  // Set mandatory initialization classes
  
  // Initialize detector geometry
  auto detector = new CaloRDetectorConstruction;
  runManager-> SetUserInitialization(detector);

  
  // Use FTF/Preco and BERT:
  // Fritiof string model (FTF) for  hadron-nucleus at Plab interactions >3GeV
  // Precompound and deexcitation (Preco)
  // Bertini Cascade (BERT) for  hadron-nucleus at Plab interactions <3GeV
  //G4VModularPhysicsList* physics = new FTFP_BERT;
  //physics->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4VUserPhysicsList* physics = new FTFP_BERT;
  runManager-> SetUserInitialization(physics);
  
  // Supress annoying "HADRONIC PROCESSES SUMMARY"
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  
  //Set Pi0 decay to photon pair
  G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* fParticleDef = fParticleTable->FindParticle("pi0");
  G4VDecayChannel* fMode =
                  new G4PhaseSpaceDecayChannel("pi0",1,2,"gamma","gamma");
  G4DecayTable* fTable = new G4DecayTable();
  fTable->Insert(fMode);
  fParticleDef->SetDecayTable(fTable);
  
  // User Action classes
  auto actionInitialization = new CaloRActionInitialization(detector);
  runManager->SetUserInitialization(actionInitialization);

  // Initialize for ParticleGun settings
  runManager-> Initialize();
  
  // Initialize visualization
  auto visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;

 //get the pointer to the User Interface manager   
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (!ui) { // batch mode
    visManager-> SetVerboseLevel(0);
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager-> ApplyCommand(command+fileName);    
  } else {  // interactive mode : define UI session

    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }

    ui-> SessionStart();
    delete ui;
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;
  delete runManager;
}

