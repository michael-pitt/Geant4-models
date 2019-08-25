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

// TODO, include how to handle with row_wise_branch (for MT jobs)
// https://gitlab.cern.ch/geant4/geant4/blob/master/examples/extended/medical/dna/microyz/plot.C
//#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif
#include "G4UImanager.hh"

#include "QGSP_FTFP_BERT.hh"
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

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " CaloR [-m macro ] [-u UIsession] [-s seedNumber] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

int main(int argc, char** argv)
{
	
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4long seed = (long) time(NULL);
  G4String macro;
  G4String session;
  #ifdef G4MULTITHREADED
  G4int nThreads = 0;
  #endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
	else if ( G4String(argv[i]) == "-s" ) {
      seed = G4UIcommand::ConvertToInt(argv[i+1]);
	  G4cout << "Using user random seed = " <<  seed  << G4endl;
    }	
    #ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = std::min(G4Threading::G4GetNumberOfCores(),G4UIcommand::ConvertToInt(argv[i+1]));
	  G4cout << "with G4MULTITHREADED: using "<<nThreads<<" of "
			 <<G4Threading::G4GetNumberOfCores()<<" threads"<<G4endl;
    }
    #endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }
  
  //choose the Random engine and set a random seed
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << "Seed: " << CLHEP::HepRandom::getTheSeed() << G4endl;

//  #ifdef G4MULTITHREADED
//  G4MTRunManager * runManager = new G4MTRunManager;
//  if ( nThreads > 0 ) { 
//    runManager->SetNumberOfThreads(nThreads);
//  }  
//  #else
  G4RunManager * runManager = new G4RunManager;
//  #endif
  
  // Set mandatory initialization classes
  
  // Initialize detector geometry
  auto detector = new CaloRDetectorConstruction;
  runManager-> SetUserInitialization(detector);

  
  // Use GSP + FTF/Preco + BERT:
  // Gluon String Plasma model (QGSP) for high energies the hadron showers
  // Fritiof string model (FTF) for  hadron-nucleus at Plab interactions >3GeV
  // Precompound and deexcitation (Preco)
  // Bertini Cascade (BERT) for  hadron-nucleus at Plab interactions <3GeV
  //G4VModularPhysicsList* physics = new QGSP_FTFP_BERT;
  //physics->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4VUserPhysicsList* physics = new QGSP_FTFP_BERT;
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
    UImanager-> ApplyCommand(command+macro);    
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

