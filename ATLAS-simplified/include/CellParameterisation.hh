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
/// \file Geant4-models/ATLAS-simplified/include/CellParameterisation.hh
/// \brief Definition of the CaloRCellParameterisation class

#ifndef CaloRCellParameterisation_H
#define CaloRCellParameterisation_H 1

#include "Constants.hh"
#include "CellMatrix.hh"

#include "globals.hh"
#include "G4VPVParameterisation.hh"

#include <array>

class G4VPhysicalVolume;

/// EM Calorimeter cell parameterisation
/// define different class for different layers

class Em1CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Em1CellParameterisation();
    virtual ~Em1CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm1Cells> fXEm1Cell;
    std::array<G4double, kNofEm1Cells> fYEm1Cell;
};

class Em1SubCellParameterisation : public G4VPVParameterisation {
  public:
    Em1SubCellParameterisation();
    virtual ~Em1SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm1SubCells> fYEm1SubCell;
    std::array<G4double, kNofEm1SubCells> fZEm1SubCell;
};

class Em2CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Em2CellParameterisation();
    virtual ~Em2CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm2Cells> fXEm2Cell;
    std::array<G4double, kNofEm2Cells> fYEm2Cell;
};

class Em2SubCellParameterisation : public G4VPVParameterisation {
  public:
    Em2SubCellParameterisation();
    virtual ~Em2SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm2SubCells> fYEm2SubCell;
    std::array<G4double, kNofEm2SubCells> fZEm2SubCell;
};

class Em3CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Em3CellParameterisation();
    virtual ~Em3CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm3Cells> fXEm3Cell;
    std::array<G4double, kNofEm3Cells> fYEm3Cell;
};

class Em3SubCellParameterisation : public G4VPVParameterisation {
  public:
    Em3SubCellParameterisation();
    virtual ~Em3SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofEm3SubCells> fYEm3SubCell;
    std::array<G4double, kNofEm3SubCells> fZEm3SubCell;
};

class Had1CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Had1CellParameterisation();
    virtual ~Had1CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad1Cells> fXHad1Cell;
    std::array<G4double, kNofHad1Cells> fYHad1Cell;
};

class Had1SubCellParameterisation : public G4VPVParameterisation {
  public:
    Had1SubCellParameterisation();
    virtual ~Had1SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad1SubCells> fYHad1SubCell;
    std::array<G4double, kNofHad1SubCells> fZHad1SubCell;
};

class Had2CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Had2CellParameterisation();
    virtual ~Had2CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad2Cells> fXHad2Cell;
    std::array<G4double, kNofHad2Cells> fYHad2Cell;
};

class Had2SubCellParameterisation : public G4VPVParameterisation {
  public:
    Had2SubCellParameterisation();
    virtual ~Had2SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad2SubCells> fYHad2SubCell;
    std::array<G4double, kNofHad2SubCells> fZHad2SubCell;
};

class Had3CellParameterisation : public G4VPVParameterisation, public CaloRCellMatrix {
  public:
    Had3CellParameterisation();
    virtual ~Had3CellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad3Cells> fXHad3Cell;
    std::array<G4double, kNofHad3Cells> fYHad3Cell;
};

class Had3SubCellParameterisation : public G4VPVParameterisation {
  public:
    Had3SubCellParameterisation();
    virtual ~Had3SubCellParameterisation();
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  private:
    std::array<G4double, kNofHad3SubCells> fYHad3SubCell;
    std::array<G4double, kNofHad3SubCells> fZHad3SubCell;
};


#endif
