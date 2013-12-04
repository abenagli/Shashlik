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
// $Id: DetectorConstruction.hh,v 1.5 2006-06-29 17:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "ConfigFile.hh"
#include "TString.h"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  DetectorConstruction(const string& configFileName);
  ~DetectorConstruction();
  
  G4double GetModule_x() const { return module_x; };
  G4double GetModule_y() const { return module_y; };
  G4double GetModule_z() const { return module_z; };
  
public:
  G4VPhysicalVolume* Construct();
  
private:
  G4VPhysicalVolume* fAbsorberPV;   // the absorber physical volume
  G4VPhysicalVolume* fCrystalPV;   // the absorber physical volume
  
  
  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;
  
  G4int    crystal_material;
  G4int    crystal_lightyield;
  G4double crystal_risetime;
  G4double crystal_abslength;
  G4double crystal_induced_abslength;
  
  G4double fiber_length;
  G4double fiber_radius;
  int NMODULES_X, NMODULES_Y;
  int NFIBERS_X, NLAYERS_Z;
  G4double spacingX, spacingZ;
  
  G4double module_xy;
  G4double module_x;
  G4double module_y;
  G4double module_z;
  G4double brass_hole_radius;
  G4int abs_material;
  
  G4double abs_d;
  G4double crystal_d;
  
  G4double win_r;
  G4double win_l;
  G4int    win_material;
  
  G4double det_d;
  G4double det_distance;
  G4int    det_material;
  
  G4double depth;
  
  void readConfigFile(string configFileName);
  
  //Materials
  void initializeMaterials();
  G4Material* AbMaterial;
  G4Material* GaMaterial;
  G4Material* ScMaterial;
  G4Material* WiMaterial;
  G4Material* DeMaterial;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
