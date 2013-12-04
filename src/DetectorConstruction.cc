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
// $Id: DetectorConstruction.cc,v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "MyMaterials.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(const string& configFileName)
{
  readConfigFile(configFileName);
  
  
  
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  initializeMaterials();
  
  expHall_x = expHall_y = expHall_z = 1*m;
  
  spacingZ = abs_d + crystal_d;
  
  module_x = module_xy;
  module_y = module_xy;
  module_z = (NLAYERS_Z) * spacingZ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << ">>>>>> DetectorConstruction::Construct()::begin <<<<<<" << G4endl;
  
  
  //------------------------------------
  //------------- Volumes --------------
  //------------------------------------
  G4RotationMatrix* piRot = new G4RotationMatrix;
  piRot->rotateX(M_PI/2.*rad);  
  
  
  // The experimental Hall
  G4VSolid* worldS = new G4Box("World",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,MyMaterials::Air(),"World",0,0,0);
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"World",0,false,0,true);
  
  
  // The calorimeter
  G4VSolid* calorS = new G4Box("Calorimeter",0.5*module_xy,0.5*module_xy,0.5*module_z);
  G4LogicalVolume* calorLV = new G4LogicalVolume(calorS,MyMaterials::Air(),"Calorimeter");
  new G4PVPlacement(0,G4ThreeVector(),calorLV,"Calorimeter",worldLV,false,0,true);
  
  
  // A layer
  G4VSolid* layerS = new G4Box("Layer",0.5*module_xy,0.5*module_xy,0.5*spacingZ);
  G4LogicalVolume* layerLV = new G4LogicalVolume(layerS,MyMaterials::Air(),"Layer");
  new G4PVReplica("Layer",layerLV,calorLV,kZAxis,NLAYERS_Z,spacingZ);
  
  
  // Absorber
  G4VSolid* absorberS = new G4Box("Absorber",0.5*module_xy,0.5*module_xy,0.5*abs_d);
  G4LogicalVolume* absorberLV = new G4LogicalVolume(absorberS,AbMaterial,"Absorber");
  fAbsorberPV = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.5*spacingZ+0.5*abs_d),absorberLV,"Absorber",layerLV,false,0,true);
  
  
  // Crystal
  G4VSolid* crystalS = new G4Box("Crystal",0.5*module_xy,0.5*module_xy,0.5*crystal_d);
  G4LogicalVolume* crystalLV = new G4LogicalVolume(crystalS,ScMaterial,"Crystal");
  fCrystalPV = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.5*spacingZ+abs_d+0.5*crystal_d),crystalLV,"Crystal",layerLV,false,0,true);
  
  

  
  
  
  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------
  
  G4Colour  white   (1.0, 1.0, 1.0) ;  // white
  G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
  G4Colour  black   (0.0, 0.0, 0.0) ;  // black
  G4Colour  red     (1.0, 0.0, 0.0) ;  // red
  G4Colour  green   (0.0, 1.0, 0.0) ;  // green
  G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
  G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
  G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
  G4Colour  brass   (0.8, 0.6, 0.4) ;  // brass
  G4Colour  brown   (0.7, 0.4, 0.1) ;  // brown
  
  G4VisAttributes* VisAttWorld = new G4VisAttributes(white);
  VisAttWorld->SetVisibility(true);
  VisAttWorld->SetForceWireframe(true);
  worldLV->SetVisAttributes(VisAttWorld);
  
  G4VisAttributes* VisAttCalor = new G4VisAttributes(yellow);
  VisAttCalor->SetVisibility(false);
  VisAttCalor->SetForceWireframe(true);
  calorLV->SetVisAttributes(VisAttCalor);
  
  G4VisAttributes* VisAttLayer = new G4VisAttributes(red);
  VisAttLayer->SetVisibility(false);
  VisAttLayer->SetForceWireframe(true);
  layerLV->SetVisAttributes(VisAttLayer);
  
  G4VisAttributes* VisAttAbsorber = new G4VisAttributes(brass);
  VisAttAbsorber->SetVisibility(true);
  VisAttAbsorber->SetForceWireframe(false);
  absorberLV->SetVisAttributes(VisAttAbsorber);
  
  G4VisAttributes* VisAttCrystal = new G4VisAttributes(blue);
  VisAttCrystal->SetVisibility(true);
  VisAttCrystal->SetForceWireframe(false);
  crystalLV->SetVisAttributes(VisAttCrystal);
  
  
  
  
  G4cout << ">>>>>> DetectorConstruction::Construct()::end <<< " << G4endl;
  return worldPV;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::readConfigFile(string configFileName)
{	
  ConfigFile config(configFileName);
  
  config.readInto(crystal_material,"crystal_material");
  config.readInto(crystal_risetime,"crystal_risetime");
  config.readInto(crystal_abslength,"crystal_abslength");
  config.readInto(crystal_lightyield,"crystal_lightyield");
  
  config.readInto(fiber_radius,"fiber_radius");
  config.readInto(fiber_length,"fiber_length");
  config.readInto(module_xy,"module_xy");
  config.readInto(abs_d,"abs_d");
  config.readInto(crystal_d,"crystal_d");
  config.readInto(NLAYERS_Z,"NLAYERS_Z");
  config.readInto(abs_material,"abs_material");
  
  
  // Crystal parameters
  /*
    G4double absorber_x = config.read<double>("absorber_x")*mm;
    G4cout << "Absorber x [mm]: " << absorber_x << G4endl;
    
    G4double absorber_y = config.read<double>("absorber_y")*mm;
    G4cout << "Absorber y [mm]: " << absorber_y << G4endl;
    
    G4double absorber_z = config.read<double>("absorber_z")*mm;
    G4cout << "Absorber z [mm]: " << absorber_z << G4endl;
    
    const int NFIBERS_X = config.read<int>("NFIBERS_X");
    G4cout << "NFIBERS_X: " << NFIBERS_X << G4endl;
    const int NLAYERS_Z = config.read<int>("NLAYERS_Z");
    G4cout << "NLAYERS_Z: " << NLAYERS_Z << G4endl;
    
    G4double spacingX = config.read<double>("spacingX")*mm;
    G4cout << "spacingX [mm]: " << spacingX << G4endl;
    G4double spacingZ = config.read<double>("spacingZ")*mm;
    G4cout << "spacingZ [mm]: " << spacingZ << G4endl;
  */
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::initializeMaterials()
{
  // define materials
  
  AbMaterial = NULL;
  if( abs_material == 1 ) AbMaterial = MyMaterials::Brass();
  else if( abs_material == 2 ) AbMaterial = MyMaterials::Tungsten();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << abs_material << G4endl;
    exit(-1);
  }
  G4cout << "Ab. material: "<< AbMaterial << G4endl;
  
  ScMaterial = NULL;
  if     ( crystal_material == 1 ) ScMaterial = MyMaterials::LSO();
  else if( crystal_material == 2 ) ScMaterial = MyMaterials::LYSO();
  else if( crystal_material == 3 ) ScMaterial = MyMaterials::LuAG_Ce();  
  else if( crystal_material == 4 ) ScMaterial = MyMaterials::LuAG_Pr();
  else if( crystal_material == 5 ) ScMaterial = MyMaterials::PWO();
  else if( crystal_material == 6 ) ScMaterial = MyMaterials::Air();
  else if( crystal_material == 7 ) ScMaterial = MyMaterials::Quartz();
  else if( crystal_material == 8 ) ScMaterial = MyMaterials::DSBCe();
  else if( crystal_material == 9 ) ScMaterial = MyMaterials::SiO2Ce();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid crystal material specifier " << crystal_material << G4endl;
    exit(-1);
  }
  G4cout << "Sc. material: "<< ScMaterial << G4endl;
  
  
  
  
  
  // modify default properties of the scintillator
  if( crystal_lightyield >= 0 )
  {
    ScMaterial->GetMaterialPropertiesTable()->RemoveConstProperty("SCINTILLATIONYIELD");
    ScMaterial->GetMaterialPropertiesTable()->AddConstProperty("SCINTILLATIONYIELD",crystal_lightyield/MeV);  
  } 
  
  if( crystal_risetime >= 0 )
  {
    ScMaterial->GetMaterialPropertiesTable()->RemoveConstProperty("FASTSCINTILLATIONRISETIME");
    ScMaterial->GetMaterialPropertiesTable()->AddConstProperty("FASTSCINTILLATIONRISETIME",crystal_risetime/ns);  
  } 
  
  if( crystal_abslength >= 0 ) 
  {
    ScMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    ScMaterial->GetMaterialPropertiesTable()->AddConstProperty("ABSLENGTH",crystal_abslength);  
  } 
  else 
  {
    for(unsigned int j = 0; j < ScMaterial->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH")->GetVectorLength();++j)
    {
      ScMaterial->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH")->Energy(j);
    }
  }
}
