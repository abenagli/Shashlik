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
// * institutes, nor the agencies providing financial support for this *
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
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  readConfigFile (configFileName) ;
  
  
  
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  initializeMaterials () ;
  
  expHall_x = expHall_y = expHall_z = 1*m ;
  
  spacing_z = abs_d + crystal_d ;
  
  module_x = module_xy ;
  module_y = module_xy ;
  module_z = (nLayers_z) * spacing_z ;
  
  fiber_length = module_z ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction ()
{}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;
  
  
  
  //------------------------------------
  //------------- Volumes --------------
  //------------------------------------
  
  
  std::vector<G4TwoVector> crystalBase ;
  fillPolygon (crystalBase, 0.5*module_xy, chamfer) ;
  
  
  // The experimental Hall
  G4VSolid* worldS = new G4Box ("World", 0.5*expHall_x, 0.5*expHall_y, 0.5*expHall_z) ;
  G4LogicalVolume* worldLV = new G4LogicalVolume (worldS, MyMaterials::Air (), "World", 0, 0, 0) ;
  G4VPhysicalVolume* worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "World", 0, false, 0, true) ;
  
  
  // The calorimeter
  G4VSolid* calorS = new G4Box ("Calorimeter", 0.5*module_xy, 0.5*module_xy, 0.5*module_z) ;
  G4LogicalVolume* calorLV = new G4LogicalVolume (calorS, MyMaterials::Air (), "Calorimeter") ;
  new G4PVPlacement (0, G4ThreeVector (), calorLV, "Calorimeter", worldLV, false, 0, true) ;
  
  
  // A layer
  G4VSolid* layerS = new G4Box ("Layer", 0.5*module_xy, 0.5*module_xy, 0.5*spacing_z) ;
  G4LogicalVolume* layerLV = new G4LogicalVolume (layerS, MyMaterials::Air (), "Layer") ;
  new G4PVReplica ("Layer", layerLV, calorLV, kZAxis, nLayers_z, spacing_z) ;
  
  
  // Crystal
  G4VSolid* crystalS = new G4ExtrudedSolid ("Crystal", crystalBase, 0.5*crystal_d, G4TwoVector (0., 0.), 1., G4TwoVector (0., 0.), 1.) ;
  G4LogicalVolume* crystalLV = new G4LogicalVolume (crystalS, ScMaterial, "Crystal") ;
  fCrystalPV = new G4PVPlacement (0, G4ThreeVector (0., 0., -0.5*spacing_z+0.5*crystal_d), crystalLV, "Crystal", layerLV, false, 0, true) ;
    
  // Absorber
  G4VSolid* absorberS = new G4ExtrudedSolid ("Absorber", crystalBase, 0.5*abs_d, G4TwoVector (0., 0.), 1., G4TwoVector (0., 0.), 1.) ;
  G4LogicalVolume* absorberLV = new G4LogicalVolume (absorberS, AbMaterial, "Absorber") ;
  fAbsorberPV = new G4PVPlacement (0, G4ThreeVector (0., 0., -0.5*spacing_z+0.5*abs_d+crystal_d), absorberLV, "Absorber", layerLV, false, 0, true) ;
  
  // Fibers
  //      The core is divided into two sub-cores, made of the same material.
  //      This does not change at all the physics, but eases the check for where the photon is passing through.
 
  G4VSolid* fiberCoreInsS = new G4Tubs ("FiberCoreIns", 0.                       , fiberCore_radius * 0.9999, 0.5*fiber_length, 0.*deg, 360.*deg) ;
  G4VSolid* fiberCoreOutS = new G4Tubs ("FiberCoreOut", fiberCore_radius * 0.9999, fiberCore_radius         , 0.5*fiber_length, 0.*deg, 360.*deg) ;
  G4VSolid* fiberCladS    = new G4Tubs ("FiberClad"   , fiberCore_radius         , fiberClad_radius         , 0.5*fiber_length, 0.*deg, 360.*deg) ;
  

  G4double offset_x = 0. ;
  G4double offset_y = 0. ;

  //PG first edge: chessboard disposition
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  G4LogicalVolume* fiberCoreInsLV_0 = new G4LogicalVolume (fiberCoreInsS, CoMaterial, "FiberCoreIns_0") ;
  G4LogicalVolume* fiberCoreOutLV_0 = new G4LogicalVolume (fiberCoreOutS, CoMaterial, "FiberCoreOut_0") ;
  G4LogicalVolume* fiberCladLV_0    = new G4LogicalVolume (fiberCladS,    ClMaterial, "FiberClad_0") ;

  fFiberCoreInsPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCoreOutPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCladPV.push_back    (std::vector <G4VPhysicalVolume*> ()) ;

  int edge = 0 ;
  std::pair<G4TwoVector, G4TwoVector> theChamfer = getChamfer (crystalBase, edge) ;

  int nTotFibers = 0 ;  
  int numberOfRadius = 1 ;
  while (1)
  {
    int fibersNumberInFirstRow = floor ( 
      ( (theChamfer.first - theChamfer.second).mag () - 2.* numberOfRadius * fiberClad_radius ) / // length of the usable line
      (2 * fiberClad_radius)                                                                     // transverse length of a single fiber
      ) ;
    if ( fibersNumberInFirstRow <= 0 ) break ;
    
    G4double offset_x = 0. ;
    G4double offset_y = 0. ;
    
    //PG put the first fiber
    G4TwoVector fiberAxisPosition = centerOfTheFirstFiber (theChamfer, fibersNumberInFirstRow, fiberClad_radius, numberOfRadius) ;
    fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_0, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
    fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_0, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
    fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_0, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
    ++nTotFibers ;
    
    //PG add the following fibers in the line
    for (int i = 1 ; i < fibersNumberInFirstRow ; ++i)
    {
      fiberAxisPosition = getNextCenter (theChamfer, fiberAxisPosition, fiberClad_radius) ;
      fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_0, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
      fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_0, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
      fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_0, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;  
      ++nTotFibers ;
    }
    
    numberOfRadius += 2 ;
  } // while

  //PG second edge: a single line
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  G4LogicalVolume* fiberCoreInsLV_1 = new G4LogicalVolume (fiberCoreInsS, CoMaterial, "FiberCoreIns_1") ;
  G4LogicalVolume* fiberCoreOutLV_1 = new G4LogicalVolume (fiberCoreOutS, CoMaterial, "FiberCoreOut_1") ;
  G4LogicalVolume* fiberCladLV_1    = new G4LogicalVolume (fiberCladS,    ClMaterial, "FiberClad_1"   ) ;

  fFiberCoreInsPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCoreOutPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCladPV.push_back    (std::vector <G4VPhysicalVolume*> ()) ;

  edge = 1 ;
  theChamfer = getChamfer (crystalBase, edge) ;

  numberOfRadius = 1 ;
  int fibersNumberInFirstRow = floor ( 
    ( (theChamfer.first - theChamfer.second).mag () - 2.* numberOfRadius * fiberClad_radius ) / // length of the usable line
    (2 * fiberClad_radius)                                                                     // transverse length of a single fiber
    ) ;

  //PG put the first fiber
  G4TwoVector fiberAxisPosition = centerOfTheFirstFiber (theChamfer, fibersNumberInFirstRow, fiberClad_radius, numberOfRadius) ;
  fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_1, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
  fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_1, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
  fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_1, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
  
  //PG add the following fibers in the line
  for (int i = 1 ; i < fibersNumberInFirstRow ; ++i)
  {
    fiberAxisPosition = getNextCenter (theChamfer, fiberAxisPosition, fiberClad_radius) ;
    fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_1, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
    fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_1, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
    fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_1, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
  }
  
  //PG third edge: the most compact disposition is possible
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  G4LogicalVolume* fiberCoreInsLV_2 = new G4LogicalVolume (fiberCoreInsS, CoMaterial, "FiberCoreIns_2") ;
  G4LogicalVolume* fiberCoreOutLV_2 = new G4LogicalVolume (fiberCoreOutS, CoMaterial, "FiberCoreOut_2") ;
  G4LogicalVolume* fiberCladLV_2    = new G4LogicalVolume (fiberCladS,    ClMaterial, "FiberClad_2"   ) ;

  fFiberCoreInsPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCoreOutPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCladPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;

  edge = 2 ;
  theChamfer = getChamfer (crystalBase, edge) ;

  fibersNumberInFirstRow = floor (
      (    (theChamfer.first - theChamfer.second).mag () - 2 * fiberClad_radius * 1.41421356237)   // length of the usable line
        /  (2 * fiberClad_radius)   // length of a single fibre
    ) ;

  numberOfRadius = 1 ;
  //PG put the first fiber
  fiberAxisPosition = centerOfTheFirstFiberPG (theChamfer, fibersNumberInFirstRow, fiberClad_radius) ;
  G4TwoVector firstFiberInRowCenter = fiberAxisPosition ;
  fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_2, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
  fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_2, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
  fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_2, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
  
  //PG add the following fibers in the first line
  for (int i = 1 ; i < fibersNumberInFirstRow ; ++i)
    {
      fiberAxisPosition = getNextCenter (theChamfer, fiberAxisPosition, fiberClad_radius) ;
      fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_2, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
      fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_2, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
      fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_2, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
    }

  G4TwoVector chamferDirection = theChamfer.second - theChamfer.first ;
  G4TwoVector chamferOrtogonal = chamferDirection ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  chamferOrtogonal.setX (chamferDirection.y ()) ;
  chamferOrtogonal.setY (-1 * chamferDirection.x ()) ;

  int safetyCounter = 0 ;
  do {
      ++safetyCounter ;
      // put the first fibre of the second row
      fiberAxisPosition = centerOfTheFirstFibreOnSecondLayer (theChamfer, fiberClad_radius, firstFiberInRowCenter) ;
      firstFiberInRowCenter = fiberAxisPosition ;

      if (chamferOrtogonal.dot (fiberAxisPosition - theChamfer.second) > 0.5 * chamfer)
        { 
          //PG no more lines of fibres can be put FIXME
          break ;
        }

      if (checkIfOutOfChamfer (fiberClad_radius, fiberAxisPosition, crystalBase, 2)) 
        {
          fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_2, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
          fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_2, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
          fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_2, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
        }

      //PG add the following fibres in the line
      for (int i = 1 ; i < fibersNumberInFirstRow ; ++i)
        {
          fiberAxisPosition = getNextCenter (theChamfer, fiberAxisPosition, fiberClad_radius) ;
          if (checkIfOutOfChamfer (fiberClad_radius, fiberAxisPosition, crystalBase, 2)) 
            {
              fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreInsLV_2, Form ("FiberCoreIns%d", edge), worldLV, false, 0, false)) ;
              fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCoreOutLV_2, Form ("FiberCoreOut%d", edge), worldLV, false, 0, false)) ;
              fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), fiberCladLV_2, Form ("FiberClad%d", edge), worldLV, false, 0, false)) ;
            }
          else
            { continue ; }  
        }
     }  
   while (safetyCounter < 20) ;
   if (safetyCounter > 20)
     {
       std::cout << "WARNING: abnormal termination of loop in filling the third chamfer" << std::endl ;
     }

  //PG fourth edge: a single large fiber
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  float bigfiberClad_radius = chamfer * 1.73205080757 * 0.125 ;
  bigfiberClad_radius *= 0.95 ;
  float bigfiberCore_radius = 0.5 * bigfiberClad_radius ;

  G4VSolid* bigfiberCoreInsS = new G4Tubs ("bigfiberCoreIns", 0., bigfiberCore_radius * 0.9999, 0.5*fiber_length, 0.*deg, 360.*deg) ;
  G4VSolid* bigfiberCoreOutS = new G4Tubs ("bigfiberCoreOut", bigfiberCore_radius * 0.9999, bigfiberCore_radius, 0.5*fiber_length, 0.*deg, 360.*deg) ;
  G4VSolid* bigfiberCladS = new G4Tubs ("bigfiberClad", bigfiberCore_radius, bigfiberClad_radius, 0.5*fiber_length, 0.*deg, 360.*deg) ;
  
  G4LogicalVolume* bigfiberCoreInsLV = new G4LogicalVolume (bigfiberCoreInsS, CoMaterial, "fiberCoreIns_3") ;
  G4LogicalVolume* bigfiberCoreOutLV = new G4LogicalVolume (bigfiberCoreOutS, CoMaterial, "fiberCoreOut_3") ;
  G4LogicalVolume* bigfiberCladLV = new G4LogicalVolume (bigfiberCladS, ClMaterial, "fiberClad_3") ;

  fFiberCoreInsPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCoreOutPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;
  fFiberCladPV.push_back (std::vector <G4VPhysicalVolume*> ()) ;

  // find the center
  
  edge = 3 ;
  theChamfer = getChamfer (crystalBase, edge) ;

  chamferDirection = theChamfer.second - theChamfer.first ;
  chamferOrtogonal = chamferDirection ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  chamferOrtogonal.setX (chamferDirection.y ()) ;
  chamferOrtogonal.setY (-1 * chamferDirection.x ()) ;
  
  fiberAxisPosition = theChamfer.first 
      + 0.5 * chamfer * chamferDirection
      + bigfiberClad_radius * chamferOrtogonal ;
  fFiberCoreInsPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), bigfiberCoreInsLV, Form ("BigFiberCoreIns%d", edge), worldLV, false, 0, false)) ;
  fFiberCoreOutPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), bigfiberCoreOutLV, Form ("BigFiberCoreOut%d", edge), worldLV, false, 0, false)) ;
  fFiberCladPV.back ().push_back (new G4PVPlacement (0, G4ThreeVector (offset_x+fiberAxisPosition.x (), offset_y+fiberAxisPosition.y (), 0.), bigfiberCladLV, Form ("BigFiberClad%d", edge), worldLV, false, 0, false)) ;
  
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
  
  G4VisAttributes* VisAttWorld = new G4VisAttributes (white) ;
  VisAttWorld->SetVisibility (true) ;
  VisAttWorld->SetForceWireframe (true) ;
  worldLV->SetVisAttributes (VisAttWorld) ;
  
  G4VisAttributes* VisAttCalor = new G4VisAttributes (yellow) ;
  VisAttCalor->SetVisibility (true) ;
  VisAttCalor->SetForceWireframe (true) ;
  calorLV->SetVisAttributes (VisAttCalor) ;
  
  G4VisAttributes* VisAttLayer = new G4VisAttributes (red) ;
  VisAttLayer->SetVisibility (false) ;
  VisAttLayer->SetForceWireframe (true) ;
  layerLV->SetVisAttributes (VisAttLayer) ;
  
  G4VisAttributes* VisAttAbsorber = new G4VisAttributes (gray) ;
  VisAttAbsorber->SetVisibility (true) ;
  VisAttAbsorber->SetForceWireframe (false) ;
  absorberLV->SetVisAttributes (VisAttAbsorber) ;
  
  G4VisAttributes* VisAttCrystal = new G4VisAttributes (blue) ;
  VisAttCrystal->SetVisibility (true) ;
  VisAttCrystal->SetForceWireframe (false) ;
  crystalLV->SetVisAttributes (VisAttCrystal) ;

  G4VisAttributes* VisAttFiberCore = new G4VisAttributes (green) ;
  VisAttFiberCore->SetVisibility (true) ;
  VisAttFiberCore->SetForceWireframe (false) ;
  fiberCoreInsLV_0->SetVisAttributes (VisAttFiberCore) ;  
  fiberCoreOutLV_0->SetVisAttributes (VisAttFiberCore) ;  
  fiberCoreInsLV_1->SetVisAttributes (VisAttFiberCore) ;  
  fiberCoreOutLV_1->SetVisAttributes (VisAttFiberCore) ;  
  fiberCoreInsLV_2->SetVisAttributes (VisAttFiberCore) ;  
  fiberCoreOutLV_2->SetVisAttributes (VisAttFiberCore) ;  
  bigfiberCoreInsLV->SetVisAttributes (VisAttFiberCore) ;  
  bigfiberCoreOutLV->SetVisAttributes (VisAttFiberCore) ;  
  
  G4VisAttributes* VisAttFiberClad = new G4VisAttributes (cyan) ;
  VisAttFiberClad->SetVisibility (true) ;
  VisAttFiberClad->SetForceWireframe (false) ;
  fiberCladLV_0->SetVisAttributes (VisAttFiberClad) ;  
  fiberCladLV_1->SetVisAttributes (VisAttFiberClad) ;  
  fiberCladLV_2->SetVisAttributes (VisAttFiberClad) ;  
  bigfiberCladLV->SetVisAttributes (VisAttFiberClad) ;  
  
  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::readConfigFile (string configFileName)
{	
  ConfigFile config (configFileName) ;
  
  config.readInto (chamfer, "chamfer") ;
  config.readInto (module_xy, "module_xy") ;
  config.readInto (nLayers_z, "nLayers_z") ;
  
  config.readInto (abs_material, "abs_material") ;
  config.readInto (abs_d, "abs_d") ;
  
  config.readInto (crystal_material, "crystal_material") ;
  config.readInto (crystal_risetime, "crystal_risetime") ;
  config.readInto (crystal_abslength, "crystal_abslength") ;
  config.readInto (crystal_lightyield, "crystal_lightyield") ;
  config.readInto (crystal_d, "crystal_d") ;
  
  config.readInto (fiberCore_material, "fiberCore_material") ;
  config.readInto (fiberCore_radius, "fiberCore_radius") ;
  config.readInto (fiberClad_material, "fiberClad_material") ;
  config.readInto (fiberClad_radius, "fiberClad_radius") ;
  config.readInto (fiber_length, "fiber_length") ;
  
  config.readInto (depth, "depth") ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::initializeMaterials ()
{
  // define materials
  
  AbMaterial = NULL ;
  if     ( abs_material == 1 ) AbMaterial = MyMaterials::Brass () ;
  else if ( abs_material == 2 ) AbMaterial = MyMaterials::Tungsten () ;
  else if ( abs_material == 3 ) AbMaterial = MyMaterials::Lead () ;
  else if ( abs_material == 4 ) AbMaterial = MyMaterials::Iron () ;
  else if ( abs_material == 5 ) AbMaterial = MyMaterials::Aluminium () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << abs_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Ab. material: "<< AbMaterial << G4endl ;
  
  
  ScMaterial = NULL ;
  if      ( crystal_material == 1 ) ScMaterial = MyMaterials::LSO () ;
  else if ( crystal_material == 2 ) ScMaterial = MyMaterials::LYSO () ;
  else if ( crystal_material == 3 ) ScMaterial = MyMaterials::LuAG_Ce () ;  
  else if ( crystal_material == 4 ) ScMaterial = MyMaterials::LuAG_Pr () ;
  else if ( crystal_material == 5 ) ScMaterial = MyMaterials::PWO () ;
  else if ( crystal_material == 6 ) ScMaterial = MyMaterials::Air () ;
  else if ( crystal_material == 7 ) ScMaterial = MyMaterials::Quartz () ;
  else if ( crystal_material == 8 ) ScMaterial = MyMaterials::DSB_Ce () ;
  else if ( crystal_material == 9 ) ScMaterial = MyMaterials::SiO2_Ce () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid crystal material specifier " << crystal_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Sc. material: "<< ScMaterial << G4endl ;
  
  
  CoMaterial = NULL ;
  if      ( fiberCore_material == 1 ) CoMaterial = MyMaterials::Quartz () ;
  else if ( fiberCore_material == 2 ) CoMaterial = MyMaterials::SiO2_Ce () ;
  else if ( fiberCore_material == 3 ) CoMaterial = MyMaterials::DSB_Ce () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fiber core material specifier " << fiberCore_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Co. material: "<< CoMaterial << G4endl ;
  
  
  ClMaterial = NULL ;
  if      ( fiberClad_material == 1 ) ClMaterial = MyMaterials::Quartz () ;
  else if ( fiberClad_material == 2 ) ClMaterial = MyMaterials::SiO2_Ce () ;
  else if ( fiberClad_material == 3 ) ClMaterial = MyMaterials::DSB_Ce () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fiber clad material specifier " << fiberClad_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cl. material: "<< ClMaterial << G4endl ;
  
  
  
  // modify default properties of the scintillator
  if ( crystal_lightyield >= 0 )
  {
    ScMaterial->GetMaterialPropertiesTable ()->RemoveConstProperty ("SCINTILLATIONYIELD") ;
    ScMaterial->GetMaterialPropertiesTable ()->AddConstProperty ("SCINTILLATIONYIELD", crystal_lightyield/MeV) ;  
  } 
  
  if ( crystal_risetime >= 0 )
  {
    ScMaterial->GetMaterialPropertiesTable ()->RemoveConstProperty ("FASTSCINTILLATIONRISETIME") ;
    ScMaterial->GetMaterialPropertiesTable ()->AddConstProperty ("FASTSCINTILLATIONRISETIME", crystal_risetime/ns) ;  
  } 
  
  if ( crystal_abslength >= 0 ) 
  {
    ScMaterial->GetMaterialPropertiesTable ()->RemoveProperty ("ABSLENGTH") ;
    ScMaterial->GetMaterialPropertiesTable ()->AddConstProperty ("ABSLENGTH", crystal_abslength) ;  
  } 
  else 
  {
    for (unsigned int j = 0 ; j < ScMaterial->GetMaterialPropertiesTable ()->GetProperty ("ABSLENGTH")->GetVectorLength () ;++j)
    {
      ScMaterial->GetMaterialPropertiesTable ()->GetProperty ("ABSLENGTH")->Energy (j) ;
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::fillPolygon (std::vector<G4TwoVector>& theBase, const float& side, const float& chamfer)
{
  // wrt the centre (0, 0)
  double delta = side - chamfer * 0.707106781188 ;
  theBase.push_back (G4TwoVector (side, delta)) ;
  theBase.push_back (G4TwoVector (delta, side)) ;
  theBase.push_back (G4TwoVector (-1 * delta, side)) ;
  theBase.push_back (G4TwoVector (-1 * side, delta)) ;
  theBase.push_back (G4TwoVector (-1 * side, -1 * delta)) ;
  theBase.push_back (G4TwoVector (-1 * delta, -1 * side)) ;
  theBase.push_back (G4TwoVector (delta, -1 * side)) ;
  theBase.push_back (G4TwoVector (side, -1 * delta)) ;
  return ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


std::pair<G4TwoVector, G4TwoVector> DetectorConstruction::getChamfer (std::vector<G4TwoVector>& theBase, const int& index)
{
  return std::pair<G4TwoVector, G4TwoVector> (
    theBase.at (2 * index), theBase.at (2*index + 1) ) ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4TwoVector DetectorConstruction::centerOfTheFirstFiber (std::pair<G4TwoVector, G4TwoVector>& theChamfer, const int& fibersNumberInRow, const float& fiberExternalRadius, const int& numberOfRadius)
{
  // assume that the chamfer coordinates are given counter-clockwise
  // so the orthogonal vector aims towards the exterior of the chamfer
  G4TwoVector chamferDirection = theChamfer.second - theChamfer.first ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  G4TwoVector chamferOrtogonal = chamferDirection ;
  chamferOrtogonal.setX (chamferDirection.y ()) ;
  chamferOrtogonal.setY (-1 * chamferDirection.x ()) ;
  G4double freeSpace = (theChamfer.second-theChamfer.first).mag () - 2.*fiberExternalRadius - 2.*fiberExternalRadius*numberOfRadius - 2.*fiberExternalRadius* (fibersNumberInRow-2) ;
  
  return theChamfer.first +
    fiberExternalRadius * numberOfRadius * chamferOrtogonal +  // go out for the length of the radius
    fiberExternalRadius * numberOfRadius * chamferDirection +  // move towards the other edge for the space where fibers cannot fit
    0.5 * freeSpace * chamferDirection ;                        // equally divide free space
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4TwoVector 
DetectorConstruction::centerOfTheFirstFiberPG (
    std::pair<G4TwoVector, G4TwoVector>&  theChamfer, 
    const int&                           fibersNumberInRow, 
    const float&                         fiberExternalRadius)
{
  // assume that the chamfer coordinates are given counter-clockwise
  // so the orthogonal vector aims towards the exterior of the chamfer
  G4TwoVector chamferDirection = theChamfer.second - theChamfer.first ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  G4TwoVector chamferOrtogonal = chamferDirection ;
  chamferOrtogonal.setX (chamferDirection.y ()) ;
  chamferOrtogonal.setY (-1 * chamferDirection.x ()) ;
  G4double freeSpace = (theChamfer.second-theChamfer.first).mag () 
                      - 2. * fiberExternalRadius * fibersNumberInRow ;
  return theChamfer.first +                                        // the starting point
         fiberExternalRadius * chamferOrtogonal +                  // go out for the length of the radius
         (0.5 * freeSpace + fiberExternalRadius) * chamferDirection ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4TwoVector DetectorConstruction::getNextCenter (std::pair<G4TwoVector, G4TwoVector>& theChamfer, G4TwoVector& thisCenter, const float& fiberExternalRadius)
{
  G4TwoVector chamferDirection = theChamfer.second - theChamfer.first ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  return thisCenter + (2 * fiberExternalRadius) * chamferDirection ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


/**
- check whether a fibre is out of the chamfer
- does NOT check whether the fibre enters the crystal volume
- assumes that the edges follow this notation in the solid front face:
      6     5
   7 /-----\ 4
     |     |
   0 \_____/ 3
      1    2
- the chamfer index goes as:      
   3 /-----\ 2
     |     |
   0 \_____/ 1      
*/
bool 
DetectorConstruction::checkIfOutOfChamfer (
    double                    radius, 
    G4TwoVector               centre, 
    std::vector<G4TwoVector>  solid, 
    int                       chamferIndex)
{
  // check wrt the side before the chamfer 
  // -------------------------------------

  //FIXME this is bugged I don't understand why  !!!
  G4TwoVector vectorFromFirstEdgeOFChamfer = centre - solid.at (2 * chamferIndex) ;

  int index = 2 * chamferIndex - 1 ;
  if (index < 0) index += 8 ;

  G4TwoVector sideDirection = solid.at (index) - solid.at (2 * chamferIndex) ;
  sideDirection *= 1. / sideDirection.mag () ; 
  G4TwoVector sideOrtogonal = sideDirection ;
  sideOrtogonal.setX (sideDirection.y ()) ;
  sideOrtogonal.setY (-1 * sideDirection.x ()) ;

  float distance = sideOrtogonal.dot (vectorFromFirstEdgeOFChamfer) ;
  if (distance < radius) return false ; 
//  if (fabs (distance) < radius) return false ;

  // check wrt the side after the chamfer 
  // -------------------------------------

  //FIXME this works    !!!!!

  G4TwoVector vectorFromSecondEdgeOFChamfer = centre - solid.at (2 * chamferIndex + 1) ;

  index = 2 * chamferIndex + 2 ;
  if (index > 7) index -= 8 ;
  sideDirection = solid.at (2 * chamferIndex + 1) - solid.at (index) ;
  sideDirection *= 1. / sideDirection.mag () ; 
  sideOrtogonal = sideDirection ;
  sideOrtogonal.setX (sideDirection.y ()) ;
  sideOrtogonal.setY (-1 * sideDirection.x ()) ;
  
  distance = sideOrtogonal.dot (vectorFromSecondEdgeOFChamfer) ;
  // since the ortogonal aims outwards with respect to the crystal model, 
  // the distance should be at least minus radius
  if (distance < radius) return false ;

  return true ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4TwoVector 
DetectorConstruction::centerOfTheFirstFibreOnSecondLayer (
    std::pair<G4TwoVector, G4TwoVector> &  theChamfer, 
    const float &                          fiberExternalRadius, 
    G4TwoVector                            previousLayerStart)
{
  // assume that the chamfer coordinates are given counter-clockwise
  // so the orthogonal vector aims towards the exterior of the chamfer

  G4TwoVector chamferDirection = theChamfer.second - theChamfer.first ;
  chamferDirection *= 1. / chamferDirection.mag () ; 
  G4TwoVector chamferOrtogonal = chamferDirection ;
  chamferOrtogonal.setX (chamferDirection.y ()) ;
  chamferOrtogonal.setY (-1 * chamferDirection.x ()) ;

  return previousLayerStart                                            // the starting point
         + fiberExternalRadius * chamferDirection                       
         + (fiberExternalRadius * 1.73205080757) * chamferOrtogonal ; 
}


