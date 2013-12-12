    
#include "SteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4UnitsTable.hh"
#include "CreateTree.hh"
#include "MyMaterials.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"



SteppingAction::SteppingAction()
{}



SteppingAction::~SteppingAction()
{}



void SteppingAction::UserSteppingAction(const G4Step * theStep)
{
//  G4Track* theTrack = theStep->GetTrack();
//  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  
  G4VPhysicalVolume* thePrePV  = thePrePoint ->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName  = ""; if( thePrePV )  thePrePVName  = thePrePV  -> GetName();
  G4String thePostPVName = ""; if( thePostPV ) thePostPVName = thePostPV -> GetName();
  
  // optical photon
  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition ())
  { 
      // Check that the step starts and ends in any of the core volumes of any fibers.
      // This includes both the inner and outer unphysical sub-volumes of the core.
      
      std::size_t pos = thePrePVName.find ("FiberCore") ;
      if (pos == std::string::npos) continue ;
      pos   = thePostPVName.find ("FiberCore") ;
      if (pos == std::string::npos) continue ;

      // Find the chamfer where the photon is traveling in.
      for (int i = 0 ; i < 4 ; ++i)
        {
          string num_s = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
          pos = thePrePVName.find (num_s) ;
          if (pos == std::string::npos) continue ;
          float length = theStep.GetStepLength () ;
          // attribute the length to the chamfer
          // sum the lengths for each photon separately
        }
   } // optical photon
  
  
}
