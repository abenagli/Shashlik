// Martin Goettlich @ DESY
//
 
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "MyMaterials.hh"
#include "CreateTree.hh"
#include "PrimaryGeneratorAction.hh"

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"


EventAction::EventAction()
{
}
 
EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  if( evtNb%1 == 0 ) 
  {
    G4cout << "---> Begin of Event: " << evtNb << G4endl;
  }
  
  CreateTree::Instance() -> Clear();
  
  // --------------------- STORE INFO FOR X_0 / R_M ----------------------------- //
  int Radial_nSteps       = 5000;
  int Longitudinal_nSteps = 1000;
  CreateTree::Instance() -> Radial_stepLength       = 1000. / Radial_nSteps;       // in mm
  CreateTree::Instance() -> Longitudinal_stepLength = 5000. / Longitudinal_nSteps; // in mm
  
  // --------------------- STORE PRIMARY VERTEX ----------------------------- //
  
  if(CreateTree::Instance() -> Init_data())
  {
    G4double InitEnergy   = evt -> GetPrimaryVertex() -> GetPrimary() -> GetTotalEnergy();
    G4ThreeVector InitPos = evt -> GetPrimaryVertex() -> GetPosition();
    G4ThreeVector InitDir = evt -> GetPrimaryVertex() -> GetPrimary() -> GetMomentumDirection();
    
    CreateTree::Instance() -> InitialEnergy = InitEnergy/GeV;
    
    CreateTree::Instance() -> InitialPositionX = InitPos[0]/mm;
    CreateTree::Instance() -> InitialPositionY = InitPos[1]/mm;
    CreateTree::Instance() -> InitialPositionZ = InitPos[2]/mm;
    
    CreateTree::Instance() -> InitalMomentumDirectionX = InitDir[0];
    CreateTree::Instance() -> InitalMomentumDirectionY = InitDir[1];
    CreateTree::Instance() -> InitalMomentumDirectionZ = InitDir[2];
  }
    
  // -------------------- INSTANCE RUN/EVENT IN TREE ---------------------- //
  CreateTree::Instance() -> Event = evt -> GetEventID();
}

void EventAction::EndOfEventAction(const G4Event* evt)
{ 
  CreateTree::Instance()->Fill();
}

