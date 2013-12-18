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


EventAction::EventAction ()
{
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


EventAction::~EventAction ()
{
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::BeginOfEventAction (const G4Event* evt)
{
  G4int evtNb = evt->GetEventID () ;
  if ( evtNb%1 == 0 ) 
  {
    G4cout << "---> Begin of Event: " << evtNb << G4endl ;
  }
  
  CreateTree::Instance ()->Clear () ;
  
  // INSTANCE RUN/EVENT IN TREE
  CreateTree::Instance ()->Event = evt->GetEventID () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::EndOfEventAction (const G4Event* evt)
{ 
  CreateTree::Instance ()->Fill () ;
}

