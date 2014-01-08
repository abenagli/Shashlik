#include "CreateTree.hh"
#include <cassert>


using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if ( fInstance )
  {
    return ;
  }

  this->fInstance = this ;
  this->fname     = name ;
  this->ftree     = new TTree (name,name) ;
  
  this->GetTree ()->Branch ("Event",                  &this->Event,                  "Event/I") ;
  this->GetTree ()->Branch ("totalPhLengthInChamfer", &this->totalPhLengthInChamfer, "totalPhLengthInChamfer[4]/F") ;
  this->GetTree ()->Branch ("numPhotonsInChamfer",    &this->numPhotonsInChamfer,    "numPhotonsInChamfer[4]/I") ;
  
  this->Clear () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::~CreateTree ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


/**
Loop on the container of the single photon total track in fibers cores.
Check the correctness of the chamfer ID assignment.
Count the number of photons for each chamfer.
The total length of photons in each chamfer is already saved in totalPhLengthInChamfer,
for historical reasons.
*/
int CreateTree::Fill () 
{ 
  for (std::map <int, std::pair<int, float> >::const_iterator iMap = fsingleGammaInfo.begin () ;
       iMap != fsingleGammaInfo.end () ;
       ++iMap)
    {
      assert (iMap->second.first < 4) ;
      assert (iMap->second.first >= 0) ;
      ++numPhotonsInChamfer[iMap->second.first] ;
    }
  return this->GetTree ()->Fill () ; 
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


bool CreateTree::Write ()
{
  TString filename = this->GetName () ;
  filename+=".root" ;
  TFile* file = new TFile (filename, "RECREATE") ;
  this->GetTree ()->Write () ;
  file->Write () ;
  file->Close () ;
  return true ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


/**
For each photon, record the full path length of that particular photon,
according to the length passed to the function.
In SteppingAction.cc, this length is calculated as the one traveled
in the core of each fiber, therefore this is the total length traveled
in the fibers cores, for a single photon, per event.
*/
void CreateTree::addPhoton (int trackId, float length, int chamferId)
{
  if (fsingleGammaInfo.find (trackId) == fsingleGammaInfo.end ())
    {
      fsingleGammaInfo[trackId] = std::pair<int, float> (chamferId, length) ;
    }
  else  
    {
      fsingleGammaInfo[trackId].second += length ;
    }
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::Clear ()
{
  Event	= 0 ;
  for (int i = 0 ; i < 4 ; ++i) 
    {
      totalPhLengthInChamfer[i] = 0. ;
      numPhotonsInChamfer[i] = 0. ;
    }
  fsingleGammaInfo.clear () ;
}
