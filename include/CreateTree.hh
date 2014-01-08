#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"



class CreateTree
{
private:
  
  TTree*  ftree ;
  TString fname ;
  std::map <int, std::pair<int, float> > fsingleGammaInfo ;
  //        photonID       chamferID  lengthInChamfer 
  // lengthInChamfer is redundant
  
public:
  
  CreateTree (TString name) ;
  ~CreateTree () ;
  
  TTree*             GetTree  () const { return ftree ; } ;
  TString            GetName  () const { return fname ; } ;
  int                Fill     () ;
  bool               Write    () ;
  void               Clear    () ;
  static CreateTree* Instance () { return fInstance ; } ;
  
  // feed the info of each single photon to the tree
  void               addPhoton (int trackId, float length, int chamferId) ;

  static CreateTree* fInstance ;
  
  int Event ;
  float totalPhLengthInChamfer[4] ;           // total photons length in chamfers
  int   numPhotonsInChamfer[4] ;              // number of photons in chamfers

} ;
