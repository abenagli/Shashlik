#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "ConfigFile.hh"



class MyMaterials
{
private:
  
public:
  MyMaterials();
  ~MyMaterials();
  
  static G4Material* Vacuum();
  static G4Material* Air();
  static G4Material* Water();
  static G4Material* Silicon();
  static G4Material* Brass();
  static G4Material* Tungsten();
  static G4Material* Quartz();
  static G4Material* OpticalGrease();
  static G4Material* LSO();
  static G4Material* LYSO();
  static G4Material* LYSO_lowLY();
  static G4Material* LuAG_Ce();
  static G4Material* LuAG_Pr();
  static G4Material* PWO();
  static G4Material* DSBCe();
  static G4Material* SiO2Ce();
  
  static G4double fromNmToEv(G4double wavelength);
  static G4double fromEvToNm(G4double energy);
  static G4double CalculateSellmeier(int size, G4double indexZero, G4double *nVec, G4double *lVec, G4double wavelength);
};
