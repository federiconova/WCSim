//
// Created by zsoldos on 24/01/18.
//

#ifndef WCSIM_WCSIMWLSPROPERTIES_HH
#define WCSIM_WCSIMWLSPROPERTIES_HH

#include <TGraph.h>
#include <TH1F.h>
#include "globals.hh"

class WCSimWLSProperties {

 public:
  virtual G4String  GetWLSPlateName()=0;
  virtual G4int     GetNumEntries()=0;
  virtual G4double* GetPhotonWL()=0;
  virtual G4double* GetPhotonEnergy()=0;
  virtual G4double* GetRIndex()=0;
  virtual G4double* GetAbs()=0;
  virtual G4double* GetEm()=0;

 public:
  TGraph* GetgAbs(){ return gAbs;};
  TGraph* GetgEm() { return gEm;};
  TGraph *gAbs;
  TGraph *gEm;
  TH1F *emissionHist;
};

class EljenEJ286 : public WCSimWLSProperties{

 public:
  EljenEJ286();
  ~EljenEJ286();

 public:
  G4String  GetWLSPlateName();
  G4int     GetNumEntries();
  G4double* GetPhotonWL();
  G4double* GetPhotonEnergy();
  G4double* GetRIndex();

  G4int     GetNumEntries_ABS();
  G4double* GetPhotonEnergy_ABS();
  G4double* GetAbs();

  G4int     GetNumEntries_EM();
  G4double* GetPhotonEnergy_EM();
  G4double* GetEm();

  void SetgAbs();
  void SetgEm();

  void SethEm();
};

class Kuraray : public WCSimWLSProperties{

 public:
  Kuraray();
  ~Kuraray();

 public:
  G4String  GetWLSPlateName();
  G4int     GetNumEntries();
  G4double* GetPhotonWL();
  G4double* GetPhotonEnergy();
  G4double* GetRIndex();

  G4int     GetNumEntries_ABS();
  G4double* GetPhotonEnergy_ABS();
  G4double* GetAbs();

  G4int     GetNumEntries_EM();
  G4double* GetPhotonEnergy_EM();
  G4double* GetEm();

  void SetgAbs();
  void SetgEm();

  void SethEm();
};

class Lablogic : public WCSimWLSProperties{

 public:
  Lablogic();
  ~Lablogic();

 public:
  G4String  GetWLSPlateName();
  G4int     GetNumEntries();
  G4double* GetPhotonWL();
  G4double* GetPhotonEnergy();
  G4double* GetRIndex();

  G4int     GetNumEntries_ABS();
  G4double* GetPhotonEnergy_ABS();
  G4double* GetAbs();

  G4int     GetNumEntries_EM();
  G4double* GetPhotonEnergy_EM();
  G4double* GetEm();

  void SetgAbs();
  void SetgEm();

  void SethEm();
};

#endif //WCSIM_WCSIMWLSPROPERTIES_HH
