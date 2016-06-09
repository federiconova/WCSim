#ifndef WCSimPrimaryGeneratorMessenger_h
#define WCSimPrimaryGeneratorMessenger_h 1

class WCSimPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"

class WCSimPrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* mpga);
  ~WCSimPrimaryGeneratorMessenger();
  
 public:
  void     SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
  
 private:
  WCSimPrimaryGeneratorAction* myAction;
  
 private: //commands
  G4UIdirectory*      mydetDirectory;
  G4UIdirectory*               CRYDir;
  G4UIcmdWithAString* genCmd;
  G4UIcmdWithAString* fileNameCmd;
  G4UIcmdWithAString*          FileCmd;
  G4UIcmdWithAString*          InputCmd;
  G4UIcmdWithoutParameter*     UpdateCmd;
  G4UIcmdWithABool*            fuseCRYCmd;
  std::string* MessInput;
  
};

#endif


