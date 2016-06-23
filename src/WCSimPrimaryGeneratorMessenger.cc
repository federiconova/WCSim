#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{


  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  //T. Akiri: Addition of laser
  genCmd->SetGuidance(" Available generators : muline, normal, laser, CRY");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  genCmd->SetCandidates("muline normal laser CRY");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  CRYDir = new G4UIdirectory("/CRY/");
  CRYDir->SetGuidance("CRY initialization");
  FileCmd = new G4UIcmdWithAString("/CRY/file",this);
  FileCmd->SetGuidance("This reads the CRY definition from a file");
  FileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  InputCmd = new G4UIcmdWithAString("/CRY/input",this);
  InputCmd->SetGuidance("CRY input lines");
  InputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/CRY/update",this);
  UpdateCmd->SetGuidance("Update CRY definition.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed the CRY definition.");
  UpdateCmd->AvailableForStates(G4State_Idle);

  MessInput = new std::string;


}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
  delete CRYDir;
  delete InputCmd;
  delete UpdateCmd;
  delete FileCmd;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    if (newValue == "muline")
    {
      myAction->SetMulineEvtGenerator(true);
      myAction->SetNormalEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetCRYEvtGenerator(false);
    }
    else if ( newValue == "normal")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetNormalEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetCRYEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetNormalEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetCRYEvtGenerator(false);
    }
    else if ( newValue == "CRY")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetNormalEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetCRYEvtGenerator(true);
    }
  }

  if( command == fileNameCmd )
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }

  if( command == InputCmd )
   {
     myAction->InputCRY();
     (*MessInput).append(newValue);
     (*MessInput).append(" ");
   }

  if( command == UpdateCmd )
   {
     myAction->UpdateCRY(MessInput);
     *MessInput = "";
   }

  if( command == FileCmd )
   { myAction->CRYFromFile(newValue); }


}

G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingNormalEvtGenerator())
      { cv = "normal"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingCRYEvtGenerator())
      { cv = "CRY"; }   
  }
  
  return cv;
}

