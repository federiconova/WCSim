#include "WCSimMultiPMTParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCSimMultiPMTParameterisation::WCSimMultiPMTParameterisation(
        G4int    noPmt,
        G4double apoth,          // Distance from the Z axis
	std::vector<G4int>& vNiC,
	std::vector<G4double>& vAlpha,
	std::vector<G4int>& vCircle,
	G4double height)
 : G4VPVParameterisation()
{
   fNoPmt	= noPmt;
   fApothema	= apoth;
   vNiCLocal	= vNiC;
   vAlphaLocal	= vAlpha;
   vCircleLocal = vCircle;
   fHeight      = height;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCSimMultiPMTParameterisation::~WCSimMultiPMTParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WCSimMultiPMTParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
   // Note: copyNo will start with zero!
  G4RotationMatrix* rotm = new G4RotationMatrix(); // Rotation matrix for each chamber
  G4double angle = ((vNiCLocal[vCircleLocal[copyNo]]-2)*CLHEP::pi/vNiCLocal[vCircleLocal[copyNo]]); // Internal angle of each polygon
  G4ThreeVector origin(0,0,0); 
  //const G4double chambSize = fApothema/std::tan(angle/2);

  origin.setRThetaPhi(fApothema,CLHEP::halfpi-vAlphaLocal[vCircleLocal[copyNo]],copyNo*(CLHEP::pi-angle));
  // TF: Positioning vector, for an offset single hemisphere (for now)
  origin.setZ(origin.getZ()+fHeight);

  // TF: should be better documented
  //rotm->rotateX(CLHEP::halfpi);
  //rotm->rotateY((CLHEP::halfpi-(copyNo*angle))); //TF: strange as Y is a symmetry axis now, so rotation doesn't matter, unless Y' is meant
  //rotm->rotateX((std::pow(-1,copyNo)*vAlphaLocal[vCircleLocal[copyNo]]));
    
  // rotation of mother volume wrt daughter, hence minus sign.
  rotm->rotateZ(-origin.getPhi());
  rotm->rotateY(-acos((origin.getZ()-fHeight)/fApothema)); //over Y', origin.getTheta() is incorrect!

  /*
  std::cout << copyNo << " " << origin.getX() << " " << origin.getY() << " " << origin.getZ() << std::endl;
  std::cout << copyNo << " " << origin.getPhi() << " " << origin.getTheta() << std::endl;
  std::cout << copyNo << " " << rotm->getPhi() << " " << rotm->getTheta() << std::endl;
  std::cout << copyNo << " " << acos((origin.getZ()-fHeight)/fApothema) << " " << origin.getTheta() << std::endl; */
  physVol->SetTranslation(origin);
  physVol->SetRotation(rotm);
  //physVol->SetRotation(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void WCSimMultiPMTParameterisation::ComputeDimensions
(G4Tubs& pmtidmulti, const G4int copyNo, const G4VPhysicalVolume*) const
{
  // Note: copyNo will start with zero!

  G4double angle = ((fNoPmt-2)*pi/fNoPmt);
  const G4double chambSize = fApothema/std::tan(angle/2);
  pmtidmulti.SetZHalfLength(chambSize/5);
  pmtidmulti.SetOuterRadius(chambSize);
  G4cout << "Chamber radius = " << chambSize << " cm" << G4endl;

}
*/
