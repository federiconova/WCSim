#include "WCSimDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SDManager.hh"
#include "WCSimWCSD.hh"
#include "WCSimPMTObject.hh"

#include "G4SystemOfUnits.hh"

//PMT logical volume construction.

WCSimDetectorConstruction::PMTMap_t WCSimDetectorConstruction::PMTLogicalVolumes;

G4LogicalVolume* WCSimDetectorConstruction::ConstructPMT(G4String PMTName, G4String CollectionName)
{
  PMTKey_t key(PMTName,CollectionName);

  PMTMap_t::iterator it = PMTLogicalVolumes.find(key);
  if (it != PMTLogicalVolumes.end()) {
      //G4cout << "Restore PMT" << G4endl;
      return it->second;
  }


  //G4cout << "Create PMT" << G4endl;


if (Vis_Choice == "RayTracer"){
    // Blue wireframe visual style
    // Used in the RayTracer visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
}

else
   { // Gray wireframe visual style
    // used in OGLSX visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  WCPMTVisAtt->SetForceWireframe(true);}

  G4double expose;
  G4double radius;
  G4double glassThickness;
  
  WCSimPMTObject *PMT = GetPMTPointer(CollectionName);
  expose = PMT->GetExposeHeight();
  radius = PMT->GetRadius();
  glassThickness = PMT->GetPMTGlassThickness();

  G4double sphereRadius = (expose*expose+ radius*radius)/(2*expose);
  G4double PMTOffset =  sphereRadius - expose;

  //All components of the PMT are now contained in a single logical volume logicWCPMT.
  //Origin is on the blacksheet, faces positive z-direction.
  
  //Need a volume to cut away excess behind blacksheet
  G4Box* solidCutOffTubs =
      new G4Box(    "cutOffTubs",
            sphereRadius+1.*cm,
            sphereRadius+1.*cm,
            PMTOffset);

  G4LogicalVolume* logicWCPMT;

   if( PMTName == "TriangularTile3inch" ){
     G4cout << " qqq building PMT " << PMTName << G4endl;

     logicWCPMT = new G4LogicalVolume(NULL, NULL, "WCPMT", 0, 0, 0);
     BuildWLSplate(radius, expose, sphereRadius, PMTOffset, solidCutOffTubs, logicWCPMT);

     }else{

    G4double PMTHolderZ[2] = {0, expose};
    G4double PMTHolderR[2] = {radius, radius};
    G4double PMTHolderr[2] = {0,0};
    G4Polycone* solidWCPMT = 
      new G4Polycone("WCPMT",                    
		     0.0*deg,
		     360.0*deg,
		     2,
		     PMTHolderZ,
		     PMTHolderr, // R Inner
		     PMTHolderR);// R Outer
    logicWCPMT =
      new G4LogicalVolume(    solidWCPMT,
			      G4Material::GetMaterial("Water"),
			      "WCPMT",
			      0,0,0);
     }

if (Vis_Choice == "RayTracer"){
// Makes the volume containing the PMT visible, solid, and forces the auxiliary edges to be viewed.
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

    logicWCPMT->SetVisAttributes(WCPMTVisAtt);}

else{
// Makes the volume containg the PMT invisible for normal visualization
    logicWCPMT->SetVisAttributes(G4VisAttributes::Invisible);}


  //Create PMT Interior
  G4Sphere* tmpSolidInteriorWCPMT =
      new G4Sphere(    "tmpInteriorWCPMT",
                       0.0*m,(sphereRadius-glassThickness),
                       0.0*deg,360.0*deg,
                       0.0*deg,90.0*deg);

  G4SubtractionSolid* solidInteriorWCPMT =
      new G4SubtractionSolid(    "InteriorWCPMT",
                    tmpSolidInteriorWCPMT,
                    solidCutOffTubs);

  // "Air" here is not true air, but a modified material
  // with n = 1 and a very short absorption length
  G4LogicalVolume* logicInteriorWCPMT =
    new G4LogicalVolume(    solidInteriorWCPMT,
                    G4Material::GetMaterial("Air"),
                    "InteriorWCPMT",
                    0,0,0);

  G4VPhysicalVolume* physiInteriorWCPMT =
      new G4PVPlacement(0,
                  G4ThreeVector(0, 0, -1.0*PMTOffset),
                  logicInteriorWCPMT,
                  "InteriorWCPMT",
                  logicWCPMT,
                  false,
                  0);

if (Vis_Choice == "RayTracer"){
// Adding color and forcing the inner portion of the PMT's to be viewed
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

  logicInteriorWCPMT->SetVisAttributes(WCPMTVisAtt);}

else {
// Making the inner portion of the detector invisible for OGLSX visualization
  logicInteriorWCPMT->SetVisAttributes(G4VisAttributes::Invisible);}


  //Create PMT Glass Face
  G4Sphere* tmpGlassFaceWCPMT =
      new G4Sphere(    "tmpGlassFaceWCPMT",
                       (sphereRadius-glassThickness),
                       sphereRadius,
                       0.0*deg,360.0*deg,
                       0.0*deg,90.0*deg);
  
  G4SubtractionSolid* solidGlassFaceWCPMT =
      new G4SubtractionSolid(    CollectionName,
                                 tmpGlassFaceWCPMT,
                                 solidCutOffTubs); 

  G4LogicalVolume *logicGlassFaceWCPMT =
    new G4LogicalVolume(    solidGlassFaceWCPMT,
                            G4Material::GetMaterial("Glass"),
                            CollectionName,
                            0,0,0);

  G4VPhysicalVolume* physiGlassFaceWCPMT =
      new G4PVPlacement(0,
                        G4ThreeVector(0, 0, -1.0*PMTOffset),
                        logicGlassFaceWCPMT,
                        CollectionName,
                        logicWCPMT,
                        false,
                        0,
                        checkOverlaps);

// For either visualization type, logicGlassFaceWCPMT will either be visible or invisible depending on which
// line is commented at the end of the respective if statements

  if (Vis_Choice == "OGLSX")
   { // Gray wireframe visual style
    // used in OGLSX visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  WCPMTVisAtt->SetForceWireframe(true);
  //logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
  logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

  if (Vis_Choice == "RayTracer"){
    // Blue wireframe visual style
    // Used in the RayTracer visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
  WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
  //logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);

  logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

  else
   { // Gray wireframe visual style
    // used in OGLSX visualizer
  G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  WCPMTVisAtt->SetForceWireframe(true);
  //logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
  logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

  // Instantiate a new sensitive detector 
  // and register this sensitive detector volume with the SD Manager. 
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDName = "/WCSim/";
  SDName += CollectionName;

  // If there is no such sensitive detector with that SDName yet,
  // make a new one
  if( ! SDman->FindSensitiveDetector(SDName, false) ) {
    
    aWCPMT = new WCSimWCSD(CollectionName,SDName,this );
    SDman->AddNewDetector( aWCPMT );
  }

  logicGlassFaceWCPMT->SetSensitiveDetector( aWCPMT );

  PMTLogicalVolumes[key] = logicWCPMT;

  //Add Logical Border Surface
  new G4LogicalBorderSurface("GlassCathodeSurface",
                             physiGlassFaceWCPMT,
                             physiInteriorWCPMT,
                             OpGlassCathodeSurface);

  return logicWCPMT;
}

void WCSimDetectorConstruction::BuildWLSplate(double PMT_radius, double PMT_height, double sphereRadius, double PMTOffset, G4Box* solidCutOffTubs, G4LogicalVolume* logicWCPMT){

  G4double PetalLength = 35.*cm;
  G4double PetalHalfThickness = 1.2*cm;
  G4double PetalHalfWidth = PMT_radius;
  G4double CladdingThickness = 1.*mm;
  
  G4cout << " qqq PetalHalfWidth " << PetalHalfWidth/m << " m, PetalLength " << PetalLength/m << " m, PetalHalfThickness " << PetalHalfThickness/m  << " m " << G4endl;
  
  G4double Trapezoid_dx1 = PetalHalfWidth; //  half-length along x at z = - dz
  G4double Trapezoid_dx2 = 0.; //  half-length along x at z = + dz
  G4double Trapezoid_dy1 = PetalHalfThickness; //  half-length along y at z = - dz
  G4double Trapezoid_dy2 = PetalHalfThickness; //  half-length along y at z = + dz
  G4double Trapezoid_dz = PetalLength/2.; //  half-length along z
  G4cout << " qqq Trapezoid: dx1 " << Trapezoid_dx1/m << " m, dx2 " << Trapezoid_dx2/m << " m, dy1 " << Trapezoid_dy1/m << " m, dy2 " << Trapezoid_dy2/m << " m, dz " << Trapezoid_dz/m << " m " << G4endl;
  
  
  G4RotationMatrix* TrapezoidRotation = new G4RotationMatrix();
  TrapezoidRotation->rotateY(180.*deg);
  G4Transform3D TrapezoidTransform(*TrapezoidRotation, G4ThreeVector(0,0,-2.*Trapezoid_dz));
  G4Trd *OuterPetal_half = new G4Trd("Outer Petal Half", Trapezoid_dx1, Trapezoid_dx2, Trapezoid_dy1, Trapezoid_dy2, Trapezoid_dz);
  G4UnionSolid *OuterPetal = new G4UnionSolid("Outer Petal", OuterPetal_half, OuterPetal_half, TrapezoidTransform);
  
  
  G4Trd *OuterCladding_half = new G4Trd("Outer Cladding Half", Trapezoid_dx1+CladdingThickness, Trapezoid_dx2+CladdingThickness, Trapezoid_dy1, Trapezoid_dy2, Trapezoid_dz);
  G4UnionSolid *OuterCladding = new G4UnionSolid("Outer Cladding", OuterCladding_half, OuterCladding_half, TrapezoidTransform);
  
  
  G4Sphere* tmp_glass_outer_surface =
    new G4Sphere(    "tmp_glass_outer_surface",
		     0.0*m,sphereRadius,
		     0.0*deg,360.0*deg,
		     0.0*deg,90.0*deg);
  
  G4SubtractionSolid* glass_outer_surface =
    new G4SubtractionSolid(    "glass_outer_surface",
			       tmp_glass_outer_surface,
			       solidCutOffTubs);


    G4RotationMatrix* PmtRotation = new G4RotationMatrix();
    PmtRotation->rotateX(90.*deg);
    G4Transform3D PmtTransform(*PmtRotation, G4ThreeVector(0,PetalHalfThickness + PMTOffset,-Trapezoid_dz));
    G4SubtractionSolid * WLSplate = new G4SubtractionSolid("wls_plate", OuterPetal, glass_outer_surface, PmtTransform);

    G4SubtractionSolid * CladdingWithoutPetal = new G4SubtractionSolid("Cladding With Floor",OuterCladding,OuterPetal,0,G4ThreeVector(0,0,0));

    G4SubtractionSolid * cladding = new G4SubtractionSolid("cladding",CladdingWithoutPetal, glass_outer_surface,PmtTransform);



    G4LogicalVolume * WLSplate_log = new G4LogicalVolume(WLSplate,G4Material::GetMaterial("WLS_PVT"),"wls_plate_log");

    G4LogicalVolume * cladding_log = new G4LogicalVolume(cladding,G4Material::GetMaterial("Pethylene"),"cladding",0,0,0);

    G4double PMTHolderZ[2] = {0, PMT_height};
    G4double PMTHolderR[2] = {PMT_radius, PMT_radius};
    G4double PMTHolderr[2] = {0,0};
    G4Polycone* solidWCPMT = 
      new G4Polycone("WCPMT",                    
		     0.0*deg,
		     360.0*deg,
		     2,
		     PMTHolderZ,
		     PMTHolderr, // R Inner
		     PMTHolderR);// R Outer
    logicWCPMT->SetSolid(solidWCPMT);
    logicWCPMT->SetMaterial(G4Material::GetMaterial("Water"));

}
