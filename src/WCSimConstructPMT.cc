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
#include "G4LogicalSkinSurface.hh"

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

     logicWCPMT = new G4LogicalVolume(NULL, NULL, "WCPMT", 0, 0, 0);
     BuildWLSplateTriangularTile(radius, expose, sphereRadius, PMTOffset, solidCutOffTubs, logicWCPMT);

   }else if( PMTName == "FastStar3inch" ){

     logicWCPMT = new G4LogicalVolume(NULL, NULL, "WCPMT", 0, 0, 0);
     BuildWLSplateFastStar(radius, expose, sphereRadius, PMTOffset, solidCutOffTubs, logicWCPMT);

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

void WCSimDetectorConstruction::BuildWLSplateTriangularTile(double PMT_radius, double PMT_height, double sphereRadius, double PMTOffset, G4Box* solidCutOffTubs, G4LogicalVolume* logicWCPMT){

  G4double PetalLength = 35.*cm;
  G4double PetalHalfThickness = 1.2*cm;
  G4double PetalHalfWidth = PMT_radius;
  G4double CladdingThickness = 1.*mm;
  
  G4double Trapezoid_dx1 = PetalHalfWidth; //  half-length along x at z = - dz
  G4double Trapezoid_dx2 = 0.; //  half-length along x at z = + dz
  G4double Trapezoid_dy1 = PetalHalfThickness; //  half-length along y at z = - dz
  G4double Trapezoid_dy2 = PetalHalfThickness; //  half-length along y at z = + dz
  G4double Trapezoid_dz = PetalLength/2.; //  half-length along z
  
  
  G4RotationMatrix* TrapezoidRotation = new G4RotationMatrix();
  TrapezoidRotation->rotateY(180.*deg);
  G4Transform3D TrapezoidTransform(*TrapezoidRotation, G4ThreeVector(0,0,-2.*Trapezoid_dz));
  G4Trd *OuterPetal_half = new G4Trd("Outer Petal Half", Trapezoid_dx1, Trapezoid_dx2, Trapezoid_dy1, Trapezoid_dy2, Trapezoid_dz);
  G4UnionSolid *OuterPetal = new G4UnionSolid("Outer Petal", OuterPetal_half, OuterPetal_half, TrapezoidTransform);
  
  
  G4Trd *OuterCladding_half = new G4Trd("Outer Cladding Half", Trapezoid_dx1+CladdingThickness, Trapezoid_dx2+CladdingThickness, Trapezoid_dy1 + CladdingThickness/2., Trapezoid_dy2 + CladdingThickness/2., Trapezoid_dz);
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
    G4Transform3D PmtTransform_for_WLSplate(*PmtRotation, G4ThreeVector(0,PetalHalfThickness + PMTOffset + CladdingThickness,-Trapezoid_dz));
    G4Transform3D PmtTransform_for_cladding(*PmtRotation, G4ThreeVector(0,PetalHalfThickness + PMTOffset + CladdingThickness/2.,-Trapezoid_dz));
    G4SubtractionSolid * WLSplate = new G4SubtractionSolid("wls_plate", OuterPetal, glass_outer_surface, PmtTransform_for_WLSplate);

    G4SubtractionSolid * CladdingWithoutPetal = new G4SubtractionSolid("Cladding With Floor",OuterCladding,OuterPetal,0,G4ThreeVector(0,-CladdingThickness/2.,0));

    G4SubtractionSolid * cladding = new G4SubtractionSolid("cladding",CladdingWithoutPetal, glass_outer_surface,PmtTransform_for_cladding);



    G4LogicalVolume * WLSplate_log = new G4LogicalVolume(WLSplate,G4Material::GetMaterial("WLS_PVT"),"wls_plate_log");

    //    G4LogicalVolume * cladding_log = new G4LogicalVolume(cladding,G4Material::GetMaterial("Pethylene"),"cladding",0,0,0);
    G4LogicalVolume * cladding_log = new G4LogicalVolume(cladding,G4Material::GetMaterial("Tyvek"),"cladding",0,0,0);




#if 1 //  F. Nova Comment this block to avoid building WLS plate and its cladding
       new G4PVPlacement(PmtRotation, G4ThreeVector(0, 0 + Trapezoid_dz, 0 + PetalHalfThickness + CladdingThickness), WLSplate_log, "wlsplate", logicWCPMT , false, 0); 
    
       new G4PVPlacement(PmtRotation,G4ThreeVector(0, 0 + Trapezoid_dz, 0 + PetalHalfThickness + CladdingThickness/2.), cladding_log,"cladding",logicWCPMT,false,0);
#endif


  //**Create logical skin surfaces
  // F. Nova Need this step to enforce cladding reflectivity
  new G4LogicalSkinSurface("cladding_surf",   cladding_log,   OpCladdingSurface);


    G4double PMTHolderZ[2] = {0, PMT_height};
    //G4double PMTHolderR[2] = {PMT_radius, PMT_radius};
    G4double PMTHolderR[2] = {PetalLength, PetalLength};
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


void WCSimDetectorConstruction::BuildWLSplateFastStar(double PMT_radius, double PMT_height, double sphereRadius, double PMTOffset, G4Box* solidCutOffTubs, G4LogicalVolume* logicWCPMT){


  G4double PetalLength = 35*cm;
  G4double PetalHalfThickness = 1.2*cm;
  G4double CladdingThickness = 1.*mm;
  
  G4double n_triangles_per_plane = 8;
  G4double angle_of_triangle = 2.*pi/n_triangles_per_plane; //  45 degrees
  G4double half_angle_of_triangle = angle_of_triangle/2.; //  45/2 degrees
  
  G4double BottomPetalProtrudingLength = PetalLength - PMT_radius*cos(half_angle_of_triangle);
  G4double BottomPetalProtrudedBase = PMT_radius*sqrt(2.*(1. - cos(angle_of_triangle)));
  G4double BottomPetalBase = BottomPetalProtrudedBase*PetalLength/BottomPetalProtrudingLength;
  G4double BottomPetalHalfWidth = BottomPetalBase/2.;
  
  G4double BottomTrapezoid_dx1 = BottomPetalHalfWidth; //  half-length along x at z = - dz
  G4double BottomTrapezoid_dx2 = 0.; //  half-length along x at z = + dz
  G4double BottomTrapezoid_dy1 = PetalHalfThickness; //  half-length along y at z = - dz
  G4double BottomTrapezoid_dy2 = PetalHalfThickness; //  half-length along y at z = + dz
  G4double BottomTrapezoid_dz = PetalLength/2.; //  half-length along z

  G4double TopPetalZ = 2.*PetalHalfThickness+2*CladdingThickness;

  G4double TopExpose = PMT_height - TopPetalZ;
  G4double TopRadius = sqrt(2.*TopExpose*sphereRadius - pow(TopExpose,2));
  G4double TopPetalProtrudingLength = PetalLength - TopRadius*cos(half_angle_of_triangle);
  G4double TopPetalProtrudedBase = TopRadius*sqrt(2.*(1. - cos(angle_of_triangle)));
  G4double TopPetalBase = TopPetalProtrudedBase*PetalLength/TopPetalProtrudingLength;
  G4double TopPetalHalfWidth = TopPetalBase/2.;
  
  
  G4double TopTrapezoid_dx1 = TopPetalHalfWidth; //  half-length along x at z = - dz
  G4double TopTrapezoid_dx2 = 0.; //  half-length along x at z = + dz
  G4double TopTrapezoid_dy1 = PetalHalfThickness; //  half-length along y at z = - dz
  G4double TopTrapezoid_dy2 = PetalHalfThickness; //  half-length along y at z = + dz
  G4double TopTrapezoid_dz = PetalLength/2.; //  half-length along z
  
  
  G4RotationMatrix* TrapezoidRotation = new G4RotationMatrix();
  TrapezoidRotation->rotateY(180.*deg);
  G4Transform3D BottomTrapezoidTransform(*TrapezoidRotation, G4ThreeVector(0,0,-2.*BottomTrapezoid_dz));
  G4Trd *OuterPetalBottom_half = new G4Trd("Outer Petal Bottom Half", BottomTrapezoid_dx1, BottomTrapezoid_dx2, BottomTrapezoid_dy1, BottomTrapezoid_dy2, BottomTrapezoid_dz);
  G4UnionSolid *OuterPetalBottom = new G4UnionSolid("Outer Petal Bottom", OuterPetalBottom_half, OuterPetalBottom_half, BottomTrapezoidTransform);


  G4Trd *OuterCladdingBottom_half = new G4Trd("Outer Cladding Bottom Half", BottomTrapezoid_dx1+CladdingThickness, BottomTrapezoid_dx2+CladdingThickness, BottomTrapezoid_dy1+CladdingThickness/2., BottomTrapezoid_dy2+CladdingThickness/2., BottomTrapezoid_dz);
  G4UnionSolid *OuterCladdingBottom = new G4UnionSolid("Outer Cladding Bottom", OuterCladdingBottom_half, OuterCladdingBottom_half, BottomTrapezoidTransform);

  G4Transform3D TopTrapezoidTransform(*TrapezoidRotation, G4ThreeVector(0,0,-2.*TopTrapezoid_dz));
  G4Trd *OuterPetalTop_half = new G4Trd("Outer Petal Top Half", TopTrapezoid_dx1, TopTrapezoid_dx2, TopTrapezoid_dy1, TopTrapezoid_dy2, TopTrapezoid_dz);
  G4UnionSolid *OuterPetalTop = new G4UnionSolid("Outer Petal Top", OuterPetalTop_half, OuterPetalTop_half, TopTrapezoidTransform);

  G4Trd *OuterCladdingTop_half = new G4Trd("Outer Cladding Top Half", TopTrapezoid_dx1+CladdingThickness, TopTrapezoid_dx2+CladdingThickness, TopTrapezoid_dy1+CladdingThickness/2., TopTrapezoid_dy2+CladdingThickness/2., TopTrapezoid_dz);
  G4UnionSolid *OuterCladdingTop = new G4UnionSolid("Outer Cladding Top", OuterCladdingTop_half, OuterCladdingTop_half, TopTrapezoidTransform);


  G4RotationMatrix* PetalRotation21 = new G4RotationMatrix();
  PetalRotation21->rotateY(-45.*deg);
  G4Transform3D BottomPetalTransform1(*PetalRotation21, G4ThreeVector(-BottomTrapezoid_dz*sin(angle_of_triangle),0,-BottomTrapezoid_dz*(1. - cos(angle_of_triangle))));
  G4UnionSolid * GrowingFlowerBottom1 = new G4UnionSolid("Growing Flower Bottom", OuterPetalBottom, OuterPetalBottom, BottomPetalTransform1 );
  G4UnionSolid * GrowingCladdingBottom1 = new G4UnionSolid("Growing Cladding Bottom", OuterCladdingBottom, OuterCladdingBottom, BottomPetalTransform1);

  G4RotationMatrix* PetalRotation22 = new G4RotationMatrix();
  PetalRotation22->rotateY(-90.*deg);
  G4Transform3D BottomPetalTransform2(*PetalRotation22, G4ThreeVector(-BottomTrapezoid_dz,0,-BottomTrapezoid_dz));
  G4UnionSolid * GrowingFlowerBottom = new G4UnionSolid("Growing Flower Bottom", GrowingFlowerBottom1, OuterPetalBottom, BottomPetalTransform2 );
  G4UnionSolid * GrowingCladdingBottom = new G4UnionSolid("Growing Cladding Bottom", GrowingCladdingBottom1, OuterCladdingBottom, BottomPetalTransform2);




  G4Transform3D TopPetalTransform1(*PetalRotation21, G4ThreeVector(-TopTrapezoid_dz*sin(angle_of_triangle),0,-TopTrapezoid_dz*(1. - cos(angle_of_triangle))));
  G4UnionSolid * GrowingFlowerTop1 = new G4UnionSolid("Growing Flower Top", OuterPetalTop, OuterPetalTop, TopPetalTransform1 );
  G4UnionSolid * GrowingCladdingTop1 = new G4UnionSolid("Growing Cladding Top", OuterCladdingTop, OuterCladdingTop, TopPetalTransform1 );

  G4Transform3D TopPetalTransform2(*PetalRotation22, G4ThreeVector(-TopTrapezoid_dz,0,-TopTrapezoid_dz));
  G4UnionSolid * GrowingFlowerTop = new G4UnionSolid("Growing Flower Top", GrowingFlowerTop1, OuterPetalTop, TopPetalTransform2 );
  G4UnionSolid * GrowingCladdingTop = new G4UnionSolid("Growing Cladding Top", GrowingCladdingTop1, OuterCladdingTop, TopPetalTransform2 );



  G4RotationMatrix* PetalRotation3 = new G4RotationMatrix();
  PetalRotation3->rotateY(-half_angle_of_triangle);
  G4Transform3D PetalTransform2(*PetalRotation3, G4ThreeVector(-TopTrapezoid_dz*sin(half_angle_of_triangle),-TopPetalZ,-BottomTrapezoid_dz+TopTrapezoid_dz*cos(half_angle_of_triangle)));
  G4UnionSolid * GrownFlower = new G4UnionSolid("Grown Flower", GrowingFlowerBottom, GrowingFlowerTop, PetalTransform2 );
  G4UnionSolid * GrownCladding = new G4UnionSolid("Grown Cladding", GrowingCladdingBottom, GrowingCladdingTop, PetalTransform2 );
  
  G4RotationMatrix* PmtRotation = new G4RotationMatrix();
  PmtRotation->rotateX(90.*deg);

  G4Transform3D PmtTransform_for_WLSplate(*PmtRotation, G4ThreeVector(0,PetalHalfThickness + PMTOffset + CladdingThickness,-BottomTrapezoid_dz));
  G4Transform3D PmtTransform_for_cladding(*PmtRotation, G4ThreeVector(0,PetalHalfThickness + PMTOffset + CladdingThickness/2.,-BottomTrapezoid_dz));

  G4Sphere* tmp_glass_outer_surface =
    new G4Sphere(    "tmp_glass_outer_surface",
		     0.0*m,sphereRadius,
		     0.0*deg,360.0*deg,
		     0.0*deg,90.0*deg);
  
  G4SubtractionSolid* glass_outer_surface =
    new G4SubtractionSolid(    "glass_outer_surface",
			       tmp_glass_outer_surface,
			       solidCutOffTubs);

  G4SubtractionSolid * WLSplate = new G4SubtractionSolid("wls_plate", GrownFlower, glass_outer_surface, PmtTransform_for_WLSplate);

  G4SubtractionSolid * CladdingWithoutPetal = new G4SubtractionSolid("Cladding With Floor",GrownCladding,GrownFlower,0,G4ThreeVector(0,-CladdingThickness/2.,0));

  G4SubtractionSolid * cladding = new G4SubtractionSolid("cladding",CladdingWithoutPetal, glass_outer_surface,PmtTransform_for_cladding);


  G4LogicalVolume * WLSplate_log = new G4LogicalVolume(WLSplate,G4Material::GetMaterial("WLS_PVT"),"wls_plate_log");
  
  //    G4LogicalVolume * cladding_log = new G4LogicalVolume(cladding,G4Material::GetMaterial("Pethylene"),"cladding",0,0,0);
  G4LogicalVolume * cladding_log = new G4LogicalVolume(cladding,G4Material::GetMaterial("Tyvek"),"cladding",0,0,0);

#if 1 //  F. Nova Comment this block to avoid building WLS plate and its cladding
       new G4PVPlacement(PmtRotation, G4ThreeVector(0, 0 + BottomTrapezoid_dz, 0 + PetalHalfThickness + CladdingThickness), WLSplate_log, "wlsplate", logicWCPMT , false, 0); 
    
       new G4PVPlacement(PmtRotation,G4ThreeVector(0, 0 + BottomTrapezoid_dz, 0 + PetalHalfThickness + CladdingThickness/2.), cladding_log,"cladding",logicWCPMT,false,0);
#endif


  //**Create logical skin surfaces
  // F. Nova Need this step to enforce cladding reflectivity
  new G4LogicalSkinSurface("cladding_surf",   cladding_log,   OpCladdingSurface);


    G4double PMTHolderZ[2] = {0, PMT_height};
    //G4double PMTHolderR[2] = {PMT_radius, PMT_radius};
    G4double PMTHolderR[2] = {PetalLength, PetalLength};
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
