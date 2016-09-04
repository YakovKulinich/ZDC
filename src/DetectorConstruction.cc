//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "SharedData.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <iostream>
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), m_checkOverlaps(true),
    m_moduleSizeX(0), m_moduleSizeY(0), m_moduleSizeZ(0),
    m_chamberSizeX(0), m_chamberSizeY(0), m_chamberSizeZ(0),
    m_matBox(0), m_matQuartz(0), m_matOil(0), m_matReflector(0),
    m_matAbsorber(0), m_matTop(0),
    m_solidWorld(0), m_logicWorld(0), m_physWorld(0),
    m_solidOil(0), m_logicOil(0), m_physOil(0),
    m_solidTop(0), m_logicTop(0), m_physTop(0),
    m_scoringVolume(0),
    m_sd(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( SharedData* sd )
  : G4VUserDetectorConstruction(), m_checkOverlaps(true),
    m_moduleSizeX(0), m_moduleSizeY(0), m_moduleSizeZ(0),
    m_chamberSizeX(0), m_chamberSizeY(0), m_chamberSizeZ(0),
    m_matBox(0), m_matQuartz(0), m_matOil(0), m_matReflector(0),
    m_matAbsorber(0), m_matTop(0),
    m_solidWorld(0), m_logicWorld(0), m_physWorld(0),
    m_solidOil(0), m_logicOil(0), m_physOil(0),
    m_solidTop(0), m_logicTop(0), m_physTop(0),
    m_scoringVolume(0),
    m_sd( sd )
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::InitializeParameters()
{

  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  m_moduleSizeX     = config->GetValue( "moduleSizeX", 90);
  m_moduleSizeY     = config->GetValue( "moduleSizeY", 240);
  m_moduleSizeZ     = config->GetValue( "moduleSizeZ", 150);

  G4double quartzThickness     = config->GetValue( "quartzThickness", 2);
  G4double moduleWallThickness = config->GetValue( "moduleWallThickness", 2);

  m_chamberSizeX  = m_moduleSizeX - 2*(moduleWallThickness + quartzThickness);
  m_chamberSizeY  = config->GetValue( "absorberHeight", 180 );
  m_chamberSizeZ  = m_moduleSizeZ - 2*moduleWallThickness;

  // Option to switch on/off checking of volumes overlaps
  //
  m_checkOverlaps = true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){

  if ( m_physWorld ) {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
  }
  
  InitializeParameters();
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //----------------------------------------------     
  // Define Materials
  //----------------------------------------------
  m_matBox       = nist->FindOrBuildMaterial("G4_Al");
  m_matQuartz    = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  m_matOil       = nist->FindOrBuildMaterial("G4_WATER");
  m_matReflector = nist->FindOrBuildMaterial("G4_Au"); 
  m_matAbsorber  = nist->FindOrBuildMaterial("G4_W");
  m_matTop       = nist->FindOrBuildMaterial("G4_Fe");

  
  //----------------------------------------------     
  // Define Material Properties
  //----------------------------------------------
  const G4int NUMENTRIES = 2;
  
  G4double ephoton         [NUMENTRIES] = {2.00*eV,4.80*eV};

  G4double rindexOil       [NUMENTRIES] = {1.48,1.48};
  G4double absorptionOil   [NUMENTRIES] = {26*m,26*m};

  G4double rindexQuartz    [NUMENTRIES] = {1.46,1.46};
  G4double absorptionQuartz[NUMENTRIES] = {46*m,46*m};

  //Fill in the Marterial properties table for each material.
  //Guide for undestanding Optical processes at http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Photo
  G4MaterialPropertiesTable *oilMPT = new G4MaterialPropertiesTable();
  oilMPT->AddProperty("RINDEX",ephoton,rindexOil,NUMENTRIES);
  oilMPT->AddProperty("ABSLENGTH",ephoton,absorptionOil,NUMENTRIES);
  m_matOil->SetMaterialPropertiesTable(oilMPT);

  G4MaterialPropertiesTable *quartzMPT = new G4MaterialPropertiesTable();
  quartzMPT->AddProperty("RINDEX",ephoton,rindexQuartz,NUMENTRIES);
  quartzMPT->AddProperty("ABSLENGTH",ephoton,absorptionQuartz,NUMENTRIES);   
  m_matQuartz->SetMaterialPropertiesTable(quartzMPT);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction :: DefineBorderProperties(){
  const G4int NUMENTRIES = 2;
  G4double ephoton[NUMENTRIES] = {2.00*eV,4.80*eV};

  //----------------------------------------------     
  // Box Skin
  //----------------------------------------------
  G4double boxReflectivity      [NUMENTRIES] = {0.4, 0.4};
  G4double boxEfficiency        [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* boxMPT
    = new G4MaterialPropertiesTable();
  boxMPT->AddProperty("REFLECTIVITY", ephoton, boxReflectivity, NUMENTRIES);
  boxMPT->AddProperty("EFFICIENCY"  , ephoton, boxEfficiency  , NUMENTRIES);
  G4OpticalSurface* boxOS =
    new G4OpticalSurface("BoxOpSurface",unified, polished, dielectric_metal);
  boxOS->SetMaterialPropertiesTable( boxMPT );

  new G4LogicalSkinSurface("boxSkinSurface", m_logicBox, boxOS );
  
  //----------------------------------------------     
  // Reflector Skin
  //----------------------------------------------
  G4double reflectorReflectivity[NUMENTRIES] = {0.95, 0.95};
  G4double reflectorEfficiency  [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* reflectorMPT
    = new G4MaterialPropertiesTable();
  reflectorMPT->AddProperty("REFLECTIVITY", ephoton, reflectorReflectivity, NUMENTRIES);
  reflectorMPT->AddProperty("EFFICIENCY"  , ephoton, reflectorEfficiency  , NUMENTRIES);
  G4OpticalSurface* reflectorOS =
    new G4OpticalSurface("ReflectorOpSurface",unified, polished, dielectric_metal);
  reflectorOS->SetMaterialPropertiesTable( reflectorMPT );

  for( auto& logicReflector : m_v_logicPanel )
    new G4LogicalSkinSurface("reflectorSkinSurface", logicReflector, reflectorOS );
    
  //----------------------------------------------     
  // Absorber Skin
  //----------------------------------------------
  G4double absorberReflectivity      [NUMENTRIES] = {0.4, 0.4};
  G4double absorberEfficiency        [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* absorberMPT
    = new G4MaterialPropertiesTable();
  absorberMPT->AddProperty("REFLECTIVITY", ephoton, absorberReflectivity, NUMENTRIES);
  absorberMPT->AddProperty("EFFICIENCY"  , ephoton, absorberEfficiency  , NUMENTRIES);
  G4OpticalSurface* absorberOS =
    new G4OpticalSurface("AbsorberOpSurface",unified, polished, dielectric_metal);
  absorberOS->SetMaterialPropertiesTable( absorberMPT );

  for( auto& logicAbsorber : m_v_logicAbsorber )
    new G4LogicalSkinSurface("absorberSkinSurface", logicAbsorber, absorberOS );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  //----------------------------------------------     
  // World
  //----------------------------------------------
  G4double world_sizeX  = 1.2 * m_moduleSizeX;   // mm
  G4double world_sizeY  = 1.2 * m_moduleSizeY;   // mm
  G4double world_sizeZ  = 1.2 * m_moduleSizeZ;   // mm

  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");
  
  printf( "Building world with %5.1f x %5.1f x %5.1f\n",
	  world_sizeX, world_sizeY, world_sizeZ );
  
  m_solidWorld =    
    new G4Box("World",                       //its name
	      0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);   //its size
  
  m_logicWorld =                         
    new G4LogicalVolume(m_solidWorld,        //its solid
                        g4Air,               //its material
                        "World");            //its name
                                   
  m_physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      m_logicWorld,          //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      m_checkOverlaps);      //overlaps checking
                     
  //----------------------------------------------     
  // Box
  //---------------------------------------------- 
  m_solidBox =    
    new G4Box("Box",                         //its name
	      0.5*m_moduleSizeX, 0.5*m_moduleSizeY, 0.5*m_moduleSizeZ);  //its size
  
  m_logicBox =                         
    new G4LogicalVolume(m_solidBox,          //its solid
                        m_matBox,            //its material
                        "Box");              //its name
                                   
  m_physBox = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      m_logicBox,            //its logical volume
                      "Box",                 //its name
                      m_logicWorld,          //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      m_checkOverlaps);      //overlaps checking
                     
  //----------------------------------------------     
  // Quartz
  //----------------------------------------------
  G4double quartzThickness     = config->GetValue( "quartzThickness", 2);
  G4double moduleWallThickness = config->GetValue( "moduleWallThickness", 2);

  std::vector< std::pair<G4double, G4double> > quartzCenters;
  quartzCenters.emplace_back( 0,  0.5*m_moduleSizeY -
			      moduleWallThickness - 0.5*quartzThickness ); // top
  quartzCenters.emplace_back( 0, -0.5*m_moduleSizeY +
			      moduleWallThickness + 0.5*quartzThickness ); // bottom
  quartzCenters.emplace_back( -0.5*m_moduleSizeX + moduleWallThickness +
			      0.5*quartzThickness, 0 ); // left
  quartzCenters.emplace_back(  0.5*m_moduleSizeX - moduleWallThickness -
			       0.5*quartzThickness, 0 ); // right
  
  std::vector< std::pair<G4double, G4double> > quartzSizeXY;
  quartzSizeXY.emplace_back( m_moduleSizeX - 2*moduleWallThickness - 2*quartzThickness,
			     quartzThickness );          // top
  quartzSizeXY.emplace_back( m_moduleSizeX - 2*moduleWallThickness - 2*quartzThickness,
			     quartzThickness );          // bottom
  quartzSizeXY.emplace_back( quartzThickness,
			     m_moduleSizeY - 2*moduleWallThickness); // left
  quartzSizeXY.emplace_back( quartzThickness,
			     m_moduleSizeY - 2*moduleWallThickness); // right
  
  G4double quartzSizeZ = m_moduleSizeZ - 2*moduleWallThickness;
  
  G4VisAttributes* quartzColor = new G4VisAttributes( G4Colour::Blue() );
  
  for( unsigned int qn = 0; qn < quartzCenters.size(); qn++ ){
    m_v_solidQuartz.
      push_back( new G4Box("Quartz",            //its name
			   0.5*quartzSizeXY.at(qn).first,
			   0.5*quartzSizeXY.at(qn).second,
			   0.5*quartzSizeZ) );  //its size
    
    m_v_logicQuartz.
      push_back( new G4LogicalVolume(m_v_solidQuartz.back(), //its solid
				     m_matQuartz,            //its material
				     "Quartz") );            //its name

    m_v_logicQuartz.back()->SetVisAttributes( quartzColor );
    
    m_v_physQuartz.
      push_back( new G4PVPlacement(0,                       //no rotation
				   G4ThreeVector(quartzCenters.at(qn).first,
						 quartzCenters.at(qn).second,
						 0),
				   m_v_logicQuartz.back(),  //its logical volume
				   "Quartz",                //its name
				   m_logicBox,              //its mother  volume
				   false,                   //no boolean operation
				   qn,                      //copy number
				   m_checkOverlaps) );      //overlaps checking
  }

  //----------------------------------------------     
  // Top
  //----------------------------------------------
  G4double absorberSizeY = m_chamberSizeY;        
  
  G4double topSizeX      = m_moduleSizeX - 2*(moduleWallThickness + quartzThickness);
  G4double topSizeY      = m_moduleSizeY - absorberSizeY - 2*(moduleWallThickness + quartzThickness);
  G4double topSizeZ      = m_moduleSizeZ - 2*moduleWallThickness;
  G4double topPosY       = 0.5*m_moduleSizeY - moduleWallThickness -
    quartzThickness - 0.5*topSizeY;
  
  m_solidTop =    
    new G4Box("Top",                         //its name
	      0.5*topSizeX, 0.5*topSizeY, 0.5*topSizeZ);  // its size
  
  m_logicTop =                         
    new G4LogicalVolume(m_solidTop,          //its solid
                        m_matTop,            //its material
                        "Top");              //its name
  
  m_physTop = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,topPosY,0),       //at (0,topPosY,0)
                      m_logicTop,            //its logical volume
                      "Top",                 //its name
                      m_logicBox,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      m_checkOverlaps);      //overlaps checking
  
  G4VisAttributes* topColor = new G4VisAttributes( G4Colour::Grey() );
  m_logicTop->SetVisAttributes( topColor );
  
  //----------------------------------------------     
  // Oil
  //----------------------------------------------
  G4double oilSizeX      = m_chamberSizeX;
  G4double oilSizeY      = m_chamberSizeY;
  G4double oilSizeZ      = m_chamberSizeZ;
  G4double oilPosY       = -0.5*m_moduleSizeY + moduleWallThickness +
    quartzThickness + 0.5*oilSizeY;
  
  m_solidOil =    
    new G4Box("Oil",                         //its name
	      0.5*oilSizeX, 0.5*oilSizeY, 0.5*oilSizeZ);  // its size
  
  m_logicOil =                         
    new G4LogicalVolume(m_solidOil,          //its solid
                        m_matOil,            //its material
                        "Oil");              //its name
  
  m_physOil = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, oilPosY, 0), //at (0,oilPosY,0)
                      m_logicOil,            //its logical volume
                      "Oil",                 //its name
                      m_logicBox,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      m_checkOverlaps);      //overlaps checking
  
  G4VisAttributes* oilColor = new G4VisAttributes( G4Colour::Magenta() );
  m_logicOil->SetVisAttributes( oilColor );

  //----------------------------------------------     
  // Build Plates with Reflectors 
  //----------------------------------------------
  std::string geoConfig = config->GetValue("geometryConfiguration","chevron");
  if( !geoConfig.compare("chevron") ) BuildChevronGeo();
  else if( !geoConfig.compare("diagonal") ) BuildDiagonalGeo();

  //----------------------------------------------     
  // Define Surface/Border Properties
  //----------------------------------------------  
  DefineBorderProperties();
  
  //----------------------------------------------     
  // SD and Scoring Volumes
  //----------------------------------------------  
  // Set Box as scoring volume
  //
  m_scoringVolume = m_logicBox;
  
  //
  //always return the physical World
  //
  return m_physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction :: BuildChevronGeo(){
  // Get config
  TEnv* config = m_sd->GetConfig();
  
  double absorberThicknessZ = config->GetValue( "absorberThicknessZ", 10 );
  double      gapThicknessZ = config->GetValue( "gapThicknessZ", 2 );
  double              theta = config->GetValue( "absorberTheta", 0.5 );

  double reflectorThickness  = config->GetValue( "reflectorThickness", 0.1 );
  double reflectorThicknessZ = reflectorThickness / cos( theta );

  double tanTheta   = TMath::Tan( theta );

  double totalThicknessZ = absorberThicknessZ + 2 * reflectorThicknessZ ;

  double yStep = totalThicknessZ + gapThicknessZ;
  
  double z0 = 0.5*m_chamberSizeZ - 0.5*m_chamberSizeX * tanTheta - totalThicknessZ;

  printf(" 0.5 m_chamberSizeZ = %f\n", 0.5 * m_chamberSizeZ );
  printf(" 0.5 m_chamberSizeX  = %f\n", 0.5 * m_chamberSizeX  );
  printf(" totalThick = %f\n", totalThicknessZ );

  std::vector< std::pair< G4double, G4double > > panelCenters;
  std::vector< G4double > panelLengths;
  
  while( z0 > -0.5 * m_chamberSizeZ ){
    printf("\n-----Panel-----\n");
    printf("zo=%5.1f ", z0 );
    
    double zCoord =
      z0 + 0.5*totalThicknessZ + 0.25 * m_chamberSizeX * tan(theta);

    double xCoord = 0.25 * m_chamberSizeX;

    panelCenters.emplace_back( xCoord, zCoord);

    printf( "(x,y) = (%5.1f,%5.1f)\n",
	    panelCenters.back().first, panelCenters.back().second );
    
    z0 -= yStep; 
  }

  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  G4double absorberSizeY       = m_chamberSizeY;
  G4double absorberSizeZ       = absorberThicknessZ * TMath::Cos( theta );
  
  // These are before rotation. Thats why the x<->z
  G4double panelSizeX         = absorberSizeZ + 2*reflectorThicknessZ;
  G4double panelSizeY         = absorberSizeY;
  G4double panelSizeZ         = 0.5 * m_chamberSizeX;

  // These are before rotation. Thats why the x<->z
  G4double absorberSizeXX     = absorberSizeZ;
  G4double absorberSizeYY     = absorberSizeY;
  G4double absorberSizeZZ     = 0.5 * m_chamberSizeX;
    
  G4VisAttributes* panelColor    = new G4VisAttributes( G4Colour::Yellow() );
  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );

  G4RotationMatrix* rotPanel = new G4RotationMatrix();
  rotPanel->rotateY( 90*deg );
  
  for( unsigned int sn = 0; sn < panelCenters.size(); sn++ ){
    
    for( int LR = 0; LR < 2; LR++ ){

      m_v_solidPanel. 
	push_back( new G4Para("Panel",               //its name
			      0.5 * panelSizeX,
			      0.5 * panelSizeY,
			      0.5 * panelSizeZ,
			      0,
			      TMath::Power( -1, LR + 1 ) * theta,
			      0 ) );     //its size
    
      m_v_solidAbsorber.
	push_back( new G4Para("Absorber",            //its name
			      0.5*absorberSizeXX,
			      0.5*absorberSizeYY,
			      0.5*absorberSizeZZ,
			      0,
			      TMath::Power( -1, LR + 1 ) * theta,
			      0) );  
    
      m_v_logicPanel.
	push_back( new G4LogicalVolume(m_v_solidPanel.back(), //its solid
				       m_matReflector,        //its material
				       "Panel") );            //its name

      m_v_logicPanel.back()->SetVisAttributes( panelColor );

      m_v_logicAbsorber.
	push_back( new G4LogicalVolume(m_v_solidAbsorber.back(), //its solid
				       m_matAbsorber,            //its material
				       "Absorber") );            //its name

      m_v_logicAbsorber.back()->SetVisAttributes( absorberColor );
    
      m_v_physPanel.
	push_back( new G4PVPlacement(rotPanel,               //no rotation
				     G4ThreeVector( TMath::Power( -1, LR + 1 ) *
						    panelCenters.at(sn).first,
						    0,
						    panelCenters.at(sn).second),
				     m_v_logicPanel.back(),  //its logical volume
				     "Panel",                //its name
				     m_logicOil,             //its mother  volume
				     false,                  //no boolean operation
				     sn,                     //copy number
				     m_checkOverlaps) );     //overlaps checking

      m_v_physAbsorber.
	push_back( new G4PVPlacement(0,                         //no rotation
				     G4ThreeVector(),
				     m_v_logicAbsorber.back(),  //its logical volume
				     "Absorber",                //its name
				     m_v_logicPanel.back(),     //its mother  volume
				     false,                     //no boolean operation
				     sn,                        //copy number
				     m_checkOverlaps) );        //overlaps checking
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction :: BuildDiagonalGeo()
{
  // Get config
  TEnv* config = m_sd->GetConfig();
 
  double absorberThicknessZ = config->GetValue( "absorberThicknessZ", 5 );
  double      gapThicknessZ = config->GetValue( "gapThicknessZ", 10 );
  double              theta = config->GetValue( "absorberTheta", 0.5 );
  double             deltaX = config->GetValue( "absorberDeltaX", 10 );
  
  double reflectorThickness  = config->GetValue( "reflectorThickness", 0.1 );
  double reflectorThicknessZ = reflectorThickness / cos( theta );

  double thickness = absorberThicknessZ * cos(theta) + 2*reflectorThicknessZ;
  double       gap =      gapThicknessZ * cos(theta);

  double tanTheta   = TMath::Tan( theta );
  double m_chamberSizeXPrime = m_chamberSizeX/TMath::Cos( theta ) - thickness * tanTheta; 
  printf( "m_chamberSizeXPrime = %f\n", m_chamberSizeXPrime );
  
  double yStep = ( thickness + gap ) / TMath::Cos( theta ); 
  printf( "yStep = %f\n", yStep );
  
  double dz    = thickness * TMath::Cos( theta );
  double dx    = thickness * TMath::Sin( theta );
  printf( "dx = %f      dz = %f\n", dx, dz );
  
  double z0    = 0.5 * m_chamberSizeZ - 0.5 * tanTheta * m_chamberSizeX - deltaX;
  printf( "z0 = %f\n", z0 );

  std::vector< std::pair< G4double, G4double > > panelCenters;
  std::vector< G4double > panelLengths; 
  
  while( z0 > -0.5 * m_chamberSizeZ ){
    printf("\n-----Panel-----\n");
    double x0 = -0.5 * m_chamberSizeX;
    
    double z1 = z0 + m_chamberSizeXPrime * TMath::Sin( theta );
    double x1 = 0.5 * m_chamberSizeX - dx;
    
    double z2 = z1 - dz;
    double x2 = 0.5 * m_chamberSizeX;

    double z3 = z0 - dz;
    double x3 = - 0.5 * m_chamberSizeX + dx;

    // test bad cases
    if( ( z1 > 0.5 * m_chamberSizeZ ) && ( z2 >= 0.5 * m_chamberSizeZ ) ){
      z1 = m_chamberSizeZ / 2;
      x1 = ( z1 - z0 ) / TMath::Tan( theta ) - m_chamberSizeX / 2;

      z2 = z1 - dz;
      x2 = x1 + dx;
    }
    else if( ( z1 >= 0.5 * m_chamberSizeZ ) && ( z2 < 0.5 * m_chamberSizeZ ) ){
      z1  = m_chamberSizeZ / 2;
      x1  = (m_chamberSizeZ/2 - z0)/TMath::Tan( theta ) - m_chamberSizeX/2;
      
      z2  = z1 - dz;
      x2  = m_chamberSizeX / 2;
    }
    else if( ( z0 >= -0.5 * m_chamberSizeZ ) && ( z3 < -0.5 * m_chamberSizeZ ) ){
      break;
    }

    double zPts[5];  double xPts[5];
    
    zPts[0] = z0;  xPts[0] = x0;
    zPts[1] = z1;  xPts[1] = x1;
    zPts[2] = z2;  xPts[2] = x2;
    zPts[3] = z3;  xPts[3] = x3;

    for ( int i = 0; i < 5; i++ ) {
      printf("( %f, %f )\n", xPts[i], zPts[i] );
    }


    panelCenters.emplace_back( 0.5*(x1 + x0 + dx), 0.5*(z1 + z0 - dz) );
    // take some small fraction off due to some rounding problems
    // that create overlap later
    panelLengths.push_back( 0.999 * TMath::Sqrt( TMath::Power( x1 - x0, 2 ) +
						 TMath::Power( z1 - z0, 2 ) ) );

    printf( "(x,y) = (%5.1f,%5.1f)    length = %5.1f\n",
	    panelCenters.back().first, panelCenters.back().second, panelLengths.back() );
    
    // decrement z0
    z0 -= yStep;
  }

  //----------------------------------------------     
  // Plates
  //---------------------------------------------- 
  G4double absorberSizeY       = m_chamberSizeY;
  G4double absorberSizeZ       = absorberThicknessZ * TMath::Cos( theta );
 
  G4double panelSizeY         = absorberSizeY;
  G4double panelSizeZ         = absorberSizeZ + 2*reflectorThickness;
    
  G4VisAttributes* panelColor    = new G4VisAttributes( G4Colour::Yellow() );
  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );

  G4RotationMatrix* rotPanel = new G4RotationMatrix();
  rotPanel->rotateY( theta*rad);
  
  for( unsigned int sn = 0; sn < panelCenters.size(); sn++ ){
    m_v_solidPanel. 
      push_back( new G4Box("Panel",               //its name
			   0.5*panelLengths.at(sn),
			   0.5*panelSizeY,  
			   0.5*panelSizeZ) );     //its size
    
    m_v_solidAbsorber.
      push_back( new G4Box("Absorber",            //its name
			   0.5*panelLengths.at(sn),
			   0.5*absorberSizeY,
			   0.5*absorberSizeZ) );  //its size

    m_v_logicPanel.
      push_back( new G4LogicalVolume(m_v_solidPanel.back(), //its solid
				     m_matReflector,        //its material
				     "Panel") );            //its name

    m_v_logicPanel.back()->SetVisAttributes( panelColor );

    m_v_logicAbsorber.
      push_back( new G4LogicalVolume(m_v_solidAbsorber.back(), //its solid
				     m_matAbsorber,            //its material
				     "Absorber") );            //its name

    m_v_logicAbsorber.back()->SetVisAttributes( absorberColor );
    
    m_v_physPanel.
      push_back( new G4PVPlacement(rotPanel,               //no rotation
				   G4ThreeVector(panelCenters.at(sn).first,
						 0,
						 panelCenters.at(sn).second),
				   m_v_logicPanel.back(),  //its logical volume
				   "Panel",                //its name
				   m_logicOil,             //its mother  volume
				   false,                  //no boolean operation
				   sn,                     //copy number
				   m_checkOverlaps) );     //overlaps checking

    m_v_physAbsorber.
      push_back( new G4PVPlacement(0,                         //no rotation
				   G4ThreeVector(),
				   m_v_logicAbsorber.back(),  //its logical volume
				   "Absorber",                //its name
				   m_v_logicPanel.back(),     //its mother  volume
				   false,                     //no boolean operation
				   sn,                        //copy number
				   m_checkOverlaps) );        //overlaps checking
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
