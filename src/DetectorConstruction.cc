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
  : G4VUserDetectorConstruction(), m_matHousing(0), m_matQuartz(0),
    m_matOil(0), m_matReflector(0), m_matAbsorber(0),
    m_solidWorld(0), m_logicWorld(0), m_physWorld(0),
    m_solidOil(0), m_logicOil(0), m_physOil(0),
    m_scoringVolume(0),
    m_sd(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( SharedData* sd )
  : G4VUserDetectorConstruction(), m_matHousing(0), m_matQuartz(0),
    m_matOil(0), m_matReflector(0), m_matAbsorber(0),
    m_solidWorld(0), m_logicWorld(0), m_physWorld(0),
    m_solidOil(0), m_logicOil(0), m_physOil(0),
    m_scoringVolume(0),
    m_sd( sd )
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::InitializeParameters()
{}

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
  m_matHousing   = nist->FindOrBuildMaterial("G4_Al");
  m_matQuartz    = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  m_matOil       = nist->FindOrBuildMaterial("G4_WATER");
  m_matReflector = nist->FindOrBuildMaterial("G4_Au"); 
  m_matAbsorber  = nist->FindOrBuildMaterial("G4_W");
  
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
  // Housing Skin
  //----------------------------------------------
  G4double housingReflectivity      [NUMENTRIES] = {0.4, 0.4};
  G4double housingEfficiency        [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* housingMPT
    = new G4MaterialPropertiesTable();
  housingMPT->AddProperty("REFLECTIVITY", ephoton, housingReflectivity, NUMENTRIES);
  housingMPT->AddProperty("EFFICIENCY"  , ephoton, housingEfficiency  , NUMENTRIES);
  G4OpticalSurface* housingOS =
    new G4OpticalSurface("HousingOpSurface",unified, polished, dielectric_metal);
  housingOS->SetMaterialPropertiesTable( housingMPT );

  new G4LogicalSkinSurface("housingSkinSurface", m_logicHousing, housingOS );
  
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
  // Set Some Values
  //----------------------------------------------
  G4double theta               = config->GetValue( "absorberTheta", 0.5 );
  
  G4double moduleSizeX         = config->GetValue( "moduleSizeX", 90);
  G4double moduleSizeY         = config->GetValue( "moduleSizeY", 180);
  G4double moduleSizeZ         = config->GetValue( "moduleSizeZ", 150);

  G4double housingThickness    = config->GetValue( "housingThickness", 2);
    
  G4double absorberThicknessZ  = config->GetValue( "absorberThicknessZ", 10 );
  G4double gapThicknessZ       = config->GetValue( "gapThicknessZ", 2 );

  G4double reflectorThicknessZ = config->GetValue( "reflectorThicknessZ", 0.1 );
  G4double totalThicknessZ     = absorberThicknessZ + 2 * reflectorThicknessZ ;

  // Option to switch on/off checking of volumes overlaps
  //
  bool checkOverlaps = true;
  
  //----------------------------------------------     
  // World
  //----------------------------------------------
  G4double worldSizeX       = 1.2 * moduleSizeX;    // mm
  G4double worldSizeY       = 1.1 * moduleSizeY;    // mm
  G4double worldSizeZ       =
    1.1 * ( moduleSizeZ + 0.5 * moduleSizeX * TMath::Tan(theta) ); // mm
    

  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");
  
  printf( "Building world with %5.1f x %5.1f x %5.1f\n",
	  worldSizeX, worldSizeY, worldSizeZ );
  
  m_solidWorld =    
    new G4Box("World",              //its name
	      0.5*worldSizeX,       //its size
	      0.5*worldSizeY,
	      0.5*worldSizeZ );   
  
  m_logicWorld =                         
    new G4LogicalVolume(m_solidWorld,      //its solid
                        g4Air,             //its material
                        "World");          //its name
                                   
  m_physWorld = 
    new G4PVPlacement(0,                   //no rotation
                      G4ThreeVector(),     //at (0,0,0)
                      m_logicWorld,        //its logical volume
                      "World",             //its name
                      0,                   //its mother  volume
                      false,               //no boolean operation
                      0,                   //copy number
                      checkOverlaps);      //overlaps checking

  //----------------------------------------------     
  // variables ending in prime are
  // for construction of G4Para
  // the sizes there have to be adjusted for
  // depending on theta
  //----------------------------------------------
  
  //----------------------------------------------     
  // Housing
  //----------------------------------------------  
  G4double moduleHalfSizeX       = 0.5 * moduleSizeX;
  G4double moduleHalfSizeXprime  = moduleHalfSizeX / TMath::Cos(theta);

  G4double moduleSizeZprime  = moduleSizeZ * TMath::Cos(theta);
  
  G4double housingSizeX      = moduleHalfSizeXprime;
  G4double housingSizeY      = moduleSizeY;
  G4double housingSizeZ      = moduleSizeZprime;

  printf("\n");  
  printf("housingSizeX = %f\n", housingSizeX);
  printf("housingSizeY = %f\n", housingSizeY);
  printf("housingSizeZ = %f\n", housingSizeZ);
  
  m_solidHousing =    
    new G4Para("Housing",               //its name
	       0.5 * housingSizeX,      //its size
	       0.5 * housingSizeY,
	       0.5 * housingSizeZ,
	       0, theta, 0 );              
  
  m_logicHousing =                         
    new G4LogicalVolume(m_solidHousing,      //its solid
                        m_matHousing,        //its material
                        "Housing");          //its name

  G4RotationMatrix* rotModuleL = new G4RotationMatrix();
  rotModuleL->rotateY( theta*rad );

  G4RotationMatrix* rotModuleR = new G4RotationMatrix();
  rotModuleR->rotateZ( 180*deg );
  rotModuleR->rotateY( theta*rad );

  G4double housingPosX = 0.5 * moduleHalfSizeX;
  
  m_physHousingL = 
    new G4PVPlacement(rotModuleL, 
                      G4ThreeVector( -housingPosX, 0, 0),
                      m_logicHousing,      //its logical volume
                      "Housing",           //its name
                      m_logicWorld,        //its mother  volume
                      false,               //no boolean operation
                      0,                   //copy number
                      checkOverlaps);      //overlaps checking
  m_physHousingR = 
    new G4PVPlacement(rotModuleR,     
                      G4ThreeVector( +housingPosX, 0, 0),
                      m_logicHousing,      //its logical volume
                      "Housing",           //its name
                      m_logicWorld,        //its mother  volume
                      false,               //no boolean operation
                      1,                   //copy number
                      checkOverlaps);      //overlaps checking
  //----------------------------------------------     
  // Oil
  //----------------------------------------------
  G4double chamberHalfSizeX      = moduleHalfSizeX - housingThickness;
  G4double chamberHalfSizeXprime = chamberHalfSizeX / TMath::Cos(theta);
  
  G4double chamberSizeY        = moduleSizeY - 2 * housingThickness;
  G4double chamberSizeZprime   = moduleSizeZprime - 2 * housingThickness;

  G4double oilSizeX      = chamberHalfSizeXprime;
  G4double oilSizeY      = chamberSizeY;
  G4double oilSizeZ      = chamberSizeZprime;

  printf("\n");  
  printf("oilSizeX = %f\n", oilSizeX);
  printf("oilSizeY = %f\n", oilSizeY);
  printf("oilSizeZ = %f\n", oilSizeZ);

  
  m_solidOil =    
    new G4Para("Oil",             //its name
	       0.5 * oilSizeX,      // its size
	       0.5 * oilSizeY,
	       0.5 * oilSizeZ,
	       0, theta, 0);  
  
  m_logicOil =                         
    new G4LogicalVolume(m_solidOil,          //its solid
                        m_matOil,            //its material
                        "Oil");              //its name

  double oilPosX = 0.5 * housingThickness / TMath::Cos(theta);
  
  m_physOil = 
    new G4PVPlacement(0,                      //no rotation
                      G4ThreeVector( oilPosX, 0, 0 ), 
                      m_logicOil,             //its logical volume
                      "Oil",                  //its name
                      m_logicHousing,         //its mother  volume
                      false,                  //no boolean operation
                      0,                      //copy number
                      checkOverlaps);       //overlaps checking
  
  G4VisAttributes* oilColor = new G4VisAttributes( G4Colour::Magenta() );
  m_logicOil->SetVisAttributes( oilColor );

  //----------------------------------------------     
  // Build Plates with Reflectors 
  //----------------------------------------------
  double dZprimePanel     =  0.5 * totalThicknessZ * TMath::Cos(theta);
  double zCoordPrimePanel =  0.5 * chamberSizeZprime - dZprimePanel;  
  double zStepPrimePanel  =  (totalThicknessZ + gapThicknessZ) * TMath::Cos(theta);

  printf(" 0.5 chamberSizeZprime = %f\n", 0.5 * chamberSizeZprime );
  printf(" 0.5 chamberSizeXprime = %f\n", chamberHalfSizeXprime  );
  printf(" totalThicknessZ  = %f\n", totalThicknessZ );
  printf(" zStepPrimePanel  = %f\n", zStepPrimePanel );
  printf(" zCoordPrimePanel = %f\n", zCoordPrimePanel );
  printf(" theta            = %f\n", theta );
  printf(" tanTheta         = %f\n", TMath::Tan( theta ) );
  
  std::vector< std::pair<G4double, G4double> > panelCoords;
 
  while( zCoordPrimePanel > -0.5 * chamberSizeZprime + dZprimePanel){
    printf("\n-----Panel-----\n");
    printf("zCoordPrimePanel=%5.1f ", zCoordPrimePanel );

    double xCoord = zCoordPrimePanel * TMath::Tan( theta );;
    
    panelCoords.emplace_back( xCoord, zCoordPrimePanel);

    printf( "(x,z) = (%5.1f,%5.1f)\n",
	    panelCoords.back().first,
	    panelCoords.back().second );
    
    zCoordPrimePanel -= zStepPrimePanel; 
  }

  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  G4double absorberSizeX  = chamberHalfSizeXprime;
  G4double absorberSizeY  = chamberSizeY;
  G4double absorberSizeZ  = absorberThicknessZ * TMath::Cos(theta);
  
  G4double panelSizeX     = chamberHalfSizeXprime;
  G4double panelSizeY     = absorberSizeY;
  G4double panelSizeZ     = totalThicknessZ * TMath::Cos(theta);
    
  G4VisAttributes* panelColor    = new G4VisAttributes( G4Colour::Yellow() );
  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );
  
  for( unsigned int cn = 0; cn < panelCoords.size(); cn++ ){

    m_v_solidPanel. 
      push_back( new G4Para("Panel",               //its name
			    0.5 * panelSizeX,
			    0.5 * panelSizeY,
			    0.5 * panelSizeZ,
			    0, theta, 0 ) );       //its size
    
    m_v_solidAbsorber.
      push_back( new G4Para("Absorber",            //its name
			    0.5*absorberSizeX,
			    0.5*absorberSizeY,
			    0.5*absorberSizeZ,
			    0, theta, 0) );  
    
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
      push_back( new G4PVPlacement(0,
				   G4ThreeVector( panelCoords.at(cn).first,
						  0,
						  panelCoords.at(cn).second ),
				   m_v_logicPanel.back(),  //its logical volume
				   "Panel",                //its name
				   m_logicOil,             //its mother  volume
				   false,                  //no boolean operation
				   cn,                     //copy number
				   checkOverlaps) );       //overlaps checking

    m_v_physAbsorber.
      push_back( new G4PVPlacement(0,                         //no rotation
				   G4ThreeVector(),
				   m_v_logicAbsorber.back(),  //its logical volume
				   "Absorber",                //its name
				   m_v_logicPanel.back(),     //its mother  volume
				   false,                     //no boolean operation
				   cn,                        //copy number
				   checkOverlaps) );          //overlaps checking
   
  }

  //----------------------------------------------     
  // Quartz Coordinates
  //----------------------------------------------
  G4double quartzThickness      = config->GetValue("quartzThickness", 2);
  G4double quartzThicknessPrime = quartzThickness / TMath::Cos(theta);
  
  G4double dZprimeQuartz     =  0.5 * gapThicknessZ * TMath::Cos(theta);
  G4double zCoordPrimeQuartz =  0.5 * chamberSizeZprime - 2*dZprimePanel - dZprimeQuartz;  
  G4double zStepPrimeQuartz  =  (totalThicknessZ + gapThicknessZ) * TMath::Cos(theta);
  
  std::vector< std::pair<G4double, G4double> > quartzCoords;
 
  while( zCoordPrimeQuartz > -0.5 * chamberSizeZprime + dZprimeQuartz){
    printf("\n-----Quartz-----\n");
    printf("zCoordPrimeQuartz=%5.1f ", zCoordPrimeQuartz );

    G4double xCoord =
      zCoordPrimeQuartz * TMath::Tan( theta ) -
      0.5*chamberHalfSizeXprime + 0.5 * quartzThicknessPrime;
    
    quartzCoords.emplace_back( xCoord, zCoordPrimeQuartz);

    printf( "(x,z) = (%5.1f,%5.1f)\n",
	    quartzCoords.back().first,
	    quartzCoords.back().second );
    
    zCoordPrimeQuartz -= zStepPrimeQuartz; 
  }

  //----------------------------------------------     
  // Quartz Placement
  //----------------------------------------------
  G4double quartzSizeX     = quartzThicknessPrime;
  G4double quartzSizeY     = chamberSizeY;
  G4double quartzSizeZ     = gapThicknessZ * TMath::Cos(theta);
    
  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Blue() );

  for( unsigned int cn = 0; cn < quartzCoords.size(); cn++ ){

    m_v_solidQuartz. 
      push_back( new G4Para("Quartz",             //its name
			    0.5 * quartzSizeX,    //its size
			    0.5 * quartzSizeY,
			    0.5 * quartzSizeZ,
			    0, theta, 0 ) );       
    
    m_v_logicQuartz.
      push_back( new G4LogicalVolume(m_v_solidQuartz.back(), //its solid
				     m_matReflector,         //its material
				     "Quartz") );            //its name

    m_v_logicQuartz.back()->SetVisAttributes( quartzColor );

    m_v_physQuartz.
      push_back( new G4PVPlacement(0,
				   G4ThreeVector( quartzCoords.at(cn).first,
						  0,
						  quartzCoords.at(cn).second ),
				   m_v_logicQuartz.back(),  //its logical volume
				   "Quartz",                //its name
				   m_logicOil,              //its mother  volume
				   false,                   //no boolean operation
				   cn,                      //copy number
				   checkOverlaps) );        //overlaps checking
  }
  
  
  //----------------------------------------------     
  // Define Surface/Border Properties
  //----------------------------------------------  
  DefineBorderProperties();
  
  //----------------------------------------------     
  // SD and Scoring Volumes
  //----------------------------------------------  
  // Set Box as scoring volume
  //
  m_scoringVolume = m_logicHousing;
  
  //
  //always return the physical World
  //
  return m_physWorld;
}
