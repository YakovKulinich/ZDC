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
// $Id: EMCal.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file EMCal.cc
/// \brief Implementation of the EMCal class

#include "EMCal.hh"
#include "QuartzSD.hh"
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
#include "G4SDManager.hh"

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

EMCal::EMCal( const std::string& name, const int cn,
	      const G4RotationMatrix* rot, const G4ThreeVector& pos,
	      G4LogicalVolume* mother, SharedData* sd)
  : m_name( name ), m_copyNumber( cn ), m_rot(rot), m_pos( pos ), m_logicMother( mother ),
    m_sd( sd ),
    m_matHousing(0)   , m_matQuartz(0),
    m_matEmitter(0)   , m_matReflector(0), m_matAbsorber(0),
    m_solidHousing(0) , m_logicHousing(0), m_physHousingL(0), m_physHousingR(0),
    m_solidChamber(0) , m_logicChamber(0), m_physChamber(0),
    m_solidPanel(0)   , m_logicPanel(0),
    m_solidAbsorber(0), m_logicAbsorber(0),
    m_solidQuartz(0)  , m_logicQuartz(0),
    m_solidEmitter(0) , m_logicEmitter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMCal::EMCal()
  : m_name(""), m_copyNumber( 0 ), m_rot(NULL), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_sd(NULL), 
    m_matHousing(0)   , m_matQuartz(0),
    m_matEmitter(0)   , m_matReflector(0), m_matAbsorber(0),
    m_solidHousing(0) , m_logicHousing(0), m_physHousingL(0), m_physHousingR(0),
    m_solidChamber(0) , m_logicChamber(0), m_physChamber(0),
    m_solidPanel(0)   , m_logicPanel(0),
    m_solidAbsorber(0), m_logicAbsorber(0),
    m_solidQuartz(0)  , m_logicQuartz(0),
    m_solidEmitter(0) , m_logicEmitter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMCal::~EMCal()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMCal::Construct(){
  DefineMaterials();
  ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMCal::DefineMaterials()
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
  m_matAbsorber  = nist->FindOrBuildMaterial("G4_W");

  std::string matReflectorName = config->GetValue("reflectorType","Au");
  if( matReflectorName == "Au" ){  
    m_matReflector = nist->FindOrBuildMaterial("G4_Au"); 
  } else if ( matReflectorName == "Al" ){
    m_matReflector = nist->FindOrBuildMaterial("G4_Al"); 
  } else{  // default
    m_matReflector = nist->FindOrBuildMaterial("G4_Au"); 
  }

  G4int natomsC = 16;
  G4int natomsH = 26;
  G4int ncomponents = 2;
  G4double density = 0.85*g/mL; 
  std::string matEmitterName = config->GetValue("emitterType","LAB");
  
  if( matEmitterName == "LAB" ){
    G4double a = 12.011 *g/mole;
    G4Element *C = new G4Element("Carbon","C",6,a);
    a = 1.008 * g/mole;
    G4Element *H = new G4Element("Hydrogen","H",1,a);
    m_matEmitter = new G4Material(matEmitterName,density,ncomponents);
    m_matEmitter->AddElement(C,natomsC);
    m_matEmitter->AddElement(H,natomsH);
  } else{  // default
    m_matEmitter = nist->FindOrBuildMaterial("G4_WATER"); 
  }
 
  //----------------------------------------------     
  // Define Material Properties
  //----------------------------------------------
  const G4int NUMENTRIES = 2;
  
  G4double ephoton         [NUMENTRIES] = {2.00*eV,4.80*eV};

  G4double rindexEmitter       [NUMENTRIES] = {1.48,1.48};
  G4double absorptionEmitter   [NUMENTRIES] = {26*m,26*m};

  G4double rindexQuartz    [NUMENTRIES] = {1.46,1.46};
  G4double absorptionQuartz[NUMENTRIES] = {46*m,46*m};

  //Fill in the Marterial properties table for each material.
  //Guide for undestanding Optical processes at http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Photo
  G4MaterialPropertiesTable *emitterMPT = new G4MaterialPropertiesTable();
  emitterMPT->AddProperty("RINDEX",ephoton,rindexEmitter,NUMENTRIES);
  emitterMPT->AddProperty("ABSLENGTH",ephoton,absorptionEmitter,NUMENTRIES);
  m_matEmitter->SetMaterialPropertiesTable(emitterMPT);

  G4MaterialPropertiesTable *quartzMPT = new G4MaterialPropertiesTable();
  quartzMPT->AddProperty("RINDEX",ephoton,rindexQuartz,NUMENTRIES);
  quartzMPT->AddProperty("ABSLENGTH",ephoton,absorptionQuartz,NUMENTRIES);   
  m_matQuartz->SetMaterialPropertiesTable(quartzMPT);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMCal :: DefineBorderProperties(){
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

  new G4LogicalSkinSurface("reflectorSkinSurface", m_logicPanel, reflectorOS );
    
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

  new G4LogicalSkinSurface("absorberSkinSurface", m_logicAbsorber, absorberOS );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMCal::ConstructDetector()
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

  // Variables used as input to G4Para
  G4double housingSizeX      = moduleHalfSizeXprime;
  G4double housingSizeY      = moduleSizeY;
  G4double housingSizeZ      = moduleSizeZprime;
  
  m_solidHousing =    
    new G4Para("Housing",               //its name
	       0.5 * housingSizeX,      //its size
	       0.5 * housingSizeY,
	       0.5 * housingSizeZ,
	       0, theta, 0 );              
  
  m_logicHousing =                         
    new G4LogicalVolume(m_solidHousing,    //its solid
                        m_matHousing,      //its material
                        "Housing");        //its name

  G4RotationMatrix* rotModuleL = new G4RotationMatrix();
  rotModuleL->rotateY( theta*rad );

  G4RotationMatrix* rotModuleR = new G4RotationMatrix();
  rotModuleR->rotateZ( 180*deg );
  rotModuleR->rotateY( theta*rad );

  G4double housingPosX = 0.5 * moduleHalfSizeX;

  G4ThreeVector housingPosVL =
    m_pos + G4ThreeVector( -housingPosX, 0, 0 );
  std::cout << "mpos " << m_pos << " total pos " << housingPosVL << std::endl;
  G4ThreeVector housingPosVR =
    m_pos + G4ThreeVector( +housingPosX, 0, 0 );
 
  
  m_physHousingL = 
    new G4PVPlacement(rotModuleL, 
                      housingPosVL,
                      m_logicHousing,      //its logical volume
                      "Housing",           //its name
                      m_logicMother,       //its mother  volume
                      false,               //no boolean operation
                      0 + m_copyNumber*2,  //copy number
                      checkOverlaps);      //overlaps checking
  
  m_physHousingR = 
    new G4PVPlacement(rotModuleR,     
                      housingPosVR,
                      m_logicHousing,      //its logical volume
                      "Housing",           //its name
                      m_logicMother,       //its mother  volume
                      false,               //no boolean operation
                      1 + m_copyNumber*2,  //copy number
                      checkOverlaps);      //overlaps checking
  std::cout << "DONE PLACING HOUSING" << std::endl;
  
  //----------------------------------------------     
  // Chamber
  //----------------------------------------------
  G4double chamberHalfSizeX      = moduleHalfSizeX - housingThickness;
  G4double chamberHalfSizeXprime = chamberHalfSizeX / TMath::Cos(theta);
  G4double chamberSizeZprime     = moduleSizeZprime - 2 * housingThickness;
  
  // Variables used as input to G4Para
  G4double chamberSizeX      = chamberHalfSizeXprime;
  G4double chamberSizeY      = moduleSizeY - 2 * housingThickness;
  G4double chamberSizeZ      = chamberSizeZprime;

  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");  
  
  m_solidChamber =    
    new G4Para("Chamber",            //its name
	       0.5 * chamberSizeX,   // its size
	       0.5 * chamberSizeY,
	       0.5 * chamberSizeZ,
	       0, theta, 0);  
  
  m_logicChamber =                         
    new G4LogicalVolume(m_solidChamber,    //its solid
                        g4Air,             //its material
                        "Chamber");        //its name

  double chamberPosX = 0.5 * housingThickness / TMath::Cos(theta);
  
  m_physChamber = 
    new G4PVPlacement(0,                    //no rotation
                      G4ThreeVector(chamberPosX, 0, 0 ), 
                      m_logicChamber,       //its logical volume
                      "Chamber",            //its name
                      m_logicHousing,       //its mother  volume
                      false,                //no boolean operation
                      0,                    //copy number
                      checkOverlaps);       //overlaps checking
  
  G4VisAttributes* chamberColor = new G4VisAttributes( G4Colour::Gray() );
  m_logicChamber->SetVisAttributes( chamberColor );

  //----------------------------------------------     
  // Build Plates with Reflectors 
  //----------------------------------------------
  G4double dZprimePanel     =  0.5 * totalThicknessZ * TMath::Cos(theta);
  G4double dZprimeEmitter   =  0.5 * gapThicknessZ * TMath::Cos(theta);

  G4double zCoordPrimePanel =  0.5 * chamberSizeZprime - 2*dZprimeEmitter- dZprimePanel;  
  G4double zStepPrimePanel  =  (totalThicknessZ + gapThicknessZ) * TMath::Cos(theta);

  std::vector< std::pair<G4double, G4double> > panelCoords;
 
  while( zCoordPrimePanel > -0.5 * chamberSizeZprime + dZprimePanel){
    G4double xCoord = zCoordPrimePanel * TMath::Tan( theta );;
    panelCoords.emplace_back( xCoord, zCoordPrimePanel);
    zCoordPrimePanel -= zStepPrimePanel; 
  }

  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  // Variables used as input to G4Para
  G4double absorberSizeX  = chamberHalfSizeXprime;
  G4double absorberSizeY  = chamberSizeY;
  G4double absorberSizeZ  = absorberThicknessZ * TMath::Cos(theta);

  // Variables used as input to G4Para
  G4double panelSizeX     = chamberHalfSizeXprime;
  G4double panelSizeY     = absorberSizeY;
  G4double panelSizeZ     = totalThicknessZ * TMath::Cos(theta);
    
  G4VisAttributes* panelColor    = new G4VisAttributes( G4Colour::Yellow() );
  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );


  m_solidPanel =
    new G4Para("Panel",               //its name
	       0.5 * panelSizeX,
	       0.5 * panelSizeY,
	       0.5 * panelSizeZ,
	       0, theta, 0  );        //its size
    
  m_solidAbsorber =
    new G4Para("Absorber",            //its name
	       0.5*absorberSizeX,
	       0.5*absorberSizeY,
	       0.5*absorberSizeZ,
	       0, theta, 0 );  
    
  m_logicPanel =
    new G4LogicalVolume(m_solidPanel,      //its solid
			m_matReflector,    //its material
			"Panel"  );        //its name

  m_logicPanel->SetVisAttributes( panelColor );

  m_logicAbsorber =
    new G4LogicalVolume(m_solidAbsorber,   //its solid
			m_matAbsorber,     //its material
			"Absorber" );      //its name
    
  m_logicAbsorber->SetVisAttributes( absorberColor );
    
  for( unsigned int cn = 0; cn < panelCoords.size(); cn++ ){
    m_v_physPanel.
      push_back( new G4PVPlacement(0,
				   G4ThreeVector( panelCoords.at(cn).first,
						  0,
						  panelCoords.at(cn).second ),
				   m_logicPanel,         //its logical volume
				   "Panel",              //its name
				   m_logicChamber,       //its mother  volume
				   false,                //no boolean operation
				   cn,                   //copy number
				   checkOverlaps ) );      //overlaps checking

    m_v_physAbsorber.
      push_back( new G4PVPlacement(0,                     //no rotation
				   G4ThreeVector(),
				   m_logicAbsorber,       //its logical volume
				   "Absorber",            //its name
				   m_logicPanel,          //its mother  volume
				   false,                 //no boolean operation
				   cn,                    //copy number
				   checkOverlaps) );      //overlaps checking
   
  } // end loop placement

  //----------------------------------------------     
  // Quartz Coordinates
  //----------------------------------------------
  G4double quartzThickness   = config->GetValue("quartzThickness", 2);
  G4double quartzSizeXPrime  = quartzThickness / TMath::Cos(theta);
  
  G4double dZprimeQuartz     =  0.5 * gapThicknessZ * TMath::Cos(theta);
  G4double zCoordPrimeQuartz =  0.5 * chamberSizeZprime - dZprimeQuartz;  
  G4double zStepPrimeQuartz  =  (totalThicknessZ + gapThicknessZ) * TMath::Cos(theta);
  
  std::vector< std::pair<G4double, G4double> > quartzCoords;
 
  while( zCoordPrimeQuartz > -0.5 * chamberSizeZprime + dZprimeQuartz){
    G4double xCoord = -0.5*chamberHalfSizeXprime + 0.5 * quartzSizeXPrime +
      zCoordPrimeQuartz * TMath::Tan( theta );
    quartzCoords.emplace_back( xCoord, zCoordPrimeQuartz);
    zCoordPrimeQuartz -= zStepPrimeQuartz; 
  }

  //----------------------------------------------     
  // Quartz Placement
  //----------------------------------------------
  G4double quartzSizeX     = quartzSizeXPrime;
  G4double quartzSizeY     = chamberSizeY;
  G4double quartzSizeZ     = gapThicknessZ * TMath::Cos(theta);
    
  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );

  m_solidQuartz = new G4Para("Quartz",             //its name
			     0.5 * quartzSizeX,    //its size
			     0.5 * quartzSizeY,
			     0.5 * quartzSizeZ,
			     0, theta, 0 );       
    
  m_logicQuartz =
    new G4LogicalVolume(m_solidQuartz,     //its solid
			m_matReflector,    //its material
			"Quartz" );        //its name

  m_logicQuartz->SetVisAttributes( quartzColor );
  
  for( unsigned int cn = 0; cn < quartzCoords.size(); cn++ ){
    m_v_physQuartz.
      push_back( new G4PVPlacement(0,
				   G4ThreeVector( quartzCoords.at(cn).first,
						  0,
						  quartzCoords.at(cn).second ),
				   m_logicQuartz,       //its logical volume
				   "Quartz",            //its name
				   m_logicChamber,      //its mother  volume
				   false,               //no boolean operation
				   cn,                  //copy number
				   checkOverlaps) );    //overlaps checking
  } // end loop placement

  //----------------------------------------------     
  // Emitter Coordinates
  //----------------------------------------------
  G4double emitterThicknessPrime = chamberHalfSizeXprime - quartzSizeXPrime;

  G4double zCoordPrimeEmitter =  0.5 * chamberSizeZprime - dZprimeEmitter;  
  G4double zStepPrimeEmitter  =  (totalThicknessZ + gapThicknessZ) * TMath::Cos(theta);
  
  std::vector< std::pair<G4double, G4double> > emitterCoords;
 
  while( zCoordPrimeEmitter > -0.5 * chamberSizeZprime + dZprimeEmitter){
    G4double xCoord = 0.5*chamberHalfSizeXprime - 0.5 * emitterThicknessPrime +
      zCoordPrimeEmitter * TMath::Tan( theta );
    emitterCoords.emplace_back( xCoord, zCoordPrimeEmitter);
    zCoordPrimeEmitter -= zStepPrimeEmitter; 
  }

  //----------------------------------------------     
  // Emitter Placement
  //----------------------------------------------
  G4double emitterSizeX     = emitterThicknessPrime;
  G4double emitterSizeY     = chamberSizeY;
  G4double emitterSizeZ     = gapThicknessZ * TMath::Cos(theta);
    
  G4VisAttributes* emitterColor  = new G4VisAttributes( G4Colour::Magenta() );

  m_solidEmitter =
    new G4Para("Emitter",             //its name
	       0.5 * emitterSizeX,    //its size
	       0.5 * emitterSizeY,
	       0.5 * emitterSizeZ,
	       0, theta, 0 );       

  m_logicEmitter =
    new G4LogicalVolume(m_solidEmitter,  //its solid
			m_matReflector,  //its material
			"Emitter" );     //its name

  m_logicEmitter->SetVisAttributes( emitterColor );
    
  
  for( unsigned int cn = 0; cn < emitterCoords.size(); cn++ ){
    m_v_physEmitter.
      push_back( new G4PVPlacement(0,
				   G4ThreeVector( emitterCoords.at(cn).first,
						  0,
						  emitterCoords.at(cn).second ),
				   m_logicEmitter,    //its logical volume
				   "Emitter",         //its name
				   m_logicChamber,    //its mother  volume
				   false,             //no boolean operation
				   cn,                //copy number
				   checkOverlaps) );  //overlaps checking
  } // end loop placement
  
  //----------------------------------------------     
  // Define Surface/Border Properties
  //----------------------------------------------  
  DefineBorderProperties();

  //----------------------------------------------     
  // SD and Scoring Volumes
  //----------------------------------------------
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  char quartzSDname[256];
  sprintf( quartzSDname, "QuartzSD%d", m_copyNumber);
  QuartzSD* aQuartzSD = new QuartzSD( quartzSDname );
  SDman->AddNewDetector( aQuartzSD );
  m_logicQuartz->SetSensitiveDetector( aQuartzSD );
}
