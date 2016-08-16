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

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
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
  : G4VUserDetectorConstruction(),
    m_scoringVolume(0),
    m_sd(NULL)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( SharedData* sd )
  : G4VUserDetectorConstruction(),
    m_scoringVolume(0),
    m_sd( sd )
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // Box dimensions
  G4double boxSizeX     = config->GetValue( "boxSizeX", 90);
  G4double boxSizeY     = config->GetValue( "boxSizeY", 240);
  G4double boxSizeZ     = config->GetValue( "boxSizeZ", 150);
  
  //----------------------------------------------     
  // World
  //----------------------------------------------
  G4double world_sizeX  = 1.2 * boxSizeX;   // mm
  G4double world_sizeY  = 1.2 * boxSizeY;   // mm
  G4double world_sizeZ  = 1.2 * boxSizeZ;   // mm
  
  printf( "Building world with %5.1f x %5.1f x %5.1f\n",
	  world_sizeX, world_sizeY, world_sizeZ );
  
  m_solidWorld =    
    new G4Box("World",                       //its name
	      0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);   //its size

  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");
  
  m_logicWorld =                         
    new G4LogicalVolume(m_solidWorld,        //its solid
                        g4Air,           //its material
                        "World");            //its name
                                   
  m_physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      m_logicWorld,          //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //----------------------------------------------     
  // Box
  //---------------------------------------------- 
  m_solidBox =    
    new G4Box("Box",                         //its name
	      0.5*boxSizeX, 0.5*boxSizeY, 0.5*boxSizeZ);  //its size

  G4Material* g4Al = nist->FindOrBuildMaterial("G4_Al");
  
  m_logicBox =                         
    new G4LogicalVolume(m_solidBox,          //its solid
                        g4Air,           //its material
                        "Box");              //its name
                                   
  m_physBox = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      m_logicBox,            //its logical volume
                      "Box",                 //its name
                      m_logicWorld,          //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //----------------------------------------------     
  // Quartz
  //----------------------------------------------
  G4double quartzThickness = config->GetValue( "quartzThickness", 2);
  G4double boxThickness = config->GetValue( "boxThickness", 2);

  std::vector< std::pair<G4double, G4double> > quartzCenters;
  quartzCenters.emplace_back( 0,  0.5*boxSizeY - boxThickness - 0.5*quartzThickness ); // top
  quartzCenters.emplace_back( 0, -0.5*boxSizeY + boxThickness + 0.5*quartzThickness ); // bottom
  quartzCenters.emplace_back( -0.5*boxSizeX + boxThickness + 0.5*quartzThickness, 0 ); // left
  quartzCenters.emplace_back(  0.5*boxSizeX - boxThickness - 0.5*quartzThickness, 0 ); // right
  
  std::vector< std::pair<G4double, G4double> > quartzSizeXY;
  quartzSizeXY.emplace_back( boxSizeX - 2*boxThickness - 2*quartzThickness,
			     quartzThickness );          // top
  quartzSizeXY.emplace_back( boxSizeX - 2*boxThickness - 2*quartzThickness,
			     quartzThickness );          // bottom
  quartzSizeXY.emplace_back( quartzThickness,
			     boxSizeY - 2*boxThickness); // left
  quartzSizeXY.emplace_back( quartzThickness,
			     boxSizeY - 2*boxThickness); // right
  
  G4double quartzSizeZ = boxSizeZ - 2*boxThickness;

  G4Material* g4SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4VisAttributes* quartzColor = new G4VisAttributes( G4Colour::Blue() );
  
  for( unsigned int qn = 0; qn < quartzCenters.size(); qn++ ){
    m_v_solidQuartz.
      push_back( new G4Box("Quartz",            //its name
			   0.5*quartzSizeXY.at(qn).first,
			   0.5*quartzSizeXY.at(qn).second,
			   0.5*quartzSizeZ) );  //its size
    
    m_v_logicQuartz.
      push_back( new G4LogicalVolume(m_v_solidQuartz.back(), //its solid
				     g4SiO2,                 //its material
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
				   checkOverlaps) );        //overlaps checking
  }

  //----------------------------------------------     
  // Oil
  //----------------------------------------------
  G4double oilSizeX     = boxSizeX - 2*boxThickness - 2*quartzThickness;
  G4double oilSizeY     = boxSizeY - 2*boxThickness - 2*quartzThickness;
  G4double oilSizeZ     = boxSizeZ - 2*boxThickness - 2*quartzThickness;

  std::cout << oilSizeX << " --- " << oilSizeY << " --- " << oilSizeZ << std::endl;
  
  m_solidOil =    
    new G4Box("Oil",                         //its name
	      0.5*oilSizeX, 0.5*oilSizeY, 0.5*oilSizeZ);  // its size
  
  m_logicOil =                         
    new G4LogicalVolume(m_solidOil,          //its solid
                        g4Air,               //its material
                        "Oil");              //its name
  
  m_physOil = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      m_logicOil,            //its logical volume
                      "Oil",                 //its name
                      m_logicBox,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  
  G4VisAttributes* oilColor = new G4VisAttributes( G4Colour::Magenta() );
  m_logicOil->SetVisAttributes( oilColor );

  
  /*
  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  double xPts[4]; double yPts[4];

  double absorberThicknessZ = config->GetValue( "absorberThicknessZ", 5 );
  double      gapThicknessZ = config->GetValue( "gapThicknessZ", 10 );

  double thickness = absorberThicknessZ * cos(theta);
  double       gap =      gapThicknessZ * cos(theta);
  */
  
  // Set Box as scoring volume
  //
  m_scoringVolume = m_logicBox;

  //
  //always return the physical World
  //
  return m_physWorld;
}

double DetectorConstruction :: ComputePlacementParameters( double center[2] ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
