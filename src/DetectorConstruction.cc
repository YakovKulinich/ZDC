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
#include "EMCal.hh"

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
  : G4VUserDetectorConstruction(), m_sd(NULL),
    m_solidWorld(NULL), m_logicWorld(NULL), m_physWorld(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( SharedData* sd )
  : G4VUserDetectorConstruction(), m_sd( sd ), 
    m_solidWorld(NULL), m_logicWorld(NULL), m_physWorld(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
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
  
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction :: DefineBorderProperties()
{}

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
    new G4LogicalVolume(m_solidWorld,     //its solid
                        g4Air,            //its material
                        "World");         //its name
                                   
  m_physWorld = 
    new G4PVPlacement(0,                  //no rotation
                      G4ThreeVector(),    //at (0,0,0)
                      m_logicWorld,       //its logical volume
                      "World",            //its name
                      0,                  //its mother  volume
                      false,              //no boolean operation
                      0,                  //copy number
                      checkOverlaps);     //overlaps checking

  //----------------------------------------------     
  // Build EMCal Module
  //----------------------------------------------
  m_emCal = new EMCal("EMCAL1", NULL, G4ThreeVector(), m_logicWorld, m_sd );
  m_emCal->Construct();
  
  //----------------------------------------------     
  // SD and Scoring Volumes
  //----------------------------------------------  
  // Set Housing as scoring volume
  
  return m_physWorld;
}
