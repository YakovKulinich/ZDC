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
// $Id: DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class G4Box;
class G4Para;
class G4CSGSolid;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class SharedData;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  DetectorConstruction( SharedData* );
  virtual ~DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  virtual void               InitializeParameters();
  virtual void               DefineMaterials();
  
  virtual void               DefineBorderProperties();
  virtual G4VPhysicalVolume* ConstructDetector();
  
  virtual void               BuildDiagonalGeo();
  virtual void               BuildChevronGeo();
  
  G4LogicalVolume*   GetScoringVolume() const { return m_scoringVolume; }
  
protected:
  bool               m_checkOverlaps;
  
  double             m_moduleSizeX;
  double             m_moduleSizeY;
  double             m_moduleSizeZ;
  
  double             m_chamberSizeX;
  double             m_chamberSizeY;
  double             m_chamberSizeZ;

  G4Material*        m_matBox;
  G4Material*        m_matQuartz;
  G4Material*        m_matOil;
  G4Material*        m_matReflector;
  G4Material*        m_matAbsorber;
  G4Material*        m_matTop;

  G4Box*             m_solidWorld;
  G4LogicalVolume*   m_logicWorld;
  G4VPhysicalVolume* m_physWorld;

  G4Box*             m_solidBox;
  G4LogicalVolume*   m_logicBox;
  G4VPhysicalVolume* m_physBox;
  
  G4Box*             m_solidOil;
  G4LogicalVolume*   m_logicOil;
  G4VPhysicalVolume* m_physOil;

  G4Box*             m_solidTop;
  G4LogicalVolume*   m_logicTop;
  G4VPhysicalVolume* m_physTop;

  std::vector<G4Box*>              m_v_solidQuartz;
  std::vector<G4LogicalVolume*>    m_v_logicQuartz;
  std::vector<G4VPhysicalVolume*>  m_v_physQuartz;
 
  std::vector<G4CSGSolid*>         m_v_solidPanel;
  std::vector<G4LogicalVolume*>    m_v_logicPanel;
  std::vector<G4VPhysicalVolume*>  m_v_physPanel;
  
  std::vector<G4CSGSolid*>         m_v_solidAbsorber;
  std::vector<G4LogicalVolume*>    m_v_logicAbsorber;
  std::vector<G4VPhysicalVolume*>  m_v_physAbsorber;

  std::vector<G4Box*>              m_v_solidTopQuartz;
  std::vector<G4LogicalVolume*>    m_v_logicTopQuartz;
  std::vector<G4VPhysicalVolume*>  m_v_physTopQuartz;  
 
protected:
  G4LogicalVolume*  m_scoringVolume;
  SharedData*       m_sd;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

