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
// $Id: EMCal.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file EMCal.hh
/// \brief Definition of the EMCal class

#ifndef EMCal_h
#define EMCal_h 1

#include "globals.hh"
#include "G4PVPlacement.hh"

#include <vector>

class G4Box;
class G4Para;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class SharedData;

/// Detector construction class to define materials and geometry.

class EMCal
{
public:
  EMCal( const std::string&, const G4RotationMatrix*,
	 const G4ThreeVector&, G4LogicalVolume*, SharedData* );
  EMCal();
  ~EMCal();
  
  virtual void  Construct();
  
  virtual void  DefineMaterials();  
  virtual void  DefineBorderProperties();

  virtual void  ConstructDetector();
  
protected:
  const std::string    m_name;
  
  const G4RotationMatrix*    m_rot;
  const G4ThreeVector        m_pos;

  G4LogicalVolume*   m_logicMother;

  SharedData*        m_sd;

protected:
  G4Material*        m_matHousing;
  G4Material*        m_matQuartz;
  G4Material*        m_matEmitter;
  G4Material*        m_matReflector;
  G4Material*        m_matAbsorber;

  G4Para*            m_solidHousing;
  G4LogicalVolume*   m_logicHousing;
  G4VPhysicalVolume* m_physHousingL;
  G4VPhysicalVolume* m_physHousingR;
  
  G4Para*            m_solidChamber;
  G4LogicalVolume*   m_logicChamber;
  G4VPhysicalVolume* m_physChamber;
 
  G4Para*                          m_solidPanel;
  G4LogicalVolume*                 m_logicPanel;
  std::vector<G4VPhysicalVolume*>  m_v_physPanel;
  
  G4Para*                          m_solidAbsorber;
  G4LogicalVolume*                 m_logicAbsorber;
  std::vector<G4VPhysicalVolume*>  m_v_physAbsorber;

  G4Para*                          m_solidQuartz;
  G4LogicalVolume*                 m_logicQuartz;
  std::vector<G4VPhysicalVolume*>  m_v_physQuartz;

  G4Para*                          m_solidEmitter;
  G4LogicalVolume*                 m_logicEmitter;
  std::vector<G4VPhysicalVolume*>  m_v_physEmitter;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

