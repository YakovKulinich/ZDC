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
//
// $Id: QuartzSD.cc,v 1.9 2006/06/29 17:48:27 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "QuartzSD.hh"
#include "SharedData.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TString.h"

#include <string>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QuartzSD::QuartzSD(G4String sdName, SharedData* sd)
  :G4VSensitiveDetector(sdName), m_sd(sd){
  collectionName.insert(sdName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QuartzSD::~QuartzSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::HistInitialize(){
  // G4String is strange...
  std::string name = GetName();
  
  // Add some histograms
  h2_rodNum_eDep = new TH2D( Form("h2_rodNum_eDep_%s", name.c_str() ),
				  ";rod number;eDep [keV]",
				  14,0,14,
				  50,0,25);
  m_sd->AddOutputHistogram( h2_rodNum_eDep );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::Initialize(G4HCofThisEvent* HCE){
  quartzCollection = new QuartzHitsCollection(SensitiveDetectorName,
					      collectionName[0]); 

  std::string name = collectionName[0];					    
  
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name ); }
  HCE->AddHitsCollection( HCID, quartzCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool QuartzSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){
  
  G4double eDep   = aStep->GetTotalEnergyDeposit();
  G4int    modNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);
  G4int    rodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);

  h2_rodNum_eDep->Fill( rodNum, eDep );
  
  QuartzHit* newHit = new QuartzHit();

  newHit->setTrackID   ( aStep->GetTrack()->GetTrackID() );
  newHit->setModNb     ( modNum );
  newHit->setRodNb     ( rodNum );
  newHit->setEdep      ( eDep );
  newHit->setPos       ( aStep->GetPostStepPoint()->GetPosition() );

  quartzCollection->insert( newHit );
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NbHits = quartzCollection->entries();
  if(verboseLevel>0) { 
    G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
	   << " hits in the calorimeter cells: " << G4endl;
    for (G4int i=0;i<NbHits;i++) (*quartzCollection)[i]->Print();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

