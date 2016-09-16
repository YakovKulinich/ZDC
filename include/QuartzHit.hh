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
// January 22, 2012
// Yakov Kulinich
// Stony Brook & BNL
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef QuartzHit_h
#define QuartzHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class QuartzHit : public G4VHit
{
public:

  QuartzHit();
  ~QuartzHit();
  QuartzHit(const QuartzHit&);
  const QuartzHit& operator=(const QuartzHit&);
  G4int operator==(const QuartzHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:
  void setTrackID    (G4int track)        { trackID = track; };
  void setRodNb      (G4int rod)          { rodNb = rod; };  
  void setEdep       (G4double de)        { edep = de; };
  void setPos        (G4ThreeVector xyz)  { pos = xyz; };
      
  G4int         getTrackID()    { return trackID; };
  G4int         getRodNb()      { return rodNb; };
  G4double      getEdep()       { return edep; };      
  G4ThreeVector getPos()        { return pos; };
      
private:
  G4int         trackID;
  G4int         rodNb;
  G4double      edep;
  G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<QuartzHit> QuartzHitsCollection;

extern G4Allocator<QuartzHit> QuartzHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* QuartzHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) QuartzHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void QuartzHit::operator delete(void *aHit)
{
  QuartzHitAllocator.FreeSingle((QuartzHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
