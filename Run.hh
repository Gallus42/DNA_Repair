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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// O.V. Belov, E.A. Krasavin, M.S. Lyashko, M. Batmunkh, N.H. Sweilam. 
// A quantitative model of the major pathways for radiation-induced DNA 
// double-strand break repair. J. Theor. Biol. 2015. 366:115-130.
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// -------------------------------------------------------------------
// January 2017
// -------------------------------------------------------------------
//
/// \file Run.hh
/// \brief Definition of the Run class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"

class G4ParticleDefinition;
class DNARepairModel;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(const DetectorConstruction* detector, double a, double b, double Nirrep, double Dz);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);  
        
    virtual void Merge(const G4Run*);
    void EndOfRun();
    
  private:
    const DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin;
    G4double  m_a = 0.0;
    G4double  m_b = 0.0;
    G4double  m_Nirrep= 0.0;
    G4double  m_Dz= 0.0;

    DNARepairModel * fDNARepairModel;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

