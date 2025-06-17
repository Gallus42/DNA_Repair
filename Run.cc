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
/// \file Run.cc
/// \brief Implementation of the Run class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Exp.hh"
#include "G4EmCalculator.hh"
#include "G4UImanager.hh"

// DNA double-strand break repair model
#include "DNARepairModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* detector, double a, double b, double Nirrep, double Dz)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),
    m_a(a), m_b(b), m_Nirrep(Nirrep), m_Dz(Dz)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary (G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin     = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(2);
  
  //Run conditions  
  //
  G4Material* material = fDetector->GetTargetMat();
  //G4double density  = material->GetDensity();   
  // Characteristics of Primary particle 
  G4String Particle = fParticle->GetParticleName();

  //Stopping Power from input Table
  G4EmCalculator emCalculator;
  G4double dEdxFull = 0.;
  if (fParticle->GetPDGCharge()!= 0.) 
  { 
    dEdxFull  = emCalculator.ComputeTotalDEDX(fEkin,fParticle,material)/(keV/um);    
  }
 
  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
  G4double Ekin = fEkin/MeV; 

 // -------------------------------------------------------------------------//
   G4cout << " \n The function of DSB yield on LET is activated!"<<G4endl;

  // Parameters of radaition
  G4double Lz;  // Linear energy transfer (keV/um)
  G4double Dz; // Radiation dose (Gy)

  // Setting the function of double-strand break yield  
  // on linear energy transfer (dimensionless)
  G4double alpha;

  // Parameters of repair efficiency
  G4double Nirrep = m_Nirrep; // Fraction of irrepairable DSBs

  if (0.0 == Nirrep)
  {
      if (Particle == "gamma") // for gamma-rays and X-rays
      {
          Nirrep = 0.01;
          Lz = 0.2;
      }
      else if (Particle == "O16" && Ekin/16==1000)
      {
          Nirrep = 0.01;
      }
      else if (Particle == "Si28" && Ekin/28==1000)
      {
          Nirrep = 0.02;
      }
      else if (Particle == "Fe56" && Ekin/56==1000)
      {
          Nirrep = 0.2;
      }
      else if (Particle == "Fe56" && Ekin/56==300)
      {
          Nirrep = 0.3;
      }
      else
      {
          Nirrep = 0.2; // For other radiation modalities or in the case of cpecific
          //cells lines (e.g mutants) Nirrep needs to be set manually. It is a free
          //parameter of the model.
      }
  }

  Lz = dEdxFull; // Automatic calculation charged particle LET in liquid water

  Dz = 0.0 == m_Dz ? 1. : Dz;               // Set the radaition dose
  G4double a = 0 == m_a ? 27.5 : m_a;       // "a" coefficient
  G4double b = 0 == m_b ? 0.00243 : m_b;    // "b" coefficient

  alpha = a*G4Exp(-b*Lz); // Function of DSB yield on particles LET

  // -----------
  // TODO: задавать параметры Dz и Nirrep через UI
  // -----------


  G4cout << "\n ===== The input parameters for 'DNARepair' simulation ======\n";
  G4cout << "\n  Primary particle                     = " << Particle
         << "\n  Kinetic energy (A*MeV)               = " << Ekin
         << "\n  Full LET (keV/um)                    = " << Lz 
         << "\n  Radiation dose (Gy)                  = " << Dz   
         << "\n  Total yield of DNA DSBs (Gy-1 /cell) = " << alpha  
         << "\n  Fraction of irreparable DSBs (%)     = " << Nirrep
         << G4endl;   

  // To call DNA DSBs repair model:
  fDNARepairModel = new DNARepairModel(Particle, Ekin, Dz, alpha, Nirrep);
 
 // -------------------------------------------------------------------------//

  // Reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
