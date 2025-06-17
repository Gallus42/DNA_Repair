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
/// This example is provided by the Geant4-DNA collaboration
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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4EmDNAPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmDNAPhysicsActivator.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmParameters.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//PhysicsList::PhysicsList() 
//: G4VModularPhysicsList()
//{
//  SetDefaultCutValue(1.0*micrometer);
//  SetVerboseLevel(1);
//  
//  // FIRST METHOD TO ACTIVATE Geant4-DNA Physics, 
//  //  using a Geant4-DNA Physics constructor only
//  //
//   // RegisterPhysics(new G4EmDNAPhysics_option4());
//  
//  // or SECOND METHOD TO ACTIVATE Geant4-DNA Physics
//  // (this includes combination with Geant4 EM Physics)
//
//  RegisterPhysics(new G4EmStandardPhysics());
//  // RegisterPhysics(new G4EmDNAPhysicsActivator());
//
//  //G4EmParameters::Instance()->AddDNA("Target","Opt0");
//
//  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
//  G4EmParameters* emParams = G4EmParameters::Instance();
//  std::ofstream out("G4EmParameters.txt");
//  G4cout.rdbuf(out.rdbuf());
//  emParams->Dump();
//   //G4cout << "Number of particles: " << G4ParticleTable::GetParticleTable()->GetTotalNum() << G4endl;
//}
PhysicsList::PhysicsList()
    : G4VModularPhysicsList()
{
    RegisterPhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4EmDNAPhysicsActivator());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructParticle()
{
    G4Gamma::GammaDefinition();
    G4VModularPhysicsList::ConstructParticle();
}

