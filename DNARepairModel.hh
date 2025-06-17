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
// J. Rad. Res. App. Sci. 8 (2015) 498-507
// O.V. Belov, M. Batmunkh, S. Incerti, O. Lkhagva, (2016) 
//   Radiation damage to neuronal cells: Simulating the energy deposition 
//  and water radiolysis in a small neural network. 
//  Physica Medica (accepted on 1 November 2016). 
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $Id: 
// 
/// \file DNARepairModel.hh
/// \brief Implementation of the DNARepairModel class

#ifndef DNARepairModel_H
#define DNARepairModel_H 1
//define if the program is running with Geant4
#define GEANT4

#ifdef GEANT4
//Specific to Geant4, globals.hh is used for G4cout
#include "globals.hh"
#else
#define G4Exp  std::exp
#define G4Pow  std::pow  
#define G4cout std::cout
#define G4cerr std::cerr 
#define G4endl std::endl 
#define G4String std::string 
#include <cfloat>
#endif
#include <iostream>
#include <fstream>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNARepairModel 
{ 
  public:  
    DNARepairModel(G4String ParticleName, G4double KineticEnergy, 
                   G4double Dz, G4double alpha, G4double Nirrep);
  virtual     
   ~DNARepairModel();  

  //typedef boost::array< double , 29 > state_type;
  typedef std::vector< double > state_type;

  void DSBRepairPathways( const state_type &Y , state_type &YP , double t );

  struct push_back_state_and_time
  {
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , 
                              std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &Y , double t )
    {
        m_states.push_back( Y );
        m_times.push_back( t );
    }
  };
     
  private:

  G4double fDz;
  G4double falpha;
  G4double fNirrep; 
  G4double fTime;
 
};

#endif

