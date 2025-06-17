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
/// \file DNARepairModel.cc
/// \brief Implementation of the DNARepairModel class

#include "DNARepairModel.hh"
#include <iostream>
#include <fstream>
// BOOST C++ library for solving ODE's system
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <limits>
#include <cmath>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>
//define if the program is running with Geant4
#define GEANT4
#include "globals.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

using namespace std;
using namespace boost::numeric::odeint;
namespace pl = std::placeholders;
//typedef std::vector< double > state_type;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNARepairModel::DNARepairModel(G4String ParticleName, G4double KineticEnergy, 
                               G4double Dz, G4double alpha, G4double Nirrep)
{

  // Recalling model parameters
 
  fDz = Dz;
  falpha = alpha;
  fNirrep = Nirrep;

  G4cout << "\n Calculation of DNA repair procress started ... (~25 s)\n" 
         << G4endl; 

  // INITIAL CONDITIONS

  G4int NbEquat = 29;   // Total number of model equations
  state_type Y(NbEquat);
    
  //---- Initial conditions for NHEJ -----
  Y[0] = falpha; 
  Y[1] = Y[2] = Y[3] = Y[4] = Y[5] = Y[6] = Y[7] = Y[8] = Y[9] = 0.; 
   
  //---- Initial conditions for HR -------
  Y[10] = Y[11] = Y[12] = Y[13] = Y[14] = Y[15] = Y[16] = Y[17] 
        = Y[18] = Y[19] = 0.; 
  
  //---- Initial conditions for SSA -----
  Y[20] = Y[21] = Y[22] = Y[23] = Y[24] = 0.; 
  
  //---- Initial conditions for Alt-NHEJ (MMEJ) -----
  Y[25] = Y[26] = Y[27] = Y[28] = Y[5] = Y[6] = Y[7] = Y[8] = Y[9] = 0.; 

  // Integration parameters
  G4double t0 = 0.0; // Starting time point (dimensionless)
  G4double t1 = 45.3; // Final time point (dimensionless)
  G4double dt = 0.001; // Intergration time step (dimensionless)

  G4double K8 = 0.552; // [h-1], scaling variable

  vector<state_type> Y_vec;
  vector<double> times;
  size_t steps;
  steps = integrate( bind(&DNARepairModel::DSBRepairPathways, ref(*this), 
                     pl::_1, pl::_2, pl::_3), Y , t0 , t1 , dt  , 
                     push_back_state_and_time( Y_vec , times ) );

  // Output options for different repair stages
  
  G4String FociName;
  //outputPath = G4String("../Outputs/"); 
  G4String OutputFileName; 

 // "Ku recruitment"
  // if (CommandLineParser::GetParser()->GetCommandIfActive("-Ku"))
   {
    FociName = G4String("Ku"); 
    G4cout << " \n Ku is activated!"<<G4endl;
 
    OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";

    remove (OutputFileName);

    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][1] << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl; 
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of "Ku recruitment"

  // "DNAPKcs recruitment"
  // else if (CommandLineParser::GetParser()->GetCommandIfActive("-DNAPKcs"))
   {
    FociName = G4String("DNAPKcs"); 
    G4cout << " \n DNAPKcs is activated!"<<G4endl;

    if (ParticleName == "gamma")
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                           + "MeV/n" + ParticleName + "Simulation.out";
     } 

    else if (ParticleName == "Fe56" && KineticEnergy/56==1000)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/56).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else 
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";
     }

    remove (OutputFileName);

    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][3]+Y_vec[i][4]+Y_vec[i][5]+Y_vec[i][6]+Y_vec[i][7]
        << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl; 
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of "DNAPKcs recruitment"

  // "Rad51 recruitment"
  // else if (CommandLineParser::GetParser()->GetCommandIfActive("-Rad51"))
   {
    FociName = G4String("Rad51"); 
    G4cout << " \n Rad51 is activated!"<<G4endl;

    if (ParticleName == "gamma")
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                           + "MeV" + ParticleName + "Simulation.out";
     } 

    else if (ParticleName == "Fe56" && KineticEnergy/56==1000)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/56).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else 
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";
     }

    remove (OutputFileName);

    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][15]+Y_vec[i][16]+Y_vec[i][17]
        << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl; 
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of "Rad51 recruitment"

   // "RPA recruitment"
  // else if (CommandLineParser::GetParser()->GetCommandIfActive("-RPA"))
   {
    FociName = G4String("RPA"); 
    G4cout << " \n RPA is activated!"<<G4endl;

    if (ParticleName == "gamma")
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                           + "MeV" + ParticleName + "Simulation.out";
     } 

    else 
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";
     }

    remove (OutputFileName);

    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][14]+Y_vec[i][15]+Y_vec[i][20]
        << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl; 
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of "RPA recruitment"

  // "gH2AX foci" 
  // else if (CommandLineParser::GetParser()->GetCommandIfActive("-gH2AX"))
   {
    FociName = G4String("gH2AX"); 
    G4cout << " \n gH2AX is activated!"<<G4endl;

    if (ParticleName == "gamma")
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                           + "MeV" + ParticleName + "Simulation.out";
     } 

    else if (ParticleName == "O16" && KineticEnergy/16==1000)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/16).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else if (ParticleName == "Si28" && KineticEnergy/28==1000)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/28).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else if (ParticleName == "Fe56" && KineticEnergy/56==300)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/56).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else if (ParticleName == "Fe56" && KineticEnergy/56==1000)
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy/56).substr(0,4)
                           + "MeVu" + ParticleName + "Simulation.out";
     }

    else 
     {
      OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";
     }

    remove (OutputFileName);
    
    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][9] << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl;  
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of "gH2AX foci"

   // User specified output
  // else if (CommandLineParser::GetParser()->
  // GetCommandIfActive("-YourEnzyme")) // Enter your enzyme name
   {
    FociName = G4String("YourEnzyme"); // Enter your enzyme name
    G4cout << " \n YourEnzyme is activated!"<<G4endl; // Enter your enzyme name

    OutputFileName = 
                  FociName +"_" + std::to_string(KineticEnergy).substr(0,4)
                             + "AMeV" + ParticleName + "Simulation.out";
    remove (OutputFileName);

    std::ofstream WriteFoci (OutputFileName, std::ios::app);
    for( size_t i=0; i<=steps; i+=16000 )
     {
       WriteFoci << times[i]*K8 << '\t' << "   "
        << Y_vec[i][1] //Enter a variable or a combination of variables you wish
        << '\t' << "   " << G4endl;
     }
    G4cout << " ... ... ... " << "\n Calculation finished! \n" <<G4endl; 
    G4cout << " Output file writen into " << OutputFileName <<"\n"<<G4endl;

   } // end of User specified output

  //G4cout << " Results have been written into 'Outputs' directory"<<G4endl;  
   {
       G4cout << " \n Exporting Y's" << G4endl; // Enter your enzyme name

       for (int y = 0; y < 29; ++y)
       {
           OutputFileName = "y/Y_" + std::to_string(y) + ".out";
           remove(OutputFileName);

           std::ofstream WriteFoci(OutputFileName, std::ios::app);
           for (size_t i = 0; i <= steps; i += 16000)
           {
               WriteFoci << times[i] * K8 << '\t' << "   "
                   << Y_vec[i][y] 
                   << '\t' << "   " << G4endl;
           }
       }

   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNARepairModel::~DNARepairModel()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNARepairModel::DSBRepairPathways( const state_type &Y , 
                                        state_type &YP , double t )
{

  // Concentrations of repair enzymes set to be constant
  G4double X1; // [Ku]
  G4double X2; // [DNAPKcsArt]
  G4double X3; // [LigIV/XRCC4/XLF]
  G4double X4; // [PNKP]
  G4double X5; // [Pol]
  G4double X6; // [H2AX]
  G4double X7; // [MRN/CtIP/ExoI/Dna2]
  G4double X8; // [ATM]
  G4double X9; // [RPA]
  G4double X10; // [Rad51/Rad51par/BRCA2]
  G4double X11; // [DNAinc]
  G4double X12; // [Rad52]
  G4double X13; // [ERCC1/XPF]
  G4double X14; // [LigIII]
  G4double X15; // [PARP1]
  G4double X16; // [Pol]
  G4double X17; // [LigI]

  X1 = X2 = X3 = X4 = X5 = X6 = X7 = X8 = X9 = X10 = X11 = X12 = X13 = X14 =
  X15 = X16 = X17 = 400000; 

  fTime = t;  // Recalling t

  //  DIMENSIONAL REACTION RATES
  //
  //------------NHEJ--------------
  G4double K1 = 11.052;             // M-1*h-1 
  G4double Kmin1 = 6.59999*1e-04;   // h-1
  G4double K2 = 18.8305*(1.08517-exp(-21.418/pow(fDz,1.822))); // M-1*h-1 
  G4double Kmin2 = 5.26*1e-01;      //h-1
  G4double K3 = 1.86;               // h-1
  G4double K4 = 1.38*1e+06;          // M-1*h-1
  G4double Kmin4 = 3.86*1e-04;      // h-1
  G4double K5 = 15.24;              // M-1*h-1
  G4double Kmin5 = 8.28;            // h-1
  G4double K6 = 18.06;              // M-1*h-1
  G4double Kmin6 = 1.33;            // h-1
  G4double K7 = 2.73*1e+05;          // M-1*h-1
  G4double Kmin7 = 3.2;             // h-1
  G4double K8 = 5.52*1e-01;         // h-1
  G4double K9 = 1.66*1e-01;         // h-1
  G4double K10 = 1.93*1e-07/fNirrep; // M
  G4double K11 = 7.50*1e-02;        // h-1
  G4double K12 = 11.1;              // h-1
  //
  //------------HR--------------
  G4double P1 = 1.75*1e+03;          // M-1*h-1
  G4double Pmin1 = 1.33*1e-04;      // h-1
  G4double P2 = 0.39192;            // h-1
  G4double Pmin2 = 2.7605512*1e+02;  // h-1
  G4double P3 = 1.37*1e+04;          // M-1*h-1
  G4double Pmin3 = 2.34;            // h-1
  G4double P4 = 3.588*1e-02;       // h-1
  G4double P5 = 1.20*1e+05;          // M-1*h-1
  G4double Pmin5 = 8.82*1e-05;      // h-1
  G4double P6 = 1.54368*1e+06;       // M-1*h-1
  G4double Pmin6 = 1.55*1e-03;      // h-1
  G4double P7 = 1.4904;              // h-1
  G4double P8 = 1.20*1e+04;          // M-1*h-1
  G4double Pmin8 = 2.49*1e-04;      // h-1
  G4double P9 = 1.104;            //h-1
  G4double P10 = 7.20*1e-03;        // h-1
  G4double P11 = 6.06*1e-04;        // h-1
  G4double P12 = 2.76*1e-01;        // h-1
  //
  //------------SSA--------------
  G4double Q1 = 1.9941*1e+05;       // M-1*h-1
  G4double Qmin1 = 1.71*1e-04;      // h-1
  G4double Q2 = 4.8052*1e+04;       // M-1*h-1
  G4double Q3 = 6*1e+03;             // M-1*h-1
  G4double Qmin3 = 6.06*1e-04;      // h-1
  G4double Q4 = 1.66*1e-06;         // h-1
  G4double Q5 = 8.40*1e+04;          // M-1*h-1
  G4double Qmin5 = 4.75*1e-04;      // h-1
  G4double Q6 = 11.58;              // h-1
  //
  //-------alt-NHEJ (MMEJ)--------
  G4double R1 = 2.39*1e+03;          // M-1*h-1
  G4double Rmin1 = 12.63;           // h-1
  G4double R2 = 4.07*1e+04;          // M-1*h-1
  G4double R3 = 9.82;               // h-1
  G4double R4 = 1.47*1e+05;          // M-1*h-1
  G4double Rmin4 = 2.72;            // h-1
  G4double R5 = 1.65*1e-01;         //h-1
  //
  // Scalling rate XX1
  G4double XX1 = 9.19*1e-07; // M
  //
  // DIMENSIONLESS REACTION RATES 
  //
  //------------NHEJ--------------
  G4double k1 = K1*XX1/K8;
  G4double kmin1 = Kmin1/K8;
  G4double k2 = K2*XX1/K8; 
  G4double kmin2 = Kmin2/K8;     
  G4double k3 = K3/K8;               
  G4double k4 = K4*XX1/K8;         
  G4double kmin4 = Kmin4/K8;      
  G4double k5 = K5*XX1/K8;             
  G4double kmin5 = Kmin5/K8;            
  G4double k6 = K6*XX1/K8;              
  G4double kmin6 = Kmin6/K8;            
  G4double k7 = K7*XX1/K8;          
  G4double kmin7 = Kmin7/K8;           
  G4double k8 = K8/K8;         
  G4double k9 = K9/K8;         
  G4double k10 = K10/XX1; 
  G4double k11 = K11/K8;       
  G4double k12 = K12/K8;            
  //
  //------------HR--------------
  G4double p1 = P1*XX1/K8;       
  G4double pmin1 = Pmin1/K8;     
  G4double p2 = P2/K8;           
  G4double pmin2 = Pmin2/K8;  
  G4double p3 = P3*XX1/K8;          
  G4double pmin3 = Pmin3/K8;        
  G4double p4 = P4/K8;       
  G4double p5 = P5*XX1/K8;       
  G4double pmin5 = Pmin5/K8;       
  G4double p6 = P6*XX1/K8;        
  G4double pmin6 = Pmin6/K8;       
  G4double p7 = P7/K8;              
  G4double p8 = P8*XX1/K8;          
  G4double pmin8 = Pmin8/K8;      
  G4double p9 = P9/K8;          
  G4double p10 = P10/K8;        
  G4double p11 = P11/K8;        
  G4double p12 = P12/K8;       
  //
  //------------SSA--------------
  G4double q1= Q1*XX1/K8;       
  G4double qmin1 = Qmin1/K8;      
  G4double q2 = Q2*XX1/K8;       
  G4double q3 = Q3*XX1/K8;        
  G4double qmin3 = Qmin3/K8;     
  G4double q4 = Q4/K8;         
  G4double q5 = Q5*XX1/K8;          
  G4double qmin5 = Qmin5/K8;     
  G4double q6 = Q6/K8;             
  //
  //-------alt-NHEJ (MMEJ)--------
  G4double r1 = R1*XX1/K8;         
  G4double rmin1 = Rmin1/K8;        
  G4double r2 = R2*XX1/K8;         
  G4double r3 = R3/K8;            
  G4double r4 = R4*XX1/K8;         
  G4double rmin4 = Rmin4/K8;           
  G4double r5 = R5/K8;        
  //------------------------------------

  // SYSTEM OF DIFFERENTIAL EQUATIONS

  // ----- NHEJ ----------
  YP[0] = fNirrep - k1*Y[0]*X1 + kmin1*Y[1] - p1*Y[0]*X1 + pmin1*Y[10]; // [DSB]

  YP[1] = k1*Y[0]*X1 - kmin1*Y[1] - k2*Y[1]*X2 + kmin2*Y[2]; // [DBS * Ku]

  YP[2] = k2*Y[1]*X2 - k3*Y[2] - kmin2*Y[2]; // [DSB * DNA-PK/Art]

  YP[3] = k3*Y[2] - k4*(Y[3]*Y[3]) + kmin4*Y[4]; // [DSB * DNA-PK/ArtP]

  YP[4] = k4*(Y[3]*Y[3]) - kmin4*Y[4] - k5*Y[4]*X3 + kmin5*Y[5]; // [Bridge]

  YP[5] = kmin6*Y[6] + k5*Y[4]*X3 - kmin5*Y[5] - k6*Y[5]*X4;
  // [Bridge * LigIV/XRCC4/XLF]

  YP[6] = -kmin6*Y[6] - k7*Y[6]*X5 +  kmin7*Y[7] + k6*Y[5]*X4; 
  // [Bridge * LigIV/XRCC4/XLF * PNKP]

  YP[7] = k7*Y[6]*X5 - k8*Y[7] - kmin7*Y[7];
  // [Bridge * LigIV/XRCC4/XLF * PNKP * Pol]

  YP[8] = r5*Y[28] + k8*Y[7] + p12*Y[18] + p11*Y[19] + q6*Y[24]; // [dsDNA]

  YP[9] = (k9*(Y[3] + Y[4] + Y[5] + Y[6] + Y[7])*X6)/(k10 + Y[3] + Y[4] + Y[5] 
          + Y[6] + Y[7]) - k11*Y[8] - k12*Y[9]; // [gH2AX foci]

  // ----- HR ---------- 
  YP[10] = p1*Y[0]*X7 - pmin1*Y[10] - p3*Y[10]*Y[11] + pmin3*Y[12]; 
  // [MRN/CtIP/ExoI/Dna2]

  YP[11] = p2*X8 - pmin2*Y[11] - p3*Y[10]*Y[11] + p4*Y[12] + pmin3*Y[12]; 
  // [ATMP]

  YP[12] = p3*Y[10]*Y[11] - p4*Y[12] - pmin3*Y[12]; 
  // [DSB * MRN/CtIP/ExoI/Dna2 * ATMP]

  YP[13] = rmin1*Y[25] + p4*Y[12] - r1*X15*Y[13] - p5*Y[13]*X9 + pmin5*Y[14]; 
  // [ssDNA]

  YP[14] = pmin6*Y[15] + p5*Y[13]*X9 - pmin5*Y[14] - p6*Y[14]*X10 - 
           q1*Y[14]*X12 + qmin1*Y[20]; // [ssDNA * RPA]

  YP[15] = -p7*Y[15] - pmin6*Y[15] + p6*Y[14]*X10;
  // [ssDNA * RPA * Rad51/Rad51par/BRCA2

  YP[16] = p7*Y[15] - p8*Y[16]*X11 + pmin8*Y[17]; // [Rad51 filament]

  YP[17] = p8*Y[16]*X11 - p9*Y[17] - pmin8*Y[17]; // [Rad51 filament * DNAinc]

  YP[18] = p9*Y[17] - p10*Y[18] - p12*Y[18]; // [D-loop]

  YP[19] = p10*Y[18] - p11*Y[19]; // [dHJ]

  // ----- SSA ---------- 
  YP[20] = q1*Y[14]*X12 - qmin1*Y[20] -  q2*(Y[20]*Y[20]);
  // [ssDNA * RPA * Rad52]

  YP[21] = q2*(Y[20]*Y[20]) - q3*Y[21]*X13 + qmin3*Y[22]; // [Flap]

  YP[22] = q3*Y[21]*X13 - q4*Y[22] - qmin3*Y[22]; // [Flap * ERCC1/XPF]

  YP[23] = q4*Y[22] - q5*Y[23]*X14 + qmin5*Y[24]; // [dsDNAnicks]

  YP[24] = q5*Y[23]*X14 - q6*Y[24] - qmin5*Y[24]; // [dsDNAnicks * LigIII]

  // ----- MMEJ ---------- 
  YP[25] = -rmin1*Y[25] - r2*Y[25]*X16 + r1*X15*Y[13]; // [ssDNA * PARP1]

  YP[26] = r2*Y[25]*X16 - r3*Y[26]; // [ssDNA * Pol]

  YP[27] = r3*Y[26] - r4*Y[27]*X17 + rmin4*Y[28]; // [MicroHomol]

  YP[28] = r4*Y[27]*X17 - r5*Y[28] - rmin4*Y[28]; // [MicroHomol * LigI]

  //---------------------------------------------------

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
