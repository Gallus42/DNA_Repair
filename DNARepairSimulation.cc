#include "DNARepairSimulation.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "LogicHandler.hh"
#include "G4Timer.hh"
#include "G4UIExecutive.hh"
#include "G4RunManager.hh"


DNARepairSimulation::DNARepairSimulation() 
{
}

DNARepairSimulation::~DNARepairSimulation()
{
}

void DNARepairSimulation::runSimulation(const std::string& macroFile, double a, double b, double Nirrep, double Dz)
{
    // DNARepair::runSimulation(macroFile);
    // return;

    G4Timer timer;
    timer.Start();

    G4UIExecutive* ui = new G4UIExecutive(LogicHandler::argc_main, LogicHandler::argv_main);
    auto runManager = new G4RunManager();

    DetectorConstruction* detector = new DetectorConstruction;
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new ActionInitialization(a, b, Nirrep, Dz));

    runManager->Initialize();

    auto visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    ui->SessionStart();
    if (!macroFile.empty()) {
        UImanager->ApplyCommand("/control/execute " + macroFile);
    }

    timer.Stop();
    G4cout << "Simulation time = " << timer.GetRealElapsed() << " s" << G4endl;

    delete visManager;
    delete runManager;
    delete ui;
}
  
