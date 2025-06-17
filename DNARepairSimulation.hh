#ifndef DNAREPAIR_SIMULATION_HH
#define DNAREPAIR_SIMULATION_HH

#include <string>

class DNARepairSimulation
{
public:
    DNARepairSimulation();
    ~DNARepairSimulation();

    void runSimulation(const std::string& macroFile, double a, double b, double Nirrep, double Dz);

private:
};

#endif
