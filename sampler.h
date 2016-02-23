#pragma once
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void writeToFile(double localEnergy);
    void closeFile();
    void computeAverages();
    double getEnergy()          { return m_energy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_standardDeviation = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    std::fstream m_outFile;
    class System* m_system = nullptr;
};
