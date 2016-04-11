#pragma once
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep, bool writeEnergiesToFile, bool writePositionsToFile);
    void printOutputToTerminal();
    void writeToFile(double localEnergy);
    void writeToFile();
    void closeFile();
    void computeAverages();
    void clean();
    double getEnergy()          { return m_energy; }
    double getWaveFunctionDerivative() { return m_waveFunctionDerivative; }
    double getWaveFunctionEnergy() { return m_waveFunctionEnergy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_standardDeviation = 0;
    double  m_waveFunctionDerivative = 0;
    double  m_waveFunctionEnergy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cumulativeWaveFunctionDerivative = 0;
    double  m_cumulativeWaveFunctionEnergy = 0;
    std::fstream m_outFile;
    std::fstream m_outFile2;
    class System* m_system = nullptr;
};
