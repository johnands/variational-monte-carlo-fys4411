#pragma once
#include <fstream>
#include <vector>

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
    std::vector<double> getWaveFunctionDerivative() { return m_waveFunctionDerivative; }
    std::vector<double> getWaveFunctionEnergy() { return m_waveFunctionEnergy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfParameters = 0;
    int     m_stepNumber = 0;
    double  m_localEnergy = 0;
    double  m_energy = 0;
    double  m_standardDeviation = 0;
    std::vector<double>  m_waveFunctionDerivative = std::vector<double>();
    std::vector<double>  m_waveFunctionEnergy = std::vector<double>();
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    std::vector<double>  m_cumulativeWaveFunctionDerivative = std::vector<double>();
    std::vector<double>  m_cumulativeWaveFunctionEnergy = std::vector<double>();
    std::fstream m_outFile;
    std::fstream m_outFile2;
    bool m_firstStep = true;
    class System* m_system = nullptr;
};
