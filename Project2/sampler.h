#pragma once
#include <fstream>
#include <vector>
#include <stdio.h>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void writeToFile(double localEnergy);
    void writeToFile();
    void openEnergyFile();
    void openPositionFile();
    void closeEnergyFile();
    void closePositionFile();
    void computeAverages();
    void clean();
    double getEnergy()          { return m_energy; }
    double getKineticEnergy()   { return m_kineticEnergy; }
    double getPotentialEnergy() { return m_potentialEnergy; }
    std::vector<double> getWaveFunctionDerivative() { return m_waveFunctionDerivative; }
    std::vector<double> getWaveFunctionEnergy() { return m_waveFunctionEnergy; }
    double getAlphaDerivative() { return m_alphaDerivative; }
    double getBetaDerivative() { return m_betaDerivative; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfParameters = 0;
    int     m_stepNumber = 0;
    double  m_localEnergy = 0;
    double  m_energy = 0;
    double  m_energySquared = 0;
    double  m_kineticEnergy = 0;
    double  m_potentialEnergy = 0;
    double  m_variance = 0;
    double  m_acceptanceRate = 0;
    std::vector<double>  m_parametersGradient = std::vector<double>();
    std::vector<double>  m_waveFunctionDerivative = std::vector<double>();
    std::vector<double>  m_waveFunctionEnergy = std::vector<double>();

    double m_alphaDerivative = 0;
    double m_betaDerivative = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double m_cumulativeKineticEnergy = 0;
    double m_cumulativePotentialEnergy = 0;
    std::vector<double>  m_cumulativeWaveFunctionDerivative = std::vector<double>();
    std::vector<double>  m_cumulativeWaveFunctionEnergy = std::vector<double>();

    std::fstream m_outFile;
    std::fstream m_outFile2;

    char m_energyFileName[50];
    char m_positionFileName[50];
    std::ofstream m_energyFileBinary;
    std::ofstream m_positionFileBinary;

    bool m_firstStep = true;
    class System* m_system = nullptr;
};
