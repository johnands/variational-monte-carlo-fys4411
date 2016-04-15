#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <iomanip>

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep, bool writeEnergiesToFile, bool writePositionsToFile) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_energy = 0;
        m_waveFunctionDerivative = 0;
        m_waveFunctionEnergy = 0;
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeWaveFunctionDerivative = 0;
        m_cumulativeWaveFunctionEnergy = 0;
    }

    double localEnergy = m_system->getHamiltonian()->
                         computeLocalEnergy(m_system->getParticles());

    double waveFunction = m_system->getWaveFunction()->evaluate(m_system->getParticles());

    double waveFunctionDerivative = m_system->getWaveFunction()->computeAlphaDerivative(m_system->getParticles());

    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergySquared += localEnergy*localEnergy;
    m_cumulativeWaveFunctionDerivative += waveFunctionDerivative / waveFunction;
    m_cumulativeWaveFunctionEnergy += (waveFunctionDerivative / waveFunction) * localEnergy;

    // store energies or positions
    if (writeEnergiesToFile) { writeToFile(localEnergy); }
    if (writePositionsToFile) { writeToFile(); }

    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    int     as = m_system->getNumberOfAcceptedSteps();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << std::setprecision(10) << m_energy << endl;
    cout << " Standard deviation : " << std::setprecision(10) << m_standardDeviation << endl;
    cout << " Acceptance rate : " << as / (double) ms*(1 - ef) << endl;
    cout << endl;
}

void Sampler::writeToFile(double localEnergy) {
    if (m_stepNumber == 0) {
        m_outFile.open("energyTest.dat", std::ios::out | std::ios::trunc);
        m_outFile.close();
    }

    // write local energy to file to do blocking in python
    m_outFile.open("energyTest.dat", std::ios::out | std::ios::app);
    m_outFile << std::setprecision(10) << localEnergy << endl;
    m_outFile.close();


}

void Sampler::writeToFile() {
    if (m_stepNumber == 0) {
        m_outFile2.open("positionsInteraction6.dat", std::ios::out | std::ios::trunc);
        m_outFile2.close();
    }

    // write local energy to file to do blocking in python
    m_outFile2.open("positionsInteraction6.dat", std::ios::out | std::ios::app);
    std::vector<Particle*> particles = m_system->getParticles();
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            m_outFile2 << particles[i]->getPosition()[dim] << "  ";
        }
        m_outFile2 << endl;
    }

    m_outFile2.close();


}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities

    m_energy                    = m_cumulativeEnergy / (double) m_stepNumber;
    m_standardDeviation         = sqrt(std::abs(m_cumulativeEnergySquared / (double) m_stepNumber - m_energy*m_energy));
    m_waveFunctionDerivative    = m_cumulativeWaveFunctionDerivative / (double) m_stepNumber;
    m_waveFunctionEnergy        = m_cumulativeWaveFunctionEnergy / (double) m_stepNumber;
}

void Sampler::clean() {
    m_energy = 0;
    m_waveFunctionDerivative = 0;
    m_waveFunctionEnergy = 0;
    m_cumulativeEnergy = 0;
    m_cumulativeEnergySquared = 0;
    m_cumulativeWaveFunctionDerivative = 0;
    m_cumulativeWaveFunctionEnergy = 0;
    m_stepNumber = 0;
}
