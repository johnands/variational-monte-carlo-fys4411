#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep, bool writeEnergiesToFile) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    double localEnergy = m_system->getHamiltonian()->                       // this gets the Hamiltonian, which itself is an instance of Hamiltonian
                         computeLocalEnergy(m_system->getParticles());
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergySquared += localEnergy*localEnergy;

    // store energies to do blocking
    if (writeEnergiesToFile) { writeToFile(localEnergy); }

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
    cout << " Energy : " << m_energy << endl;
    cout << " Standard deviation : " << m_standardDeviation << endl;
    cout << " Acceptance rate : " << as / (double) ms*(1 - ef) << endl;
    cout << endl;
}

void Sampler::writeToFile(double localEnergy) {
    if (m_stepNumber == 0) { m_outFile.open("energy.dat", std::ios::out); }

    // write local energy to file to do blocking in python
    m_outFile << localEnergy << endl;


}

void Sampler::closeFile() {
    m_outFile.close();
}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities

    //int numberOfSampledSteps = m_numberOfMetropolisSteps*(1 - m_system->getEquilibrationFraction())-1;
    //m_energy = m_cumulativeEnergy / (double) numberOfSampledSteps;
    m_energy = m_cumulativeEnergy / (double) m_stepNumber;
    m_standardDeviation = sqrt(m_cumulativeEnergySquared / (double) m_stepNumber - m_energy*m_energy);
}
