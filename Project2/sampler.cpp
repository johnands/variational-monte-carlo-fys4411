#include <iostream>
//#include <mpi.h>
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
    m_numberOfParameters = m_system->getWaveFunction()->getNumberOfParameters();
    m_waveFunctionDerivative.resize(m_numberOfParameters);
    m_waveFunctionEnergy.resize(m_numberOfParameters);
    m_cumulativeWaveFunctionDerivative.resize(m_numberOfParameters);
    m_cumulativeWaveFunctionEnergy.resize(m_numberOfParameters);
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_firstStep) {
        m_energy = 0;
        m_energySquared = 0;
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        for (int i=0; i < m_numberOfParameters; i++) {
            m_waveFunctionDerivative[i] = 0;
            m_waveFunctionEnergy[i] = 0;
            m_cumulativeWaveFunctionDerivative[i] = 0;
            m_cumulativeWaveFunctionEnergy[i] = 0;
        }

        m_localEnergy           = m_system->getHamiltonian()->
                                  computeLocalEnergy(m_system->getParticles());
        if (m_system->getOptimizeParameters()) {
            m_parametersGradient = m_system->getWaveFunction()->
                                   computeParametersGradient(m_system->getParticles());
        }
        m_firstStep = false;
    }

    // energy and derivative w.r.t varational parameters is the same if step was not accepted
    if (acceptedStep) {
        m_localEnergy = m_system->getHamiltonian()->
                        computeLocalEnergy(m_system->getParticles());

        if (m_system->getOptimizeParameters()) {
            m_parametersGradient = m_system->getWaveFunction()->
                                   computeParametersGradient(m_system->getParticles());
        }
    }

    // sample
    if (m_system->getOptimizeParameters()) {
        for (int i=0; i < m_numberOfParameters; i++) {
            m_cumulativeWaveFunctionDerivative[i] += m_parametersGradient[i];
            m_cumulativeWaveFunctionEnergy[i] += m_parametersGradient[i] * m_localEnergy;
        }
    }

    m_cumulativeEnergy  += m_localEnergy;
    m_cumulativeEnergySquared += m_localEnergy*m_localEnergy;

    // store energies or positions
    if (m_system->getWriteEnergiesToFile())   { writeToFile(m_localEnergy); }
    if (m_system->getWritePositionsToFile())  { writeToFile(); }

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
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << std::setprecision(10) << m_energy << endl;
    cout << " Variance : " << std::setprecision(10) << m_variance << endl;
    cout << " Acceptance rate : " << m_acceptanceRate << endl;
    cout << " Number of accepted steps " << as << endl;
    cout << endl;
}

void Sampler::writeToFile(double localEnergy) {
    if (m_stepNumber == 0) {
        m_outFile.open("energyTest2.dat", std::ios::out | std::ios::trunc);
        m_outFile.close();
    }

    // write local energy to file to do blocking in python
    m_outFile.open("energyTest2.dat", std::ios::out | std::ios::app);
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
    m_energySquared             = m_cumulativeEnergySquared / (double) m_stepNumber;
    m_variance                  = (m_energySquared - m_energy*m_energy) / (double) m_stepNumber;

    for (int i=0; i < m_numberOfParameters; i++) {
        m_waveFunctionDerivative[i]    = m_cumulativeWaveFunctionDerivative[i] / (double) m_stepNumber;
        m_waveFunctionEnergy[i]        = m_cumulativeWaveFunctionEnergy[i] / (double) m_stepNumber;
    }

    m_acceptanceRate = m_system->getNumberOfAcceptedSteps() /
                       (double) m_system->getNumberOfMetropolisSteps();

    if (m_system->getParallel()) {
        double reducedEnergy = 0;
        double reducedVariance = 0;
        double reducedAcceptanceRate = 0;
        //MPI_Reduce(&m_energy, &reducedEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //MPI_Reduce(&m_variance, &reducedVariance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //MPI_Reduce(&m_acceptanceRate, &reducedAcceptanceRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        int size = m_system->getSize();
        if (m_system->getRank() == 0){
            m_energy = reducedEnergy / size;
            m_variance = reducedVariance / size;
            m_acceptanceRate = reducedAcceptanceRate / size;
        }
    }
}

void Sampler::clean() {
    m_energy = 0;
    m_cumulativeEnergy = 0;
    m_cumulativeEnergySquared = 0;
    m_stepNumber = 0;
    for (int i=0; i < m_numberOfParameters; i++) {
        m_waveFunctionDerivative[i] = 0;
        m_waveFunctionEnergy[i] = 0;
        m_cumulativeWaveFunctionDerivative[i] = 0;
        m_cumulativeWaveFunctionEnergy[i] = 0;
    }
}
