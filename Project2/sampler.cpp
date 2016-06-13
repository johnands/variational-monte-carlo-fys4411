#include <iostream>
#include <mpi.h>
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

struct Pos {
    double x, y;
} pos;

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_firstStep) {
        m_energy = 0;
        m_energySquared = 0;
        m_kineticEnergy = 0;
        m_potentialEnergy = 0;
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeKineticEnergy = 0;
        m_cumulativePotentialEnergy = 0;
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

        if ( m_system->getWriteEnergiesToFile() ) {
            openEnergyFile();
        }
        if ( m_system->getWritePositionsToFile() ) {
            openPositionFile();
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

    m_cumulativeEnergy          += m_localEnergy;
    m_cumulativeEnergySquared   += m_localEnergy*m_localEnergy;
    m_cumulativeKineticEnergy   += m_system->getHamiltonian()->getKineticEnergy();
    m_cumulativePotentialEnergy += m_system->getHamiltonian()->getPotentialEnergy();

    // store energies or positions
    //if (m_system->getWriteEnergiesToFile())   { writeToFile(m_localEnergy); }
    //if (m_system->getWritePositionsToFile())  { writeToFile(); }

    if ( m_energyFileBinary.is_open() ) {
        m_energyFileBinary.write(reinterpret_cast<const char*>(&m_localEnergy), sizeof(double));
    }

    if ( m_positionFileBinary.is_open() ) {
        for (int i=0; i < m_system->getNumberOfParticles(); i++){
            pos.x = m_system->getParticles()[i]->getPosition()[0];
            pos.y = m_system->getParticles()[i]->getPosition()[1];
            m_positionFileBinary.write(reinterpret_cast<char*>(&pos), sizeof(Pos));
        }
    }

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
    cout << " Kinetic energy : " << std::setprecision(10) << m_kineticEnergy << endl;
    cout << " Potential energy : " << std::setprecision(10) << m_potentialEnergy << endl;
    cout << " Variance : " << std::setprecision(10) << m_variance << endl;
    cout << " Acceptance rate : " << m_acceptanceRate << endl;
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

    m_energy                    = m_cumulativeEnergy            / (double) m_stepNumber;
    m_energySquared             = m_cumulativeEnergySquared     / (double) m_stepNumber;
    m_kineticEnergy             = m_cumulativeKineticEnergy     / (double) m_stepNumber;
    m_potentialEnergy           = m_cumulativePotentialEnergy   / (double) m_stepNumber;

    m_variance                  = (m_energySquared - m_energy*m_energy) / (double) m_stepNumber;

    for (int i=0; i < m_numberOfParameters; i++) {
        m_waveFunctionDerivative[i]    = m_cumulativeWaveFunctionDerivative[i] / (double) m_stepNumber;
        m_waveFunctionEnergy[i]        = m_cumulativeWaveFunctionEnergy[i] / (double) m_stepNumber;
    }

    m_acceptanceRate = m_system->getNumberOfAcceptedSteps() /
                       (double) m_system->getNumberOfMetropolisSteps();

    if ( m_system->getWriteEnergiesToFile() )   { closeEnergyFile(); }
    if ( m_system->getWritePositionsToFile() )  { closePositionFile(); }

    if ( m_system->getParallel() ) {
        double reducedEnergy = 0;
        double reducedEnergySquared = 0;
        double reducedKineticEnergy = 0;
        double reducedPotentialEnergy = 0;
        double reducedAcceptanceRate = 0;

        //MPI_Reduce(&m_cumulativeEnergy, &reducedEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m_energy, &reducedEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m_energySquared, &reducedEnergySquared,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m_kineticEnergy, &reducedKineticEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m_potentialEnergy, &reducedPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m_acceptanceRate, &reducedAcceptanceRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        int numberOfProcessors = m_system->getSize();
        if (m_system->getRank() == 0){
            m_energy            = reducedEnergy          / numberOfProcessors;
            m_kineticEnergy     = reducedKineticEnergy   / numberOfProcessors;
            m_potentialEnergy   = reducedPotentialEnergy / numberOfProcessors;
            m_variance          = (reducedEnergySquared  / numberOfProcessors - m_energy*m_energy) /
                                  (m_stepNumber*numberOfProcessors);
            m_acceptanceRate    = reducedAcceptanceRate  / numberOfProcessors;
        }
    }
}

void Sampler::openEnergyFile() {

    std::sprintf( m_energyFileName, "dataFiles/energiesN%dw%dMs%d.bin",
                  (int) m_system->getNumberOfParticles(),
                  (int) (m_system->getWaveFunction()->getOmega()*100),
                  (int) log10(m_numberOfMetropolisSteps) );

    m_energyFileBinary.open(m_energyFileName, std::ios::out);
    if ( m_energyFileBinary.is_open() ) {
    }
}

void Sampler::openPositionFile() {

    std::sprintf( m_positionFileName, "dataFiles/positionsN%dw%dMs%d.bin",
                  (int) m_system->getNumberOfParticles(),
                  (int) (m_system->getWaveFunction()->getOmega()*100),
                  (int) log10(m_numberOfMetropolisSteps) );

    m_positionFileBinary.open(m_positionFileName, std::ios::out | std::ios::binary);
}

void Sampler::closeEnergyFile() {

    m_energyFileBinary.close();
}

void Sampler::closePositionFile() {

    m_positionFileBinary.close();
}

void Sampler::clean() {
    m_energy = 0;
    m_localEnergy = 0;
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_cumulativeEnergy = 0;
    m_cumulativeEnergySquared = 0;
    m_cumulativeKineticEnergy = 0;
    m_cumulativePotentialEnergy = 0;
    m_stepNumber = 0;
    m_acceptanceRate = 0;
    m_firstStep = true;
    m_alphaDerivative = 0;
    m_betaDerivative = 0;
    for (int i=0; i < m_numberOfParameters; i++) {
        m_waveFunctionDerivative[i] = 0;
        m_waveFunctionEnergy[i] = 0;
        m_cumulativeWaveFunctionDerivative[i] = 0;
        m_cumulativeWaveFunctionEnergy[i] = 0;
    }
}
