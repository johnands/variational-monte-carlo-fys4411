#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <iostream>

Hamiltonian::Hamiltonian(System* system, bool useNumerical, bool useInteraction) {
    m_system = system;
    m_useNumerical = useNumerical;
    m_useInteraction = useInteraction;
}

double Hamiltonian::computeLocalEnergy(std::vector<Particle*> particles) {
    double kinetic = 0;
    if (m_useNumerical) {
        kinetic = computeKineticEnergy(particles);
    } else {
        kinetic = computeAnalyticalKineticEnergy(particles);
    }
    return computePotentialEnergy(particles) + kinetic;
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles) {
    // the kinetic energy is the same for all Hamiltonians
    // it is therefore moved to the super-class

    double laplacian = 0;
    double step = 0.001;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions(); j++) {

            particles[i]->adjustPosition(step, j);
            double plusStep = m_system->getWaveFunction()->evaluate(particles);

            particles[i]->adjustPosition(-2*step, j);
            double minusStep = m_system->getWaveFunction()->evaluate(particles);

            particles[i]->adjustPosition(step, j);
            double zeroStep = m_system->getWaveFunction()->evaluate(particles);

            laplacian += (plusStep - 2*zeroStep + minusStep) / (zeroStep*step*step);
        }
    }
    //std::cout << laplacian << std::endl;
    return -0.5*laplacian;   
}
