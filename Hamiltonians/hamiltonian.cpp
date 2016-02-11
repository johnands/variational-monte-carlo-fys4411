#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles) {
    // the kinetic energy is the same for all Hamiltonians
    // it is therefore moved to the super-class

    double laplacian = 0;
    double step = 0.01;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions(); j++) {

            particles[i]->adjustPosition(step, j);
            double plusStep = m_system->getWaveFunction()->evaluate(particles);

            particles[i]->adjustPosition(-2*step, j);
            double minusStep = m_system->getWaveFunction()->evaluate(particles);

            particles[i]->adjustPosition(step, j);
            double zeroStep = m_system->getWaveFunction()->evaluate(particles);

            laplacian += (plusStep - 2*zeroStep + minusStep) / (step*step);
        }
    }
    return -0.5*laplacian;
}
