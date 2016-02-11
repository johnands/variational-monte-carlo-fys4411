#include "harmonicoscillator.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
    Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computePotentialEnergy(std::vector<Particle*> particles) {
    // compute potential energy
    double potentialEnergy = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        double r2 = 0;
        for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
            r2 += pow(particles[i]->getPosition()[j], 2);
        }
        potentialEnergy += r2;
    }
    potentialEnergy *= 0.5*m_omega*m_omega;
    return potentialEnergy;
}

/*double HarmonicOscillator::computeKineticEnergy(std::vector<Particle*> particles) {
    double kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
    return kineticEnergy;
}*/

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    return computePotentialEnergy(particles) + computeKineticEnergy(particles);
}


/* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     *

// compute potential energy
double potentialEnergy = 0;
for (int i=0; i < m_system->getNumberOfParticles(); i++) {
    double r2 = 0;
    for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
        r2 += pow(particles[i]->getPosition()[j], 2);
    }
    potentialEnergy += r2;
}
potentialEnergy *= 0.5*m_omega*m_omega;


double kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);

//double localEnergy = kineticEnergy / m_system->getWaveFunction()->evaluate(particles) + potentialEnergy;
double localEnergy = kineticEnergy + potentialEnergy;
return localEnergy;
}*/

