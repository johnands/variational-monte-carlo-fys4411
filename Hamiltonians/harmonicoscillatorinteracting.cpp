#include "harmonicoscillatorinteracting.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillatorInteracting::HarmonicOscillatorInteracting(class System* system, double omega, bool useNumerical) :
        Hamiltonian(system, useNumerical) {
    m_omega = omega;
}


double HarmonicOscillatorInteracting::computeAnalyticalKineticEnergy(std::vector<class Particle*> particles) {
    return 0;
}

double HarmonicOscillatorInteracting::computePotentialEnergy(std::vector<class Particle*> particles) {
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
