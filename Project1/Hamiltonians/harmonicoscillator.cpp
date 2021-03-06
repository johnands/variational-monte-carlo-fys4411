#include "harmonicoscillator.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, bool useNumerical) :
    Hamiltonian(system, useNumerical) {
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

double HarmonicOscillator::computeAnalyticalKineticEnergy(std::vector<Particle*> particles) {
    double kineticEnergy = -0.5*m_system->getWaveFunction()->computeLaplacian(particles);
    return kineticEnergy;
}


