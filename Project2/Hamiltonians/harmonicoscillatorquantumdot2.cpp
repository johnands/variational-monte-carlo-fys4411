#include "harmonicoscillatorquantumdot2.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillatorQuantumDot2::HarmonicOscillatorQuantumDot2(class System* system, double omega,
                                                             bool useNumerical, bool useInteraction) :
        Hamiltonian(system, useNumerical, useInteraction) {
    assert(m_system->getNumberOfParticles() == 2);
    assert(m_system->getNumberOfDimensions() == 2);
    m_omega = omega;
}


double HarmonicOscillatorQuantumDot2::computeAnalyticalKineticEnergy(std::vector<Particle*> particles) {

    return -0.5*m_system->getWaveFunction()->computeLaplacian(particles);
}

double HarmonicOscillatorQuantumDot2::computePotentialEnergy(std::vector<Particle*> particles) {

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double r1_2 = x1*x1 + y1*y1;
    double r2_2 = x2*x2 + y2*y2;
    double r12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));

    if (m_useInteraction) {
        return 0.5*m_omega*m_omega*(r1_2 + r2_2) + 1.0/r12;
    }
    else {
        return 0.5*m_omega*m_omega*(r1_2 + r2_2);
    }
}
