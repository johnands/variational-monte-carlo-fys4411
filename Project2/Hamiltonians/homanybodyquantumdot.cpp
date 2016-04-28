#include "homanybodyquantumdot.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HOManyBodyQuantumDot::HOManyBodyQuantumDot(class System* system, double omega, bool useNumerical) :
        Hamiltonian(system, useNumerical) {
    assert(m_system->getNumberOfDimensions() == 2);
    m_omega = omega;
}


double HOManyBodyQuantumDot::computeAnalyticalKineticEnergy(std::vector<Particle*> particles) {

    return -0.5*m_system->getWaveFunction()->computeLaplacian(particles);
}

double HOManyBodyQuantumDot::computePotentialEnergy(std::vector<Particle*> particles) {

    double singleParticle = 0;
    double interacting = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {

        // single-particle potential energy
        double r2 = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            r2 += pow(particles[i]->getPosition()[dim], 2);
        }
        singleParticle += r2;

        // interacting potential energy
        for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
            double rij2 = 0;
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
            }
            interacting += 1.0 / sqrt(rij2);
        }
    }

    singleParticle *= 0.5*m_omega*m_omega;

    return singleParticle + interacting;
}
