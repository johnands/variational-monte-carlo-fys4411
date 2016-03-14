#include "harmonicoscillatorinteracting.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillatorInteracting::HarmonicOscillatorInteracting(class System* system, double omega, double a, double gamma, bool useNumerical) :
        Hamiltonian(system, useNumerical) {
    m_omega = omega;
    m_gamma = gamma;
    m_a     = a;
}


double HarmonicOscillatorInteracting::computeAnalyticalKineticEnergy(std::vector<Particle*> particles) {
    double kineticEnergy = -0.5*m_system->getWaveFunction()->computeLaplacian(particles);
    return kineticEnergy;
}

double HarmonicOscillatorInteracting::computePotentialEnergy(std::vector<Particle*> particles) {

   double potentialEnergy = 0;
   for (int i=0; i < m_system->getNumberOfParticles(); i++) {

       // single-particle potential energy
       double r2 = 0;
       for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
           if (dim == 2) { r2 += pow(m_gamma*particles[i]->getPosition()[dim], 2); }
           else { r2 += pow(particles[i]->getPosition()[dim], 2); }
       }
       potentialEnergy += r2;

       // interacting potential energy
       for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
           double rij2 = 0;
           for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
               rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
           }
           if (sqrt(rij2) <= m_a) { potentialEnergy += 1e10; std::cout << "yes" << std::endl; }
       }
   }

   potentialEnergy *= 0.5*m_omega*m_omega;

   return potentialEnergy;
}
