#include "interactinggaussian.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

InteractingGaussian::InteractingGaussian(System* system, double alpha, double beta, double a) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_a = a;
}

double InteractingGaussian::evaluate(std::vector<Particle*> particles) {
    // return value of product of gaussian one-particle wavefunctions

    double alpha = m_parameters[0];
    double beta  = m_parameters[1];
    double gaussian = 1;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {

        // one-particle wavefunctions
        double r2 = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) { r2 += beta*pow(particles[i]->getPosition()[dim], 2); }
            else { r2 += pow(particles[i]->getPosition()[dim], 2); }
        }
        gaussian *= exp(-alpha*r2); // one-particle wavefunction

        // correlation function
        double rij2 = 0;
        for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
            }
            if (sqrt(rij2) > m_a) { gaussian *= 1 - m_a / (sqrt(rij2)); }
            else { gaussian *= 0; }
        }

    }

    return gaussian;
}



double InteractingGaussian::computeLaplacian(std::vector<Particle*> particles) {
    return 0;
}

std::vector<double> InteractingGaussian::computeGradient(std::vector<Particle*> particles) {

}
