#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<Particle*> particles) {
     // return value of product of gaussian one-particle wavefunctions

     double alpha = m_parameters[0];
     double gaussian = 1;
     for (int i=0; i < m_system->getNumberOfParticles(); i++) {
         double r2 = 0;
         for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
             r2 += pow(particles[i]->getPosition()[j], 2);
         }
         gaussian *= exp(-alpha*r2);
     }

     /*
     // 1d, one particle
     double x = particles[0]->getPosition()[0];
     double alpha = m_parameters[0];
     return exp(-alpha*x*x);*/

     return gaussian;
}

double SimpleGaussian::computeLaplacian(std::vector<Particle*> particles) {
    // Calculate double derivative of trial wavefunction divided by trial wavefunction

    double doubleDerivative = 0;
    double alpha = m_parameters[0];
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        double r2 = 0;
        for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
            r2 += pow(particles[i]->getPosition()[j], 2);
            //doubleDerivative += 2*alpha*exp(-alpha*x*x)*(2*alpha*x*x - 1);
        }
        doubleDerivative += 2*alpha*(2*alpha*r2 - m_system->getNumberOfDimensions());
    }
    return doubleDerivative;

    /*// 1d, one particle
    double x = particles[0]->getPosition()[0];
    double alpha = m_parameters[0];
    return 2*alpha*exp(-alpha*x*x)*(2*alpha*x*x - 1);*/
}

std::vector<double> SimpleGaussian::computeGradient(std::vector<Particle*> particles) {
    // calculate gradient of trial wavefunction divided by trial wavefunction
    // used to calculate drift velocity / quantum force

    std::vector<double> firstDerivative;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    firstDerivative.resize(numberOfParticles*numberOfDimensions);

    double alpha = m_parameters[0];

    for (int i=0; i < numberOfParticles; i++) {
        for (int j=0; j < numberOfDimensions; j++) {
        firstDerivative[i+j] += -2*alpha*particles[i]->getPosition()[j];
        }
    }
    return firstDerivative;
}
