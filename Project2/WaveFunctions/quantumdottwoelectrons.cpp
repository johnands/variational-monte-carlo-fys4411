#include "quantumdottwoelectrons.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

QuantumDotTwoElectrons::QuantumDotTwoElectrons(System* system, double alpha, double beta, double omega, double a) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_omega = omega;
    m_a = a;
}

double QuantumDotTwoElectrons::evaluate(std::vector<Particle*> particles) {
    // return value of wavefunction

    double alpha = m_parameters[0];
    double beta  = m_parameters[1];

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double r1_2 = x1*x1 + y1*y1;
    double r2_2 = x2*x2 + y2*y2;
    double r12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));

    return exp(-0.5*alpha*m_omega*(r1_2 + r2_2))*exp(r12/(1 + beta*r12));
}

double QuantumDotTwoElectrons::computeLaplacian(std::vector<Particle *> particles) {
    // Calculate laplacian of trial wavefunction divided by trial wavefunction

    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double r1_2 = x1*x1 + y1*y1;
    double r2_2 = x2*x2 + y2*y2;
    double r12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    double K = 1.0/(1 + beta*r12)*(1 + beta*r12);

    double term1 = 2*m_a*(r12*(m_a + beta) + 1)*K*K / r12;
    double term2 = -m_a*(r1_2 + r2_2 - 2*r12*r12)*K / r12;
    double term3 = 2*alpha*m_omega*(alpha*m_omega*(r1_2 + r2_2) - 4);

    return term1 + term2 + term3;
}

std::vector<double> QuantumDotTwoElectrons::computeGradient(std::vector<Particle *> particles) {
    // calculate gradient of trial wavefunction divided by trial wavefunction
    // used to calculate drift velocity / quantum force

    std::vector<double> gradient(4);

    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double dx = x1 - x2;
    double dy = y1 - y2;

    double r12 = sqrt(dx*dx + dy*dy);

    double K1 = alpha*m_omega;
    double K2 = m_a / (r12*(1 + beta*r12)*(1 + beta*r12));

    gradient[0] = -K1*x1 + K2*dx;
    gradient[1] = -K1*y1 + K2*dy;
    gradient[2] = -K1*x2 - K2*dx;
    gradient[3] = -K1*y2 - K2*dy;

    return gradient;
}

std::vector<double> QuantumDotTwoElectrons::computeParametersGradient(std::vector<Particle *> particles) {
    // calculate gradient of wave function w.r.t. the variational parameters
    // divided by the wave function

    std::vector<double> gradient(2);
    double beta = m_parameters[1];

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double r1_2 = x1*x1 + y1*y1;
    double r2_2 = x2*x2 + y2*y2;
    double r12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));

    gradient[0] = -0.5*m_omega*(r1_2 + r2_2);
    gradient[1] = (-m_a*r12*r12) / ((1 + beta*r12)*(1 + beta*r12));

    return gradient;
}
