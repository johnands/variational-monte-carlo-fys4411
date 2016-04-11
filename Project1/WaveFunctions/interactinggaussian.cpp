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
    double r2sum = 0;
    double correlation = 1;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {

        // one-particle wavefunctions
        double r2 = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) { r2 += beta*pow(particles[i]->getPosition()[dim], 2); }
            else { r2 += pow(particles[i]->getPosition()[dim], 2); }
        }
        r2sum += r2;

        // correlation function
        double rij2 = 0;
        for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
            }
            if (sqrt(rij2) > m_a) { correlation *= 1 - m_a / (sqrt(rij2)); }
            else { correlation *= 0; }
        }
    }

    return exp(-alpha*r2sum)*correlation;
}


double InteractingGaussian::computeLaplacian(std::vector<Particle*> particles) {
    // compute Laplacian of trial wavefunction divided by trial wavefunction

    double alpha = m_parameters[0];
    double beta = m_parameters[1];
    int numberOfParticles = m_system->getNumberOfParticles();

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;
    double term4 = 0;

    for (int k=0; k < numberOfParticles; k++) {

        double xk = particles[k]->getPosition()[0];
        double yk = particles[k]->getPosition()[1];
        double zk = particles[k]->getPosition()[2];

        double rk2 = xk*xk + yk*yk + beta*beta*zk*zk;

        term1 += 2*alpha*(2*alpha*rk2 - 2 - beta);

        for (int j=0; j < numberOfParticles; j++) {
            if (k != j) {

                double xkj = xk - particles[j]->getPosition()[0];
                double ykj = yk - particles[j]->getPosition()[1];
                double zkj = zk - particles[j]->getPosition()[2];

                double rkj = sqrt(xkj*xkj + ykj*ykj + zkj*zkj);

                double du_drkj = m_a / (rkj*(rkj - m_a));

                double gradient_dot_rkj =  -2*m_a*(xk*xkj + yk*ykj + beta*zk*zkj);

                term2 += gradient_dot_rkj * du_drkj / rkj;

                double du_drkj2 = (m_a*(m_a - 2*rkj)) / (rkj*rkj*(rkj - m_a)*(rkj - m_a));
                term4 += du_drkj2 + (2.0/rkj)*du_drkj;

                for (int i=0; i < numberOfParticles; i++) {
                    if (k != i) {
                        double xki = xk - particles[i]->getPosition()[0];
                        double yki = yk - particles[i]->getPosition()[1];
                        double zki = zk - particles[i]->getPosition()[2];

                        double rki = sqrt(xki*xki + yki*yki + zki*zki);

                        double du_drki = m_a / (rki*(rki - m_a));

                        double rki_dot_rkj = xki*xkj + yki*ykj + zki*zkj;

                        term3 += rki_dot_rkj * du_drki*du_drkj / (rki*rkj);

                    }
                }
            }
        }
    }

    return term1 + 2*term2 + term3 + term4;
}

std::vector<double> InteractingGaussian::computeGradient(std::vector<Particle*> particles) {

}

double InteractingGaussian::computeAlphaDerivative(std::vector<Particle *> particles) {
    // compute derivative of wavefunction w.r.t. alpha

    double alpha = m_parameters[0];
    double beta  = m_parameters[1];
    double r2sum = 0;
    double correlation = 1;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {

        // one-particle wavefunctions
        double r2 = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) { r2 += beta*pow(particles[i]->getPosition()[dim], 2); }
            else { r2 += pow(particles[i]->getPosition()[dim], 2); }
        }
        r2sum += r2; // one-particle wavefunction

        // correlation function
        double rij2 = 0;
        for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
            }
            if (sqrt(rij2) > m_a) { correlation *= 1 - m_a / (sqrt(rij2)); }
            else { correlation *= 0; }
        }

    }

    return -r2sum*exp(-alpha*r2sum)*correlation;
}
