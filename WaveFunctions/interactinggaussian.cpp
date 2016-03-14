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

    double beta = m_parameters[1];
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;
    double term4 = 0;
    for (int k=0; k < numberOfParticles; k++) {

        // term 1
        double r2 = 0;
        for (int dim=0; dim < numberOfDimensions; dim++) {
            if (dim == 2) { r2 += pow(beta*particles[k]->getPosition()[dim], 2); }
            else { r2 += pow(particles[k]->getPosition()[dim], 2); }
        }
        term1 += 2*m_a*(2*m_a*r2 - 2 - beta);

        // term 2
        std::vector<double> r_kj;
        for (int dim=0; dim < numberOfDimensions; dim++) {
            r_kj.push_back(0.0);
        }
        for (int j=0; j< numberOfParticles; j++) {
            if (k != j) {
                double r2 = 0;
                for (int dim=0; dim < numberOfDimensions; dim++) {
                    double dim_kj = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
                    r_kj[dim] += dim_kj;      // add to vector
                    r2 += dim_kj*dim_kj;      // distance
                }
                for (int dim=0; dim < numberOfDimensions; dim++) {
                    r_kj[dim] *= m_a / (r2*(sqrt(r2) - m_a));
                }
            }
        }
        for (int dim=0; dim < numberOfDimensions; dim++) {
            if (dim == 2) { term2 += -2*m_a*beta*particles[k]->getPosition()[dim]*r_kj[dim]; }
            else { term2 += -2*m_a*particles[k]->getPosition()[dim]*r_kj[dim]; }
        }

        // term 3
        for (int i=0; i< numberOfParticles; i++) {
            for (int j=0; j< numberOfParticles; j++) {
                if (k != i && k != j) {
                    double rki2 = 0;
                    double rkj2 = 0;
                    double dot = 0;
                    for (int dim=0; dim < numberOfDimensions; dim++) {
                        double dim_ki = particles[k]->getPosition()[dim] - particles[i]->getPosition()[dim];
                        double dim_kj = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
                        rki2 += dim_ki*dim_ki;     // r_ki^2
                        rkj2 += dim_kj*dim_kj;     // r_kj^2
                        dot += dim_ki*dim_kj;      // r_ki*r_kj
                    }
                    for (int dim=0; dim < numberOfDimensions; dim++) {
                        term3 += (dot*m_a*m_a) /
                                 (sqrt(rki2)*sqrt(rkj2)*rkj2*(sqrt(rkj2) - m_a)*(sqrt(rkj2) - m_a));
                    }
                }
            }
        }

        // term 4
        for (int j=0; j< numberOfParticles; j++) {
            if (k != j) {
                double rkj2 = 0;
                for (int dim=0; dim < numberOfDimensions; dim++) {
                    double dim_kj = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
                    rkj2 += dim_kj*dim_kj;     // r_kj^2
                }
                for (int dim=0; dim < numberOfDimensions; dim++) {
                    double first = (m_a*(m_a - 2*sqrt(rkj2))) / (rkj2*(sqrt(rkj2) - m_a)*(sqrt(rkj2) - m_a));
                    double second = (2*m_a) / (rkj2*(sqrt(rkj2) - m_a));
                    term4 += first + second;
                }
            }
        }
    }


    return term1 + term2 + term3 + term4;
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
