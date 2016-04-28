#include "manybodyquantumdot.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <armadillo>
#include <algorithm>

using std::cout;
using std::endl;

ManyBodyQuantumDot::ManyBodyQuantumDot(System* system, double alpha, double beta, double omega, double a) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_omega = omega;
    m_omegaSqrt = sqrt(omega);
    m_a = a;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfParticlesHalf = m_numberOfParticles / 2;
    setUpSlater();
}

double ManyBodyQuantumDot::computeRatio(std::vector<Particle *> particles, int i) {
    // compute ratio used in Metropolis algorithm

    // store particle number for later use in Gradient and Laplacian
    m_i = i;

    double ratioSD = 0;

    // get new position of particle i
    double xNew = particles[i]->getNewPosition()[0];
    double yNew = particles[i]->getNewPosition()[1];

    // spin-up slater
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratioSD += singleParticleWaveFunctions(nx, ny, xNew, yNew)*m_slaterSpinUpInverse(j,i);
        }
    }
    // spin-down slater
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratioSD += singleParticleWaveFunctions(nx, ny, xNew, yNew) *
                       m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
        }
    }

    // compute jastrow factor
    double exponent = 0;
    double beta = m_parameters[1];
    for (int j=0; j < m_numberOfParticles; j++) {
        if (j != i) {
            double r_jiNew = 0;
            double r_jiOld = 0;
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                r_jiNew += (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]) *
                           (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]);
                r_jiOld += (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]) *
                           (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]);
            }
            exponent += r_jiNew / (1 + beta*r_jiNew);
            exponent -= r_jiOld / (1 + beta*r_jiOld);
        }
    }
    double ratioJastrow = exp(m_a*exponent);

    m_ratio = ratioSD*ratioJastrow;
    return m_ratio;
}

void ManyBodyQuantumDot::setUpSlater() {
    // set up Slater matrix and invert

    m_quantumNumbers = arma::zeros<arma::mat>(10, 2);

    // quantum numbers for up to 20 particles
    m_quantumNumbers(0,0) = 0; m_quantumNumbers(0,1) = 0;
    m_quantumNumbers(1,0) = 1; m_quantumNumbers(1,1) = 0;
    m_quantumNumbers(2,0) = 0; m_quantumNumbers(2,1) = 1;
    m_quantumNumbers(3,0) = 2; m_quantumNumbers(3,1) = 0;
    m_quantumNumbers(4,0) = 1; m_quantumNumbers(4,1) = 1;
    m_quantumNumbers(5,0) = 0; m_quantumNumbers(5,1) = 2;
    m_quantumNumbers(6,0) = 3; m_quantumNumbers(6,1) = 0;
    m_quantumNumbers(7,0) = 2; m_quantumNumbers(7,1) = 1;
    m_quantumNumbers(8,0) = 1; m_quantumNumbers(8,1) = 2;
    m_quantumNumbers(9,0) = 0; m_quantumNumbers(9,1) = 3;

    // fill the matrix elements, which are the one-particle harm. osc. wavefunctions
    m_slaterSpinUp = arma::zeros<arma::mat>(m_numberOfParticlesHalf, m_numberOfParticlesHalf);
    m_slaterSpinDown = arma::zeros<arma::mat>(m_numberOfParticlesHalf, m_numberOfParticlesHalf);

    for (int i=0; i < m_numberOfParticlesHalf; i++) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double xUp = m_system->getParticles()[i]->getPosition()[0];
            double yUp = m_system->getParticles()[i]->getPosition()[1];
            double xDown = m_system->getParticles()[i+m_numberOfParticlesHalf]->getPosition()[0];
            double yDown = m_system->getParticles()[i+m_numberOfParticlesHalf]->getPosition()[1];
            m_slaterSpinUp(i,j) = singleParticleWaveFunctions(nx, ny, xUp, yUp);
            m_slaterSpinDown(i,j) = singleParticleWaveFunctions(nx, ny, xDown, yDown);
        }
    }

    m_slaterSpinUpInverse = arma::inv(m_slaterSpinUp);
    m_slaterSpinDownInverse = arma::inv(m_slaterSpinDown);
}

double ManyBodyQuantumDot::singleParticleWaveFunctions(int nx, int ny, double x, double y) {
    // evaluate single particle wave functions

    return hermitePolynomials(nx, x)*hermitePolynomials(ny, y) *
           exp(-0.5*m_omega*(x*x + y*y));
}

std::vector<double> ManyBodyQuantumDot::singleParticleWFGradient(int nx, int ny, double x, double y) {

    std::vector<double> gradient(2);
    gradient[0] =  exp(-0.5*m_omega*(x*x + y*y)) * hermitePolynomials(ny, y) *
                   ( hermitePolynomialsDerivative1(nx, x) - hermitePolynomials(nx, x)*m_omega*x);
    gradient[1] =  exp(-0.5*m_omega*(x*x + y*y)) * hermitePolynomials(nx, x) *
                   ( hermitePolynomialsDerivative1(ny, y) - hermitePolynomials(ny, y)*m_omega*y);
    return gradient;
}

double ManyBodyQuantumDot::singleParticleWFLaplacian(int nx, int ny, double x, double y) {

    return exp(-0.5*m_omega*(x*x + y*y)) *
           ( - 2*m_omega*x*hermitePolynomials(ny, y)*hermitePolynomialsDerivative1(nx, x)
             - 2*m_omega*y*hermitePolynomials(nx, x)*hermitePolynomialsDerivative1(ny, y)
             + m_omega*hermitePolynomials(nx, x)*hermitePolynomials(ny, y) * (m_omega*(x*x + y*y) - 2)
             + hermitePolynomials(ny, y)*hermitePolynomialsDerivative2(nx, x)
             + hermitePolynomials(nx, x)*hermitePolynomialsDerivative2(ny, y) );
}


double ManyBodyQuantumDot::hermitePolynomials(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 1;
    }

    else if (energyLevel == 1) {
        return 2*m_omegaSqrt*position;
    }

    else if (energyLevel == 2) {
        return 4*m_omega*position*position;
    }

    else if (energyLevel == 3) {
        return 8*m_omega*m_omegaSqrt*position*position*position - 12*m_omegaSqrt*position;
    }

    else if (energyLevel == 4) {
        return 16*m_omega*m_omega*position*position*position*position -
               48*m_omega*position*position + 12;
    }

    else {
        cout << "Energy level should not exceed n = 4" << endl;
        exit(0);
    }
}

double ManyBodyQuantumDot::hermitePolynomialsDerivative1(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 0;
    }

    else if (energyLevel == 1) {
        return 2;
    }

    else if (energyLevel == 2) {
        return 8*m_omegaSqrt*position;
    }

    else if (energyLevel == 3) {
        return 24*m_omega*position*position - 12;
    }

    else if (energyLevel == 4) {
        return 64*m_omega*m_omegaSqrt*position*position*position -
               96*m_omegaSqrt*position;
    }

    else {
        cout << "Energy level should not exceed n = 4" << endl;
        exit(0);
    }
}

double ManyBodyQuantumDot::hermitePolynomialsDerivative2(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 0;
    }

    else if (energyLevel == 1) {
        return 0;
    }

    else if (energyLevel == 2) {
        return 8;
    }

    else if (energyLevel == 3) {
        return 48*m_omegaSqrt*position;
    }

    else if (energyLevel == 4) {
        return 192*m_omega*position*position - 96;
    }

    else {
        cout << "Energy level should not exceed n = 4" << endl;
        exit(0);
    }
}

void ManyBodyQuantumDot::updateRowSlater(std::vector<Particle*> particles, int i) {
    // update row corresponding to particle i in Slater matrix

    // maybe update Slater and inverse Slater here??
    // need both current Slater and new slater to compute new inverse Slater
    // or make matrices NewSlaterUp and NewSlaterDown that are
    // temporary so that I have access to both new and old Slater matrices
    // then I must set SlaterUp = NewSlaterUp before new cycle

    // get new position of particle i
    double xNew = particles[i]->getNewPosition()[0];
    double yNew = particles[i]->getNewPosition()[1];

    // update either spin-up matrix or spin-down matrix
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            // compute new inverse here?
            m_slaterSpinUp(i,j) = singleParticleWaveFunctions(nx, ny, xNew, yNew);
        }
    }
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            m_slaterSpinDown(i-m_numberOfParticlesHalf,j) = singleParticleWaveFunctions(nx, ny, xNew, yNew);
        }
    }
    updateSlaterInverse(particles, i);
}

void ManyBodyQuantumDot::updateSlaterInverse(std::vector<Particle*> particles, int i) {
    // update row corresponding to particle i in inverse Slater matrix

    // get new position of particle i
    //double xNew = particles[i]->getNewPosition()[0];
    //double yNew = particles[i]->getNewPosition()[1];

    // compute Sj for all columns except column i
    std::vector<double> S(m_numberOfParticlesHalf);
    std::fill(S.begin(), S.end(), 0);
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            if (j != i) {
                for (int l=0; l < m_numberOfParticlesHalf; l++) {
                    S[j] += m_slaterSpinUp(i,l)*m_slaterSpinUpInverse(l,j);
                }
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinUpInverse(k,j) = m_slaterSpinUpInverse(k,j) -
                                                 S[j]*m_slaterSpinUpInverse(k,i) / m_ratio;
                }
            }
            else {
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinUpInverse(k,i) = m_slaterSpinUpInverse(k,i) / m_ratio;
                }
            }
        }
    }
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            // need to adjust i. for N = 6, i=4 is equal to i-3=1 in spin down matrix
            int iDown = i - m_numberOfParticlesHalf;
            if (j != i) {
                for (int l=0; l < m_numberOfParticlesHalf; l++) {
                    S[j] += m_slaterSpinDown(iDown,l)*m_slaterSpinDownInverse(l,j);
                }
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinDownInverse(k,j) = m_slaterSpinDownInverse(k,j) -
                                                 S[j]*m_slaterSpinDownInverse(k,iDown)
                                                 / m_ratio;
                }
            }
            else {
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinDownInverse(k,iDown) = m_slaterSpinDownInverse(k,iDown) / m_ratio;
                }
            }
        }
    }
}

double ManyBodyQuantumDot::evaluate(std::vector<Particle*> particles) {

}

double ManyBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles) {
    // Calculate laplacian of Slater matrix

    double laplacianSlater = 0;

    // spin-up slater
    if (m_i < m_numberOfParticlesHalf) {
        // find laplacian for all particles
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double x = particles[i]->getPosition()[0];
                double y = particles[i]->getPosition()[1];
                // singleParticleWFLaplacian returns the full Laplacian of the specific
                // single-particle wave function
                laplacianSlater += singleParticleWFLaplacian(nx, ny, x, y) * m_slaterSpinUpInverse(j,i);
            }
        }
    }

    // spin-down slater
    else {
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double x = particles[i]->getPosition()[0];
                double y = particles[i]->getPosition()[1];
                laplacianSlater += singleParticleWFLaplacian(nx, ny, x, y) * m_slaterSpinDownInverse(j,i);
            }
        }
    }

    // compute jastrow factor
    double term1 = 0;
    double term2 = 0;
    double beta = m_parameters[1];
    for (int i=0; i < m_numberOfParticles; i++) {

        for (int k=0; k < m_numberOfParticles; k++) {

            if (k != i) {
                double x_i = particles[i]->getPosition()[0];
                double y_i = particles[i]->getPosition()[1];

                double x_k = particles[k]->getPosition()[0];
                double y_k = particles[k]->getPosition()[1];

                double r_ik = (x_i - x_k)*(x_i - x_k) + (y_i - y_k)*(y_i - y_k);
                double factor1 = 1.0 / (1 + beta*r_ik);

                term2 += (2*m_a*factor1*factor1 / r_ik) - 2*m_a*beta*factor1*factor1*factor1;

                for (int j=0; j < m_numberOfParticles; j++) {

                    if (j != m_i) {
                        double x_j = particles[j]->getPosition()[0];
                        double y_j = particles[j]->getPosition()[1];

                        double r_ij = (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j);
                        double dot_prod = (x_i - x_j)*(x_i - x_k) + (y_i - y_j)*(y_i - y_j);
                        double factor2 = 1.0 / (1 + beta*r_ij);

                        term1 += (dot_prod*m_a*m_a*factor1*factor1*factor2*factor2) / (r_ij*r_ik);
                    }
                }
            }
        }
    }

    return laplacianSlater*(term1 + term2);
}

std::vector<double> ManyBodyQuantumDot::computeGradient(std::vector<Particle *> particles) {
    // calculate gradient of trial wavefunction divided by trial wavefunction
    // used to calculate drift velocity / quantum force

    std::vector<double> gradient(2);

    std::vector<double> ratioSD(2);
    ratioSD[0] = 0; ratioSD[1] = 0;

    // get new position of particle i
    double x_i = particles[m_i]->getNewPosition()[0];
    double y_i = particles[m_i]->getNewPosition()[1];

    // spin-up slater
    if (m_i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x_i, y_i);
            ratioSD[0] += grad[0] * m_slaterSpinUpInverse(j,m_i);
            ratioSD[1] += grad[1] * m_slaterSpinUpInverse(j,m_i);
        }
    }
    // spin-down slater
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x_i, y_i);
            ratioSD[0] += grad[0] * m_slaterSpinDownInverse(j,m_i-m_numberOfParticlesHalf);
            ratioSD[1] += grad[1] * m_slaterSpinDownInverse(j,m_i-m_numberOfParticlesHalf);
        }
    }

    // compute jastrow factor
    std::vector<double> ratioJastrow(2);
    ratioJastrow[0] = 0; ratioJastrow[1] = 0;
    double beta = m_parameters[1];
    for (int j=0; j < m_numberOfParticles; j++) {
        if (j != m_i) {
            double x_j = particles[j]->getPosition()[0];
            double y_j = particles[j]->getPosition()[1];
            double r_ij = (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j);
            double factor = 1.0 / ( r_ij*(1 + beta*r_ij)*(1 + beta*r_ij) );

            ratioJastrow[0] += (x_i - x_j)*m_a*factor;
            ratioJastrow[1] += (y_i - y_j)*m_a*factor;
        }
    }

    gradient[0] = ratioSD[0]*ratioJastrow[0];
    gradient[1] = ratioSD[1]*ratioJastrow[1];

    return gradient;
}

std::vector<double> ManyBodyQuantumDot::computeParametersGradient(std::vector<Particle *> particles) {
    // calculate gradient of wave function w.r.t. the variational parameters
    // divided by the wave function


}
