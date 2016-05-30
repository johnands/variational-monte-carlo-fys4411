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

ManyBodyQuantumDot::ManyBodyQuantumDot(System* system, double alpha, double beta, double omega) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_omega                 = omega;
    m_omegaSqrt             = sqrt(omega);
    m_omegaAlpha            = omega*alpha;
    m_omegaAlphaSqrt        = sqrt(m_omegaAlpha);
    m_alphaSqrtInv          = 1.0/sqrt(alpha);
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfParticlesHalf = m_numberOfParticles / 2;
    m_gradientUp.resize(2); m_gradientDown.resize(2);
    setUpSlater();
}

double ManyBodyQuantumDot::computeRatio(std::vector<Particle *> particles, int i) {
    // compute ratio used in Metropolis algorithm

    // store particle number for later use in Laplacian
    m_i = i;

    // Slater determinant ratio
    double ratioSD = 0;

    // get new position of chosen particle
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
    double exponent1 = 0;
    double exponent2 = 0;
    double beta = m_parameters[1];
    for (int j=0; j < i; j++) {
        double r_jiNew = 0;
        double r_jiOld = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            r_jiNew += (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]) *
                       (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]);
            r_jiOld += (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]) *
                       (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]);
        }
        r_jiNew = sqrt(r_jiNew);
        r_jiOld = sqrt(r_jiOld);
        exponent1 += ( m_a(i,j)*r_jiNew ) / ( 1 + beta*r_jiNew );
        exponent2 += ( m_a(i,j)*r_jiOld ) / ( 1 + beta*r_jiOld );

    }
    for (int j=i+1; j < m_numberOfParticles; j++) {
        double r_jiNew = 0;
        double r_jiOld = 0;
        for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
            r_jiNew += (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]) *
                       (particles[j]->getPosition()[dim] - particles[i]->getNewPosition()[dim]);
            r_jiOld += (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]) *
                       (particles[j]->getPosition()[dim] - particles[i]->getPosition()[dim]);
        }
        r_jiNew = sqrt(r_jiNew);
        r_jiOld = sqrt(r_jiOld);
        exponent1 += ( m_a(i,j)*r_jiNew ) / ( 1 + beta*r_jiNew );
        exponent2 += ( m_a(i,j)*r_jiOld ) / ( 1 + beta*r_jiOld );

    }

    double ratioJastrow = exp(exponent1 - exponent2);

    m_ratioSD = ratioSD; m_system->setRatioSD(ratioSD);
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

    m_a = arma::zeros<arma::mat>(m_numberOfParticles, m_numberOfParticles);

    // a = 1 for parallel spins, a = 1/3 for antiparallel spins
    if (m_system->getUseJastrow()) {
        for (int i=0; i < m_numberOfParticles; i++) {
            for (int j=0; j < m_numberOfParticles; j++) {
                if (i < m_numberOfParticlesHalf) {
                    if (j < m_numberOfParticlesHalf) {
                        m_a(i,j) = 1.0/3;
                    }
                    else {
                        m_a(i,j) = 1;
                    }
                }
                else {
                    if (j < m_numberOfParticlesHalf) {
                        m_a(i,j) = 1;
                    }
                    else {
                        m_a(i,j) = 1.0/3;
                    }
                }
            }
        }
    }

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

    // invert
    m_slaterSpinUpInverse = m_slaterSpinUp.i();
    m_slaterSpinDownInverse = m_slaterSpinDown.i();

    // write out initial matrices
    /*cout << "Initial Slater Up: " << m_slaterSpinUp << endl;
    cout << "Initial Slater Down: " << m_slaterSpinDown << endl;

    cout << "Initial Slater Up Inverse: " << m_slaterSpinUpInverse << endl;
    cout << "Initial Slater Down Inverse: " << m_slaterSpinDownInverse << endl;*/
}

double ManyBodyQuantumDot::singleParticleWaveFunctions(int nx, int ny, double x, double y) {

    return hermitePolynomials(nx, x)*hermitePolynomials(ny, y) *
           exp(-0.5*m_omegaAlpha*(x*x + y*y));
}

std::vector<double> ManyBodyQuantumDot::singleParticleWFGradient(int nx, int ny, double x, double y) {

    std::vector<double> gradient(2);
    double r2 = x*x + y*y;
    double hermiteX = hermitePolynomials(nx, x);
    double hermiteY = hermitePolynomials(ny, y);
    gradient[0] =  exp(-0.5*m_omegaAlpha*r2) * hermiteY *
                   ( hermitePolynomialsDerivative1(nx, x) - hermiteX*m_omegaAlpha*x);
    gradient[1] =  exp(-0.5*m_omegaAlpha*r2) * hermiteX *
                   ( hermitePolynomialsDerivative1(ny, y) - hermiteY*m_omegaAlpha*y);

    return gradient;
}

double ManyBodyQuantumDot::singleParticleWFLaplacian(int nx, int ny, double x, double y) {
    // return full Laplacian for each single-particle wave function

    double r2 = x*x + y*y;
    double hermiteX = hermitePolynomials(nx, x);
    double hermiteY = hermitePolynomials(ny, y);
    return exp(-0.5*m_omegaAlpha*r2) *
           ( - 2*m_omegaAlpha*x*hermiteY*hermitePolynomialsDerivative1(nx, x)
             - 2*m_omegaAlpha*y*hermiteX*hermitePolynomialsDerivative1(ny, y)
             + m_omegaAlpha*hermiteX*hermiteY * (m_omegaAlpha*r2 - 2)
             + hermiteY*hermitePolynomialsDerivative2(nx, x)
             + hermiteX*hermitePolynomialsDerivative2(ny, y) );
}

double ManyBodyQuantumDot::singleParticleWFParameters(int nx, int ny, double x, double y) {
    // return derivative of single-particle wave functions w.r.t. alpha

    double r2 = x*x + y*y;
    double hermiteX = hermitePolynomials(nx, x);
    double hermiteY = hermitePolynomials(ny, y);
    return exp(-0.5*m_omegaAlpha*r2) *
           ( - 0.5*m_omega*r2*hermiteX*hermiteY
             + hermitePolynomialsParametersDerivative(nx, x)*hermiteY
             + hermiteX*hermitePolynomialsParametersDerivative(ny, y) );
}

double ManyBodyQuantumDot::hermitePolynomials(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 1;
    }

    else if (energyLevel == 1) {
        return 2*m_omegaAlphaSqrt*position;
    }

    else if (energyLevel == 2) {
        return 4*m_omegaAlpha*position*position - 2;
    }

    else if (energyLevel == 3) {
        return 8*m_omegaAlpha*m_omegaAlphaSqrt*position*position*position - 12*m_omegaAlphaSqrt*position;
    }

    else {
        cout << "Energy level should not exceed n = 3" << endl;
        exit(0);
    }
}

double ManyBodyQuantumDot::hermitePolynomialsDerivative1(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 0;
    }

    else if (energyLevel == 1) {
        return 2*m_omegaAlphaSqrt;
    }

    else if (energyLevel == 2) {
        return 8*m_omegaAlpha*position;
    }

    else if (energyLevel == 3) {
        return 24*m_omegaAlpha*m_omegaAlphaSqrt*position*position - 12*m_omegaAlphaSqrt;
    }

    else {
        cout << "Energy level should not exceed n = 3" << endl;
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
        return 8*m_omegaAlpha;
    }

    else if (energyLevel == 3) {
        return 48*m_omegaAlphaSqrt*m_omegaAlpha*position;
    }

    else {
        cout << "Energy level should not exceed n = 3" << endl;
        exit(0);
    }
}

double ManyBodyQuantumDot::hermitePolynomialsParametersDerivative(int energyLevel, double position) {
    // derivative of Hermite polynomials w.r.t. alpha

    if (energyLevel == 0) {
        return 0;
    }

    else if (energyLevel == 1) {
        return m_omegaSqrt*m_alphaSqrtInv*position;
    }

    else if (energyLevel == 2) {
        return 4*m_omega*position*position;
    }

    else if (energyLevel == 3) {
        return 12*m_omegaAlphaSqrt*m_omega*position*position*position -
               6*m_omegaSqrt*m_alphaSqrtInv*position;
    }

    else {
        cout << "Energy level should not exceed n = 3" << endl;
        exit(0);
    }
}

void ManyBodyQuantumDot::updateSlaterInverse(std::vector<Particle*> particles, int i) {
    // update inverse Slater matrix using Sherman and Morris algorithm

    // get new position of chosen particle
    double xNew = particles[i]->getNewPosition()[0];
    double yNew = particles[i]->getNewPosition()[1];

    // spin-up
    if (i < m_numberOfParticlesHalf) {
        arma::mat slaterUpInverseOld = m_slaterSpinUpInverse;
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            if (j != i) {
                double sum = 0;
                for (int l=0; l < m_numberOfParticlesHalf; l++) {
                    int nx = m_quantumNumbers(l,0);
                    int ny = m_quantumNumbers(l,1);
                    sum += singleParticleWaveFunctions(nx, ny, xNew, yNew) *
                           slaterUpInverseOld(l,j);
                }
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinUpInverse(k,j) = slaterUpInverseOld(k,j) -
                                                 sum * ( slaterUpInverseOld(k,i) / m_ratioSD );
                }
            }
            else {
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinUpInverse(k,i) = slaterUpInverseOld(k,i) / m_ratioSD;
                }
            }
        }
    }
    // spin-down
    else {
        int iDown = i - m_numberOfParticlesHalf;
        arma::mat slaterDownInverseOld = m_slaterSpinDownInverse;
        for (int j=0; j < m_numberOfParticlesHalf; j++) {      
            if (j != iDown) {
                double sum = 0;
                for (int l=0; l < m_numberOfParticlesHalf; l++) {
                    int nx = m_quantumNumbers(l,0);
                    int ny = m_quantumNumbers(l,1);
                    sum += singleParticleWaveFunctions(nx, ny, xNew, yNew) *
                           slaterDownInverseOld(l,j);
                }
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinDownInverse(k,j) = slaterDownInverseOld(k,j) -
                                                   sum * ( slaterDownInverseOld(k,iDown) / m_ratioSD );
                }
            }
            else {
                for (int k=0; k < m_numberOfParticlesHalf; k++) {
                    m_slaterSpinDownInverse(k,iDown) = slaterDownInverseOld(k,iDown) / m_ratioSD;
                }
            }
        }
    }
}

double ManyBodyQuantumDot::evaluate(std::vector<Particle*> particles) {
    // compute the actual determinants to do numerical differentiation

    arma::mat slaterUp = arma::zeros<arma::mat>(m_numberOfParticlesHalf, m_numberOfParticlesHalf);
    arma::mat slaterDown = arma::zeros<arma::mat>(m_numberOfParticlesHalf, m_numberOfParticlesHalf);

    // slater part
    for (int i=0; i < m_numberOfParticlesHalf; i++) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double xUp = m_system->getParticles()[i]->getPosition()[0];
            double yUp = m_system->getParticles()[i]->getPosition()[1];
            double xDown = m_system->getParticles()[i+m_numberOfParticlesHalf]->getPosition()[0];
            double yDown = m_system->getParticles()[i+m_numberOfParticlesHalf]->getPosition()[1];
            slaterUp(i,j) = singleParticleWaveFunctions(nx, ny, xUp, yUp);
            slaterDown(i,j) = singleParticleWaveFunctions(nx, ny, xDown, yDown);
        }
    }

    // correlation part
    double exponent = 0;
    double beta = m_parameters[1];
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=i+1; j < m_numberOfParticles; j++) {
            double r_ij = 0;
            for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                r_ij += (particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim]) *
                        (particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim]);
            }
            r_ij = sqrt(r_ij);
            exponent += ( m_a(i,j)*r_ij ) / ( 1 + beta*r_ij );
        }
    }
    double jastrow = exp(exponent);

    // compute determinants
    double determinantSlaterUp = arma::det(slaterUp);
    double determinantSlaterDown = arma::det(slaterDown);

    return determinantSlaterUp*determinantSlaterDown*jastrow;
}

double ManyBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles) {
    // evaluate full Laplacian of trial wave function divided by trial wave function

    if (m_firstStepLaplacian) {
        // calculate both spin up and down for the first step
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double xUp = particles[i]->getPosition()[0];
                double yUp = particles[i]->getPosition()[1];
                double xDown = particles[i+m_numberOfParticlesHalf]->getPosition()[0];
                double yDown = particles[i+m_numberOfParticlesHalf]->getPosition()[1];
                m_laplacianUp += singleParticleWFLaplacian(nx, ny, xUp, yUp) * m_slaterSpinUpInverse(j,i);
                m_laplacianDown += singleParticleWFLaplacian(nx, ny, xDown, yDown) * m_slaterSpinDownInverse(j,i);
            }
        }
        m_firstStepLaplacian = false;
    }

    // spin-up slater
    if (m_i < m_numberOfParticlesHalf) {
        m_laplacianUp = 0;
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double x = particles[i]->getPosition()[0];
                double y = particles[i]->getPosition()[1];
                m_laplacianUp += singleParticleWFLaplacian(nx, ny, x, y) * m_slaterSpinUpInverse(j,i);
            }
        }
    }

    // spin-down slater
    else {
        m_laplacianDown = 0;
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double x = particles[i+m_numberOfParticlesHalf]->getPosition()[0];
                double y = particles[i+m_numberOfParticlesHalf]->getPosition()[1];
                m_laplacianDown += singleParticleWFLaplacian(nx, ny, x, y) * m_slaterSpinDownInverse(j,i);
            }
        }
    }

    // Laplacian jastrow
    double jastrowLaplacian = 0;
    double beta = m_parameters[1];
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> gradJastrow = gradientJastrow(particles, i);
        jastrowLaplacian += gradJastrow[0]*gradJastrow[0] + gradJastrow[1]*gradJastrow[1];

        double x_i = particles[i]->getPosition()[0];
        double y_i = particles[i]->getPosition()[1];

        for (int k=0; k < i; k++) {
            double x_k = particles[k]->getPosition()[0];
            double y_k = particles[k]->getPosition()[1];
            double r_ki = sqrt( (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i) );

            double factor = 1.0 / (1 + beta*r_ki);
            jastrowLaplacian += ( (m_a(k,i)*factor*factor) / r_ki ) -
                                2*m_a(k,i)*beta*factor*factor*factor;
        }
        for (int k=i+1; k < m_numberOfParticles; k++) {
            double x_k = particles[k]->getPosition()[0];
            double y_k = particles[k]->getPosition()[1];
            double r_ki = sqrt( (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i) );

            double factor = 1.0 / (1 + beta*r_ki);
            jastrowLaplacian += ( (m_a(k,i)*factor*factor) / r_ki ) -
                                2*m_a(k,i)*beta*factor*factor*factor;
        }
    }

    // cross term
    double slaterJastrow = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> gradSlater = gradientSlater(particles, i);
        std::vector<double> gradJastrow = gradientJastrow(particles, i);
        slaterJastrow += gradSlater[0]*gradJastrow[0];
        slaterJastrow += gradSlater[1]*gradJastrow[1];
    }

    return m_laplacianUp + m_laplacianDown + jastrowLaplacian + 2*slaterJastrow;
}

std::vector<double> ManyBodyQuantumDot::gradientSlater(std::vector<Particle *> particles, int i) {
    // return gradient of Slater part of wave function w.r.t. particle i
    // the gradient is returned as a two-dimensional std vector

    std::vector<double> gradient(2);

    // spin-down gradient is zero if i is in spin-up matrix and vice versa
    m_gradientUp[0]   = 0;    m_gradientUp[1]   = 0;
    m_gradientDown[0] = 0;    m_gradientDown[1] = 0;

    // get position of particle i
    double x = particles[i]->getPosition()[0];
    double y = particles[i]->getPosition()[1];

    // spin-up slater
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x, y);
            m_gradientUp[0] += grad[0] * m_slaterSpinUpInverse(j,i);
            m_gradientUp[1] += grad[1] * m_slaterSpinUpInverse(j,i);
        }
    }
    // spin-down slater
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x, y);
            m_gradientDown[0] += grad[0] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
            m_gradientDown[1] += grad[1] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
        }
    }

    // calculate total gradient
    gradient[0] = m_gradientUp[0] + m_gradientDown[0];
    gradient[1] = m_gradientUp[1] + m_gradientDown[1];

    return gradient;
}

std::vector<double> ManyBodyQuantumDot::gradientJastrow(std::vector<Particle *> particles, int i) {
    // return gradient of Jastrow part of wave function w.r.t. particle i
    // the gradient is returned as a two-dimensional std vector

    std::vector<double> ratioJastrow(2);
    ratioJastrow[0] = 0; ratioJastrow[1] = 0;
    double beta = m_parameters[1];

    // get position of particle i
    double x_i = particles[i]->getPosition()[0];
    double y_i = particles[i]->getPosition()[1];

    for (int j=0; j < i; j++) {
        double x_j = particles[j]->getPosition()[0];
        double y_j = particles[j]->getPosition()[1];
        double r_ij = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
        double factor = 1.0 / ( r_ij*(1 + beta*r_ij)*(1 + beta*r_ij) );

        ratioJastrow[0] += (x_i - x_j)*m_a(i,j)*factor;
        ratioJastrow[1] += (y_i - y_j)*m_a(i,j)*factor;
    }
    for (int j=i+1; j < m_numberOfParticles; j++) {
        double x_j = particles[j]->getPosition()[0];
        double y_j = particles[j]->getPosition()[1];
        double r_ij = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
        double factor = 1.0 / ( r_ij*(1 + beta*r_ij)*(1 + beta*r_ij) );

        ratioJastrow[0] += (x_i - x_j)*m_a(i,j)*factor;
        ratioJastrow[1] += (y_i - y_j)*m_a(i,j)*factor;
    }

    return ratioJastrow;
}

std::vector<double> ManyBodyQuantumDot::computeGradient(std::vector<Particle *> particles, int particle) {
    // this function is used in importance sampling to calculate the quantum force
    // returns a two-dimensional std vector that is the gradient of the trial wave function
    // w.r.t. the chosen particle

    std::vector<double> gradient(2);

    // compute Slater and Jastrow gradients for chosen particle
    std::vector<double> gradSlater = gradientSlater(particles, particle);
    std::vector<double> gradJastrow = gradientJastrow(particles, particle);

    // divide by the Slater determinant ratio to compute new quantum force
    if (m_system->getUseOldQuantumForce()) {
        gradient[0] = gradSlater[0] + gradJastrow[0];
        gradient[1] = gradSlater[1] + gradJastrow[1];
    }
    else {
        gradient[0] = gradSlater[0]/m_ratioSD + gradJastrow[0];
        gradient[1] = gradSlater[1]/m_ratioSD + gradJastrow[1];
    }

    return gradient;
}

std::vector<double> ManyBodyQuantumDot::computeParametersGradient(std::vector<Particle *> particles) {
    // calculate gradient of wave function w.r.t. the variational parameters
    // divided by the wave function

    std::vector<double> gradient(2);

    // slater part
    // compute trace of product of inverse Slater matrix and Slater
    // matrix differentiated w.r.t. the variational parameters
    double slaterUp = 0;
    double slaterDown = 0;
    for (int i=0; i < m_numberOfParticlesHalf; i++) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double xUp = particles[i]->getPosition()[0];
            double yUp = particles[i]->getPosition()[1];
            double xDown = particles[i+m_numberOfParticlesHalf]->getPosition()[0];
            double yDown = particles[i+m_numberOfParticlesHalf]->getPosition()[1];
            slaterUp   += singleParticleWFParameters(nx, ny, xUp, yUp)     * m_slaterSpinUpInverse(j,i);
            slaterDown += singleParticleWFParameters(nx, ny, xDown, yDown) * m_slaterSpinDownInverse(j,i);
        }
    }

    // jastrow
    double jastrow = 0;
    double beta = m_parameters[1];
    for (int i=0; i < m_numberOfParticles; i++) {
        double x_i = particles[i]->getPosition()[0];
        double y_i = particles[i]->getPosition()[1];

        for (int j=i+1; j < m_numberOfParticles; j++) {
            double x_j = particles[j]->getPosition()[0];
            double y_j = particles[j]->getPosition()[1];

            double r_ij = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
            double factor = 1.0 / (1 + beta*r_ij);

            jastrow -= m_a(i,j)*r_ij*r_ij*factor*factor;
        }
    }

    // derivative of wave function w.r.t. alpha
    gradient[0] = slaterUp + slaterDown;

    // derivative of wave function w.r.t. beta
    gradient[1] = jastrow;

    return gradient;

}

void ManyBodyQuantumDot::setParameters(std::vector<double> parameters) {
    // update parameter values when doing optimization

    // store new parameters vector
    m_parameters = parameters;

    double alpha = parameters[0];

    // update constants
    m_omegaAlpha            = m_omega*alpha;
    m_omegaAlphaSqrt        = sqrt(m_omegaAlpha);
    m_alphaSqrtInv          = 1.0 / sqrt(alpha);
}
