#include "manybodyquantumdot.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <armadillo>

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

    double ratioSD = 0;

    // get new position of particle i
    double xNew = particles[i]->getNewPosition()[0];
    double yNew = particles[i]->getNewPosition()[1];

    // spin-up slater
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratio += singleParticleWaveFunctions(nx, ny, xNew, yNew)*m_slaterSpinUpInverse(j,i);
        }
    }
    // spin-down slater
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratio += singleParticleWaveFunctions(nx, ny, xNew, yNew) *
                    m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
        }
    }

    // compute jastrow factor
    double exponent = 0;
    double beta = m_parameters[1];
    for (int j=0; j < m_numberOfParticles; j++) {
        if (!(j == i)) {
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

    return ratioSD*ratioJastrow;
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

double ManyBodyQuantumDot::hermitePolynomials(int energyLevel, double position) {

    if (energyLevel == 0) {
        return 1;
    }

    else if (energyLevel == 1) {
        return 2*m_omegaSqrt*position;
    }

    else if (energyLevel == 2) {
        return 4*m_omegaSqrt*position*position;
    }

    else if (energyLevel == 3) {
        return 8*m_omegaSqrt*position*position*position - 12*m_omegaSqrt*position;
    }

    else if (energyLevel == 4) {
        return 16*m_omegaSqrt*position*position*position*position -
               48*m_omegaSqrt*position*position + 12;
    }

    else {
        cout << "The energy level is too high" << endl;
        exit(0);
    }
}

void ManyBodyQuantumDot::updateRowSlater(int i) {
    // update row corresponding to particle i in Slater matrix

    // maybe update Slater and inverse Slater here??
    // need both current Slater and new slater to compute new inverse Slater
    // or make matrices NewSlaterUp and NewSlaterDown that are
    // temporary so that I have access to both new and old Slater matrices
    // then I must set SlaterUp = NewSlaterUp before new cycle

    // get new position of particle i
    double xNew = m_system->getParticles()[i]->getNewPosition()[0];
    double yNew = m_system->getParticles()[i]->getNewPosition()[1];

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
            m_slaterSpinDown(i,j) = singleParticleWaveFunctions(nx, ny, xNew, yNew);
        }
    }
    updateRowSlaterInverse(i);
}

void ManyBodyQuantumDot::updateRowSlaterInverse(int i) {
    // update row corresponding to particle i in inverse Slater matrix

    // get new position of particle i
    double xNew = m_system->getParticles()[i]->getNewPosition()[0];
    double yNew = m_system->getParticles()[i]->getNewPosition()[1];

    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            if (j == i) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                m_slaterSpinUp(i,j) = singleParticleWaveFunctions(nx, ny, xNew, yNew);
            }
            else {

            }
        }
    }
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            m_slaterSpinDown(i,j) = singleParticleWaveFunctions(nx, ny, xNew, yNew);
        }
    }

}

double ManyBodyQuantumDot::evaluate(std::vector<Particle*> particles) {

}

double ManyBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles) {
    // Calculate laplacian of trial wavefunction divided by trial wavefunction

}

std::vector<double> ManyBodyQuantumDot::computeGradient(std::vector<Particle *> particles) {
    // calculate gradient of trial wavefunction divided by trial wavefunction
    // used to calculate drift velocity / quantum force


}

std::vector<double> ManyBodyQuantumDot::computeParametersGradient(std::vector<Particle *> particles) {
    // calculate gradient of wave function w.r.t. the variational parameters
    // divided by the wave function


}
