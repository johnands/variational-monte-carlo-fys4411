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

    double ratioSD = 0;

    // get new position of particle i
    double xNew = particles[i]->getNewPosition()[0];
    double yNew = particles[i]->getNewPosition()[1];

    //cout << "xNew " << xNew << endl;
    //cout << "yNew " << yNew << endl;

    // spin-up slater
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratioSD += singleParticleWaveFunctions(nx, ny, xNew, yNew)*m_slaterSpinUpInverse(j,i);
            //cout << m_slaterSpinUpInverse(j,i) << endl;
        }
    }
    // spin-down slater
    else {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            ratioSD += singleParticleWaveFunctions(nx, ny, xNew, yNew) *
                       m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
            //cout << m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf) << endl;
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

    /*for (int j=0; j < m_numberOfParticles; j++) {
        if (j != i) {
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
    }*/
    double ratioJastrow = exp(exponent1 - exponent2);
    //cout << "Slater " << ratioSD << endl;
    //cout << "Jastrow " << ratioJastrow << endl;

    m_ratioSD = ratioSD;
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

    cout << m_a << endl;

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

    m_slaterSpinUpInverse = m_slaterSpinUp.i();
    m_slaterSpinDownInverse = m_slaterSpinDown.i();
    //cout << m_slaterSpinUpInverse(0,0) << endl;

    cout << "Initial Slater Up: " << m_slaterSpinUp << endl;
    cout << "Initial Slater Down: " << m_slaterSpinDown << endl;

    cout << "Initial Slater Up Inverse: " << m_slaterSpinUpInverse << endl;
    cout << "Initial Slater Down Inverse: " << m_slaterSpinDownInverse << endl;
}

double ManyBodyQuantumDot::singleParticleWaveFunctions(int nx, int ny, double x, double y) {
    // evaluate single particle wave functions

    return hermitePolynomials(nx, x)*hermitePolynomials(ny, y) *
           exp(-0.5*m_omegaAlpha*(x*x + y*y));
}

std::vector<double> ManyBodyQuantumDot::singleParticleWFGradient(int nx, int ny, double x, double y) {

    std::vector<double> gradient(2);
    double r2 = x*x + y*y;
    gradient[0] =  exp(-0.5*m_omegaAlpha*r2) * hermitePolynomials(ny, y) *
                   ( hermitePolynomialsDerivative1(nx, x) - hermitePolynomials(nx, x)*m_omegaAlpha*x);
    gradient[1] =  exp(-0.5*m_omegaAlpha*r2) * hermitePolynomials(nx, x) *
                   ( hermitePolynomialsDerivative1(ny, y) - hermitePolynomials(ny, y)*m_omegaAlpha*y);

    return gradient;
}

double ManyBodyQuantumDot::singleParticleWFLaplacian(int nx, int ny, double x, double y) {

    double r2 = x*x + y*y;
    return exp(-0.5*m_omegaAlpha*r2) *
           ( - 2*m_omegaAlpha*x*hermitePolynomials(ny, y)*hermitePolynomialsDerivative1(nx, x)
             - 2*m_omegaAlpha*y*hermitePolynomials(nx, x)*hermitePolynomialsDerivative1(ny, y)
             + m_omegaAlpha*hermitePolynomials(nx, x)*hermitePolynomials(ny, y) * (m_omegaAlpha*r2 - 2)
             + hermitePolynomials(ny, y)*hermitePolynomialsDerivative2(nx, x)
             + hermitePolynomials(nx, x)*hermitePolynomialsDerivative2(ny, y) );
}

double ManyBodyQuantumDot::singleParticleWFParameters(int nx, int ny, double x, double y) {

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
        //return 2;
    }

    else if (energyLevel == 2) {
        return 8*m_omegaAlpha*position;
        //return 8*m_omegaSqrt*position;
    }

    else if (energyLevel == 3) {
        return 24*m_omegaAlpha*m_omegaAlphaSqrt*position*position - 12*m_omegaAlphaSqrt;
        //return 24*m_omega*position*position - 12;
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
        //return 8;
    }

    else if (energyLevel == 3) {
        return 48*m_omegaAlphaSqrt*m_omegaAlpha*position;
        //return 48*m_omegaSqrt*position;
    }

    else {
        cout << "Energy level should not exceed n = 3" << endl;
        exit(0);
    }
}

double ManyBodyQuantumDot::hermitePolynomialsParametersDerivative(int energyLevel, double position) {

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
    // update row corresponding to particle i in inverse Slater matrix

    // get new position of particle i
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
            // need to adjust i. for N = 6, i=4 is equal to i-3=1 in spin down matrix         
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

    //cout << "Inverse: " << m_slaterSpinUpInverse << endl;
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

    double determinantSlaterUp = arma::det(slaterUp);
    double determinantSlaterDown = arma::det(slaterDown);

    //cout << "Evaluate: " << slaterUp.i() << endl;

    return determinantSlaterUp*determinantSlaterDown*jastrow;
}

double ManyBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles) {
    // must calculate laplacian for both spin up and spin down first step
    // after that I only need to calculate one of them, they must therefore be stored
    // this function calculates total laplacian, i.e. w.r.t. all particles

    if (m_firstStepLaplacian) {
        // calculate both spin up and down if first step
        for (int i=0; i < m_numberOfParticlesHalf; i++) {
            for (int j=0; j < m_numberOfParticlesHalf; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                double xUp = particles[i]->getPosition()[0];
                double yUp = particles[i]->getPosition()[1];
                double xDown = particles[i+m_numberOfParticlesHalf]->getPosition()[0];
                double yDown = particles[i+m_numberOfParticlesHalf]->getPosition()[1];
                // singleParticleWFLaplacian returns the full Laplacian of the specific
                // single-particle wave function
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
                // singleParticleWFLaplacian returns the full Laplacian of the specific
                // single-particle wave function
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

        /*for (int k=0; k < m_numberOfParticles; k++) {

            if (k != i) {
                double x_k = particles[k]->getPosition()[0];
                double y_k = particles[k]->getPosition()[1];
                double r_ki = sqrt( (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i) );

                double factor = 1.0 / (1 + beta*r_ki);
                jastrowLaplacian += ( (m_a(k,i)*factor*factor) / r_ki ) -
                                    2*m_a(k,i)*beta*factor*factor*factor;
            }
        }*/
    }

    // cross term
    double slaterJastrow = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> gradSlater = gradientSlater(particles, i);
        std::vector<double> gradJastrow = gradientJastrow(particles, i);
        slaterJastrow += gradSlater[0]*gradJastrow[0];
        slaterJastrow += gradSlater[1]*gradJastrow[1];
    }
    //cout << m_laplacianUp + m_laplacianDown + jastrowLaplacian + 2*crossTerm << endl;
    //cout << m_laplacianUp << endl;
    //cout << m_laplacianDown << endl;
    //cout << m_laplacianUp + m_laplacianDown << endl;
    //cout << slaterJastrow << endl;
    //cout << jastrowLaplacian << endl;


    return m_laplacianUp + m_laplacianDown + jastrowLaplacian + 2*slaterJastrow;
}

std::vector<double> ManyBodyQuantumDot::gradientSlater(std::vector<Particle *> particles, int i) {
    // calculates gradient of spin up/spin down slater w.r.t. one particle only,
    // but for both dimensions
    // i is the particle index that we calculate gradient for

    std::vector<double> gradient(2);

    /*if (m_firstStepGradient) {
        // need to calculate both spin up and down first step
        // make vectors and fill with zeros
        m_gradientUp.resize(2); m_gradientDown.resize(2);
        m_gradientUp[0] = 0; m_gradientDown[1] = 0;
        m_gradientDown[0] = 0; m_gradientDown[1] = 0;

        // calculate gradients for both spin up and spin down
        // gradient w.r.t. each particle is a dot product between the gradient of
        // the single-particle wave functions evaluated for particle i and the
        // i-th column of the inverse Slater matrix
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double x = particles[i]->getPosition()[0];
            double y = particles[i]->getPosition()[1];
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x, y);

            if (i < m_numberOfParticlesHalf) {
                m_gradientUp[0] += grad[0] * m_slaterSpinUpInverse(j,i);
                m_gradientUp[1] += grad[1] * m_slaterSpinUpInverse(j,i);
            }
            else {
                m_gradientDown[0] += grad[0] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
                m_gradientDown[1] += grad[1] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
            }
        }
        m_firstStepGradient = false;
    }*/

    m_gradientUp[0] = 0; m_gradientUp[1] = 0;
    m_gradientDown[0] = 0; m_gradientDown[1] = 0;

    // spin-up slater
    if (i < m_numberOfParticlesHalf) {
        for (int j=0; j < m_numberOfParticlesHalf; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double x = particles[i]->getPosition()[0];
            double y = particles[i]->getPosition()[1];
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
            double x = particles[i]->getPosition()[0];
            double y = particles[i]->getPosition()[1];
            std::vector<double> grad = singleParticleWFGradient(nx, ny, x, y);
            m_gradientDown[0] += grad[0] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
            m_gradientDown[1] += grad[1] * m_slaterSpinDownInverse(j,i-m_numberOfParticlesHalf);
        }
    }

    gradient[0] = m_gradientUp[0] + m_gradientDown[0];
    gradient[1] = m_gradientUp[1] + m_gradientDown[1];

    return gradient;
}

std::vector<double> ManyBodyQuantumDot::gradientJastrow(std::vector<Particle *> particles, int i) {
    // compute graident jastrow ratio w.r.t. to particle number i

    // compute jastrow factor
    std::vector<double> ratioJastrow(2);
    ratioJastrow[0] = 0; ratioJastrow[1] = 0;
    double beta = m_parameters[1];
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

    /*for (int j=0; j < m_numberOfParticles; j++) {
        if (j != i) {
            double x_j = particles[j]->getPosition()[0];
            double y_j = particles[j]->getPosition()[1];
            double r_ij = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
            double factor = 1.0 / ( r_ij*(1 + beta*r_ij)*(1 + beta*r_ij) );

            ratioJastrow[0] += (x_i - x_j)*m_a(i,j)*factor;
            ratioJastrow[1] += (y_i - y_j)*m_a(i,j)*factor;
        }
    }*/

    return ratioJastrow;
}

std::vector<double> ManyBodyQuantumDot::computeGradient(std::vector<Particle *> particles, int particle) {
    // this function is used in importance sampling to calculate the drift force
    // returns a two-dimensional vector that is the gradient of the trial wave function w.r.t to particle i
    // divided by the trial wave function according to formula 16.10

    /*std::vector<double> gradient(2*m_numberOfParticles);

    std::vector<double> gradSlater(2*m_numberOfParticles);
    std::vector<double> gradJastrow(2*m_numberOfParticles);
    for (int j=0; j < m_numberOfParticles; j++) {
        std::vector<double> grad1 = gradientSlater(particles, j);
        std::vector<double> grad2 = gradientJastrow(particles, j);
        gradSlater[2*j]    = grad1[0];
        gradSlater[2*j+1]  = grad1[1];
        gradJastrow[2*j]   = grad2[0];
        gradJastrow[2*j+1] = grad2[1];
    }

    if (m_quantumForceOld) {
        for (int j=0; j < 2*m_numberOfParticles; j++) {
        gradient[j] = gradSlater[j] + gradJastrow[j];
        }
        m_quantumForceOld = false;
    }
    else {
    for (int j=0; j < 2*m_numberOfParticles; j++) {
        gradient[j] = gradSlater[j]/m_ratioSD + gradJastrow[j];
        m_quantumForceOld = true;
    }

    }*/


    std::vector<double> gradient(2);

    /*std::vector<double> gradSlater(2);
    std::vector<double> gradJastrow(2);
    gradSlater[0] = 0; gradSlater[1] = 0;
    gradJastrow[0] = 0; gradJastrow[1] = 0;
    for (int j=0; j < m_numberOfParticles; j++) {
        std::vector<double> grad1 = gradientSlater(particles, j);
        std::vector<double> grad2 = gradientJastrow(particles, j);
        gradSlater[0] += grad1[0];
        gradSlater[1] += grad1[1];
        gradJastrow[0] += grad2[0];
        gradJastrow[1] += grad2[1];
    }



    //cout << m_ratioSD << endl;
    //cout << m_ratioSD << endl;
    if (m_quantumForceOld) {
        gradient[0] = (gradSlater[0] + gradJastrow[0]);
        gradient[1] = (gradSlater[1] + gradJastrow[1]);
        m_quantumForceOld = false;
    }
    else {
        gradient[0] = (gradSlater[0]/m_ratioSD + gradJastrow[0]);
        gradient[1] = (gradSlater[1]/m_ratioSD + gradJastrow[1]);
        m_quantumForceOld = true;
    }*/


    std::vector<double> gradSlater = gradientSlater(particles, particle);
    std::vector<double> gradJastrow = gradientJastrow(particles, particle);
    if (m_quantumForceOld) {
        gradient[0] = gradSlater[0] + gradJastrow[0];
        gradient[1] = gradSlater[1] + gradJastrow[1];
        m_quantumForceOld = false;
    }
    else {
        gradient[0] = gradSlater[0]/m_ratioSD + gradJastrow[0];
        gradient[1] = gradSlater[1]/m_ratioSD + gradJastrow[1];
        m_quantumForceOld = true;
    }
    /*cout << gradient[0] << endl;
    cout << gradient[1] << endl;
    cout << gradient[2] << endl;
    cout << gradient[3] << endl;*/

    return gradient;
}



std::vector<double> ManyBodyQuantumDot::computeParametersGradient(std::vector<Particle *> particles) {
    // calculate gradient of wave function w.r.t. the variational parameters
    // divided by the wave function

    std::vector<double> gradient(2);

    // slater part
    // compute trace of product of inverse Slater matrix and Slater matrix differentiated w.r.t.
    // the variational parameters
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
            // singleParticleWFLaplacian returns the full Laplacian of the specific
            // single-particle wave function
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

    gradient[0] = slaterUp + slaterDown;
    gradient[1] = jastrow;

    return gradient;

}

void ManyBodyQuantumDot::setParameters(std::vector<double> parameters) {
    // overloading setParameters in WaveFunction

    m_parameters = parameters;
    double alpha = parameters[0];
    m_omegaAlpha            = m_omega*alpha;
    m_omegaAlphaSqrt        = sqrt(m_omegaAlpha);
    m_alphaSqrtInv          = 1.0/sqrt(alpha);
}
