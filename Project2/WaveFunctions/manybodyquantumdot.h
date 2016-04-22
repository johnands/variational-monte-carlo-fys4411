#pragma once
#include "wavefunction.h"
#include <armadillo>

class ManyBodyQuantumDot : public WaveFunction {
public:
    ManyBodyQuantumDot(class System* system, double alpha, double beta, double omega, double a);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    std::vector<double> computeGradient(std::vector<class Particle*> particles);
    std::vector<double> computeParametersGradient(std::vector<class Particle*> particles);
    void setUpSlater();
    double singleParticleWaveFunctions(int nx, int ny, double x, double y);
    double hermitePolynomials(int energyLevel, double position);
    void updateRowSlater(int i);
    void updateRowSlaterInverse(int i);
    double computeRatio(std::vector<class Particle*> particles, int particle);

private:
    double m_omega;
    double m_omegaSqrt;
    double m_a;
    int m_numberOfParticles;
    int m_numberOfParticlesHalf;
    arma::mat m_quantumNumbers;
    arma::mat m_slaterSpinUp;
    arma::mat m_slaterSpinDown;
    arma::mat m_slaterSpinUpInverse;
    arma::mat m_slaterSpinDownInverse;
};
