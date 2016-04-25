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
    std::vector<double> singleParticleWFGradient(int nx, int ny, double x, double y);
    double singleParticleWFLaplacian(int nx, int ny, double x, double y);
    double hermitePolynomials(int energyLevel, double position);
    double hermitePolynomialsDerivative1(int energyLevel, double position);
    double hermitepolynomialsDerivative2(int energyLevel, double position);
    void updateRowSlater(int i);
    void updateRowSlaterInverse(int i);
    double computeRatio(std::vector<class Particle*> particles, int particle);

private:
    double m_omega;
    double m_omegaSqrt;
    double m_a;
    int m_numberOfParticles;
    int m_numberOfParticlesHalf;
    double m_ratio;
    arma::mat m_quantumNumbers;
    arma::mat m_slaterSpinUpOld;
    arma::mat m_slaterSpinDownOld;
    arma::mat m_slaterSpinUpNew;
    arma::mat m_slaterSpinDownNew;
    arma::mat m_slaterSpinUpInverse;
    arma::mat m_slaterSpinDownInverse;
};
