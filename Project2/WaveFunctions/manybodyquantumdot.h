#pragma once
#include "wavefunction.h"
#include <armadillo>

class ManyBodyQuantumDot : public WaveFunction {
public:
    ManyBodyQuantumDot(class System* system, double alpha, double beta, double omega);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    std::vector<double> computeGradient(std::vector<class Particle*> particles, int particle);
    std::vector<double> computeParametersGradient(std::vector<class Particle*> particles);
    void setUpSlater();
    double singleParticleWaveFunctions(int nx, int ny, double x, double y);
    std::vector<double> singleParticleWFGradient(int nx, int ny, double x, double y);
    double singleParticleWFLaplacian(int nx, int ny, double x, double y);
    double hermitePolynomials(int energyLevel, double position);
    double hermitePolynomialsDerivative1(int energyLevel, double position);
    double hermitePolynomialsDerivative2(int energyLevel, double position);
    void updateRowSlater(std::vector<class Particle*> particles, int i);
    void updateSlaterInverse(std::vector<class Particle*> particles, int i);
    double computeRatio(std::vector<class Particle*> particles, int particle);
    std::vector<double> gradientSlater(std::vector<class Particle*> particles, int i);
    std::vector<double> gradientJastrow(std::vector<class Particle*> particles, int i);

private:
    double m_omega;
    double m_omegaSqrt;
    arma::mat m_a;
    int m_numberOfParticles;
    int m_numberOfParticlesHalf;
    double m_ratio;
    int m_i;
    arma::mat m_quantumNumbers;
    arma::mat m_slaterSpinUp;
    arma::mat m_slaterSpinDown;
    arma::mat m_slaterSpinUpInverse;
    arma::mat m_slaterSpinDownInverse;
    double m_laplacianUp = 0;
    double m_laplacianDown = 0;
    std::vector<double> m_gradientUp = std::vector<double>();
    std::vector<double> m_gradientDown = std::vector<double>();
    bool m_firstStepLaplacian = true;
    bool m_firstStepGradient = true;
};
