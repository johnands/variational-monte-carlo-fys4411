#pragma once
#include "wavefunction.h"

class QuantumDotTwoElectrons : public WaveFunction {
public:
    QuantumDotTwoElectrons(class System* system, double alpha, double beta, double omega, double a);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    std::vector<double> computeGradient(std::vector<class Particle*> particles, int particle);
    std::vector<double> computeParametersGradient(std::vector<class Particle *> particles);

private:
    double m_omega;
    double m_a;
};


