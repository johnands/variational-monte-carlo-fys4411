#pragma once
#include "wavefunction.h"

class InteractingGaussian : public WaveFunction {
public:
    InteractingGaussian(class System* system, double alpha, double beta, double a);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    std::vector<double> computeGradient(std::vector<class Particle*> particles);
    double computeAlphaDerivative(std::vector<class Particle *> particles);

private:
    double m_a = 0;
};
