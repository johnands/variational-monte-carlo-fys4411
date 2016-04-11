#pragma once
#include <vector>
#include "hamiltonian.h"

class HarmonicOscillatorInteracting : public Hamiltonian {
public:
    HarmonicOscillatorInteracting(class System* system, double omega, double a, double gamma, bool useNumerical);
    double computeAnalyticalKineticEnergy(std::vector<class Particle*> particles);
    double computePotentialEnergy(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
    double m_gamma = 0;
    double m_a     = 0;
};
