#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, bool useNumerical, bool useInteraction);
    double computePotentialEnergy(std::vector<class Particle*> particles);
    double computeAnalyticalKineticEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};

