#pragma once
#include <vector>
#include "hamiltonian.h"

class HarmonicOscillatorQuantumDot2 : public Hamiltonian {
public:
    HarmonicOscillatorQuantumDot2(class System* system, double omega, bool useNumerical, bool useInteraction);
    double computeAnalyticalKineticEnergy(std::vector<class Particle*> particles);
    double computePotentialEnergy(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
};
