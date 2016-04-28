#pragma once
#include <vector>
#include "hamiltonian.h"

class HOManyBodyQuantumDot : public Hamiltonian {
public:
    HOManyBodyQuantumDot(class System* system, double omega, bool useNumerical);
    double computeAnalyticalKineticEnergy(std::vector<class Particle*> particles);
    double computePotentialEnergy(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
};
