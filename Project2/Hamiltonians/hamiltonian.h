#pragma once
#include <vector>


class Hamiltonian {
public:
    Hamiltonian(class System* system, bool useNumerical, bool useInteraction);
    double computeLocalEnergy(std::vector<class Particle*> particles);
    double computeKineticEnergy(std::vector< class Particle*> particles);
    virtual double computeAnalyticalKineticEnergy(std::vector<class Particle*> particle) = 0;
    virtual double computePotentialEnergy(std::vector<class Particle*> particles) = 0;

protected:
    bool m_useNumerical = true;
    bool m_useInteraction = true;
    class System* m_system = nullptr;
};

