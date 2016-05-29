#pragma once
#include <vector>


class Hamiltonian {
public:
    Hamiltonian(class System* system, bool useNumerical, bool useInteraction);
    double computeLocalEnergy(std::vector<class Particle*> particles);
    double computeKineticEnergy(std::vector< class Particle*> particles);
    virtual double computeAnalyticalKineticEnergy(std::vector<class Particle*> particle) = 0;
    virtual double computePotentialEnergy(std::vector<class Particle*> particles) = 0;

    double getKineticEnergy() { return m_kineticEnergy; }
    double getPotentialEnergy() { return m_potentialEnergy; }

protected:
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    bool m_useNumerical = true;
    bool m_useInteraction = true;
    class System* m_system = nullptr;
};

