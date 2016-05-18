#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    //std::vector<double> computeGradient(std::vector<class Particle*> particles);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual void setParameters(std::vector<double> parameters) = 0;
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeGradient(std::vector<class Particle*> particles, int particle) = 0;
    virtual std::vector<double> computeParametersGradient(std::vector<class Particle*> particles) = 0;
    virtual double computeRatio(std::vector<class Particle*> particles, int particle) {}
    virtual void updateSlaterInverse(std::vector<class Particle*> particles, int i) {}

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

