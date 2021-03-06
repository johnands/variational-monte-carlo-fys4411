#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    //std::vector<double> computeGradient(std::vector<class Particle*> particles);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    void setAlpha( double alpha );
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeGradient(std::vector<class Particle*> particles) = 0;
    virtual double computeAlphaDerivative(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

