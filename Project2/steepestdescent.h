#pragma once
#include <vector>

class SteepestDescent
{
public:
    SteepestDescent(class System* system, double stepLength,
                    int numberOfMetropolisSteps, int maxNumberOfIterations);
    void optimize(std::vector<double> parameters);

private:
    class System* m_system = nullptr;
    double m_stepLengthOptimize = 0;
    int m_numberOfMetropolisSteps;
    int m_maxNumberOfIterations;
};

