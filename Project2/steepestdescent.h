#pragma once
#include <vector>

class SteepestDescent
{
public:
    SteepestDescent(class System* system, double stepLength);
    void optimize(std::vector<double> parameters);

private:
    class System* m_system = nullptr;
    double m_stepLengthOptimize = 0;
};

