#pragma once

class SteepestDescent
{
public:
    SteepestDescent(class System* system, double stepLength);
    void optimize(double initialAlpha);

private:
    class System* m_system = nullptr;
    double m_stepLengthOptimize = 0;
};

