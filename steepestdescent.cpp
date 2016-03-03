#include "steepestdescent.h"
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "InitialStates/randomuniform.h"
#include "sampler.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

SteepestDescent::SteepestDescent(System* system, double stepLengthOptimize)
{
    m_system = system;
    m_stepLengthOptimize = stepLengthOptimize;
}

void SteepestDescent::optimize(double initialAlpha) {

    int maxNumberOfSteps = 30;
    double tolerance = 0.01;
    double oldAlpha = initialAlpha;
    for (int i=0; i < maxNumberOfSteps; i++) {

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set value of alpha
        m_system->getWaveFunction()->setAlpha(oldAlpha);

        // run metropolis steps
        m_system->runMetropolisSteps((int) 1e5, false, false);

        // compute derivative of exp. value of local energy w.r.t. alpha
        double localEnergyDerivative = 2 * ( m_system->getSampler()->getWaveFunctionEnergy() -
                                             m_system->getSampler()->getWaveFunctionDerivative() *
                                             m_system->getSampler()->getEnergy() );

        cout << "localEnergyDerivative = " << localEnergyDerivative << endl;

        // compute new alpha
        double newAlpha = oldAlpha - m_stepLengthOptimize*localEnergyDerivative;
        cout << "newAlhpa = " << newAlpha << endl;
        cout << "oldAlpha = " << oldAlpha << endl;
        cout << std::abs(newAlpha-oldAlpha) << endl;
        if ( std::abs(newAlpha - oldAlpha) < tolerance ) {
            break;
        }
        // before new iteration
        oldAlpha = newAlpha;
    }
    cout << "Optimal alpha = " << oldAlpha << endl;
}
