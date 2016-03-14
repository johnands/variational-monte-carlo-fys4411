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
    double tolerance = 1e-10;
    double oldAlpha = initialAlpha;
    double oldEnergy = 1e10;
    for (int i=0; i < maxNumberOfSteps; i++) {

        if (i > 0) {
            oldEnergy = m_system->getSampler()->getEnergy();
        }

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set value of alpha
        m_system->getWaveFunction()->setAlpha(oldAlpha);

        // run metropolis steps
        m_system->runMetropolisSteps((int) 1e5, false, false, false);

        double newEnergy = m_system->getSampler()->getEnergy();

        // compute derivative of exp. value of local energy w.r.t. alpha
        double localEnergyDerivative = 2 * ( m_system->getSampler()->getWaveFunctionEnergy() -
                                             m_system->getSampler()->getWaveFunctionDerivative() *
                                             newEnergy );

        if (newEnergy > oldEnergy) {
            m_stepLengthOptimize /= 2.0;
            cout << "New step length: " << m_stepLengthOptimize << endl;
        }

        // compute new alpha
        double newAlpha = oldAlpha - m_stepLengthOptimize*localEnergyDerivative;
        cout << "newAlhpa = " << newAlpha << endl;
        cout << "oldAlpha = " << oldAlpha << endl;
        /*if ( std::abs(newAlpha - oldAlpha) < tolerance ) {
            break;
        }*/
        if ( localEnergyDerivative < tolerance ) {
            break;
        }

        // before new iteration
        oldAlpha = newAlpha;
    }
    cout << "Optimal alpha = " << oldAlpha << endl;

    // run many Metropolis steps with the optimal alpha

    // make initial state
    m_system->getInitialState()->setupInitialState();

    // set value of alpha
    m_system->getWaveFunction()->setAlpha(oldAlpha);

    // run metropolis steps
    m_system->runMetropolisSteps((int) 1e6, false, false, false);
}
