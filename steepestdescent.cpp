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
    double tolerance = 1e-6;
    double alpha = initialAlpha;
    double oldEnergy = 1e10;
    for (int i=0; i < maxNumberOfSteps; i++) {

        if (i > 0) {
            oldEnergy = m_system->getSampler()->getEnergy();
        }

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set value of alpha
        m_system->getWaveFunction()->setAlpha(alpha);

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
        else {
            // compute new alpha
            alpha -= m_stepLengthOptimize*localEnergyDerivative;
        }

        cout << "newAlhpa = " << alpha << endl;

        if ( localEnergyDerivative < tolerance ) {
            break;
        }

    }
    cout << "Optimal alpha = " << alpha << endl;

    // run many Metropolis steps with the optimal alpha

    // make initial state
    m_system->getInitialState()->setupInitialState();

    // set value of alpha
    m_system->getWaveFunction()->setAlpha(alpha);

    // run metropolis steps
    m_system->runMetropolisSteps((int) 1e5, false, false, false);
}
