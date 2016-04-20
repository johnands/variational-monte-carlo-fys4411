#include "steepestdescent.h"
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "InitialStates/randomuniform.h"
#include "sampler.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using std::cout;
using std::endl;

SteepestDescent::SteepestDescent(System* system, double stepLengthOptimize) {
    m_system = system;
    m_stepLengthOptimize = stepLengthOptimize;
}

void SteepestDescent::optimize(std::vector<double> parameters) {

    int maxNumberOfSteps = 30;
    int numberOfParameters = parameters.size();
    double tolerance = 1e-6;
    double oldEnergy = 1e10;
    for (int i=0; i < maxNumberOfSteps; i++) {

        if (i > 0) {
            oldEnergy = m_system->getSampler()->getEnergy();
        }

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set parameter values
        m_system->getWaveFunction()->setParameters(parameters);

        // run metropolis steps
        m_system->runMetropolisSteps((int) 1e5, false, false, false);

        double newEnergy = m_system->getSampler()->getEnergy();
        cout << "New Energy: " << newEnergy << endl;

        // compute gradient of exp. value of local energy w.r.t. the variational parameters
        std::vector<double> localEnergyGradient(numberOfParameters);
        for (int i=0; i < numberOfParameters; i++) {
            localEnergyGradient[i] = 2 * ( m_system->getSampler()->getWaveFunctionEnergy()[i] -
                                           m_system->getSampler()->getWaveFunctionDerivative()[i]  *
                                           newEnergy );
        }

        if (newEnergy > oldEnergy) {
            m_stepLengthOptimize /= 2.0;
            cout << "New step length: " << m_stepLengthOptimize << endl;
        }
        else {
            // compute new paramters
            for (int i=0; i < numberOfParameters; i++) {
                parameters[i] -= m_stepLengthOptimize*localEnergyGradient[i];
            }
        }

        for (int i=0; i < numberOfParameters; i++) {
            cout << " Parameter " << i+1 << " : " << parameters.at(i) << endl;
        }

        if (std::all_of(parameters.begin(), parameters.end(), [&tolerance](double i){ return i < tolerance; } ))
            break;
    }
    for (int i=0; i < numberOfParameters; i++) {
        cout << " Optimal parameter " << i+1 << " : " << parameters.at(i) << endl;
    }

    // run many Metropolis steps with the optimal parameters

    // make initial state
    m_system->getInitialState()->setupInitialState();

    // set value of parameters
    m_system->getWaveFunction()->setParameters(parameters);

    // run metropolis steps
    m_system->runMetropolisSteps((int) 1e6, false, false, false);
}
