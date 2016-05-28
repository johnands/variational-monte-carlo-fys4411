#include "steepestdescent.h"
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "InitialStates/randomuniform.h"
#include "sampler.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <valarray>

using std::cout;
using std::endl;

SteepestDescent::SteepestDescent(System* system, double stepLengthOptimize) {
    m_system = system;
    m_stepLengthOptimize = stepLengthOptimize;
}

void SteepestDescent::optimize(std::vector<double> parameters) {

    int maxNumberOfSteps = 100;
    int stepNumber = 0;
    int numberOfParameters = parameters.size();
    double tolerance = 1e-6;

    m_system->getWaveFunction()->setParameters(parameters);
    // run one time before optimizing
    m_system->runMetropolisSteps((int) 1e5);
    std::vector<double> localEnergyGradient(numberOfParameters);
    //localEnergyGradient[0] = m_system->getSampler()->getAlphaDerivative();
    //localEnergyGradient[1] = m_system->getSampler()->getBetaDerivative();
    for (int j=0; j < numberOfParameters; j++) {
        localEnergyGradient[j] = 2 * ( m_system->getSampler()->getWaveFunctionEnergy()[j] -
                                       m_system->getSampler()->getWaveFunctionDerivative()[j]  *
                                       m_system->getSampler()->getEnergy() );
    }
    double oldAbsoluteGradient = 0;
    for (int j=0; j < numberOfParameters; j++) {
        oldAbsoluteGradient += localEnergyGradient[j]*localEnergyGradient[j];
    }
    cout << "gradient0: " << localEnergyGradient[0] << endl;
    cout << "gradient1: " << localEnergyGradient[1] << endl;


    while (oldAbsoluteGradient > tolerance && stepNumber < maxNumberOfSteps) {

        cout << "****** Iteration " << stepNumber << " ******" << endl;

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set parameter values
        m_system->getWaveFunction()->setParameters(parameters);

        // run metropolis steps
        m_system->runMetropolisSteps((int) 1e5);

        // compute gradient of exp. value of local energy w.r.t. the variational parameters
        std::vector<double> localEnergyGradient(numberOfParameters);
        for (int j=0; j < numberOfParameters; j++) {
            localEnergyGradient[j] = 2 * ( m_system->getSampler()->getWaveFunctionEnergy()[j] -
                                           m_system->getSampler()->getWaveFunctionDerivative()[j]  *
                                           m_system->getSampler()->getEnergy() );
        }
        //localEnergyGradient[0] = m_system->getSampler()->getAlphaDerivative();
        //localEnergyGradient[1] = m_system->getSampler()->getBetaDerivative();
        cout << "gradient0: " << localEnergyGradient[0] << endl;
        cout << "gradient1: " << localEnergyGradient[1] << endl;

        // calculate absolute value of gradient
        double newAbsoluteGradient = 0;
        for (int j=0; j < numberOfParameters; j++) {
            newAbsoluteGradient += localEnergyGradient[j]*localEnergyGradient[j];
        }
        //cout << "Gradient: " << newAbsoluteGradient << endl;

        // update alpha and beta if new gradient is less than old
        // if not, reduce step length
        if (newAbsoluteGradient < oldAbsoluteGradient) {
            for (int j=0; j < numberOfParameters; j++) {
                parameters[j] -= m_stepLengthOptimize*localEnergyGradient[j];
            }
        }
        else {
            m_stepLengthOptimize *= 0.8;
            cout << "New step length: " << m_stepLengthOptimize << endl;
        }

        for (int j=0; j < numberOfParameters; j++) {
            cout << " Parameter " << j+1 << " : " << parameters.at(j) << endl;
        }
        cout << endl;

        // before new iteration
        oldAbsoluteGradient = newAbsoluteGradient;
        stepNumber++;

        /*if (std::all_of(localEnergyGradient.begin(), localEnergyGradient.end(), [&tolerance](double j)
        { return std::fabs(j) < tolerance; } ))
            //cout << "Number of steps run: " << i << endl;
            break;*/
    }

    for (int j=0; j < numberOfParameters; j++) {
        cout << " Optimal parameter " << j+1 << " : " << parameters.at(j) << endl;
    }
    cout << endl;

    // run many Metropolis steps with the optimal parameters

    // make initial state
    m_system->getInitialState()->setupInitialState();

    // set value of parameters
    m_system->getWaveFunction()->setParameters(parameters);

    // run metropolis steps
    m_system->runMetropolisSteps((int) 1e6);
}
