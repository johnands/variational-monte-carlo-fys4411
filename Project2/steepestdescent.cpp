#include "steepestdescent.h"
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "InitialStates/randomuniform.h"
#include "sampler.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

SteepestDescent::SteepestDescent(System* system, double stepLengthOptimize,
                                 int numberOfMetropolisSteps, int maxNumberOfIterations) {
    m_system = system;
    m_stepLengthOptimize = stepLengthOptimize;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_maxNumberOfIterations = maxNumberOfIterations;
}

void SteepestDescent::optimize(std::vector<double> parameters) {

    int stepNumber = 0;
    int numberOfParameters = parameters.size();
    double tolerance = 1e-5;

    // run one time before optimizing
    m_system->getWaveFunction()->setParameters(parameters);
    m_system->runMetropolisSteps((int) 1e5);

    std::vector<double> localEnergyGradient(numberOfParameters);
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


    while (m_stepLengthOptimize > tolerance && stepNumber < m_maxNumberOfIterations) {

        cout << "****** Iteration " << stepNumber << " ******" << endl;

        // make initial state
        m_system->getInitialState()->setupInitialState();

        // set parameter values
        m_system->getWaveFunction()->setParameters(parameters);

        // run metropolis steps
        m_system->runMetropolisSteps((int) m_numberOfMetropolisSteps);

        // compute gradient of exp. value of local energy w.r.t. the variational parameters
        std::vector<double> localEnergyGradient(numberOfParameters);
        for (int j=0; j < numberOfParameters; j++) {
            localEnergyGradient[j] = 2 * ( m_system->getSampler()->getWaveFunctionEnergy()[j] -
                                           m_system->getSampler()->getWaveFunctionDerivative()[j]  *
                                           m_system->getSampler()->getEnergy() );
        }

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
    }

    // write the optimal parameters to screen
    for (int j=0; j < numberOfParameters; j++) {
        cout << " Optimal parameter " << j+1 << " : " << parameters.at(j) << endl;
    }
    cout << endl;
}


/*
  -- System info --
 Optimal parameter 1 : 0.7961730521
 Optimal parameter 2 : 0.8189088609

 Time elapsed: 41.783 s
 step number 899999


 Number of particles  : 20
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 2
 Parameter 1 : 0.79617
 Parameter 2 : 0.81891

  -- Results --
 Energy : 156.0638026
 Variance : 1.30429743e-06
 Acceptance rate : 0.531308
 Number of accepted steps 531308

Press <RETURN> to close this window...
*/
