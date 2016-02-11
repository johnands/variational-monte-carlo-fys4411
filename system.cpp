#include "system.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <time.h>

using std::cout;
using std::endl;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle
    int dimension = Random::nextInt(m_numberOfDimensions);  // choose random dimension
    double change = Random::nextDouble()*m_stepLength*2-1;  // propose change of particle's position

    // get old wavefunction
    double waveFuncOld = m_waveFunction->evaluate(m_particles);

    // adjust position
    m_particles[particle]->adjustPosition(change, dimension);

    // get new wavefunction
    double waveFuncNew = m_waveFunction->evaluate(m_particles);

    // accept/reject new position using Metropolis algorithm
    double ratio = pow(waveFuncNew, 2) / pow(waveFuncOld, 2);
    //cout << "ratio: " << ratio << endl;
    if (ratio >= Random::nextDouble()) {
        //cout << "yes" << endl;
        //cout << m_particles[particle]->getPosition()[0] << endl;
        return true;
    }
    else {
        // correct position change
        //cout << "no" << endl;
        m_particles[particle]->adjustPosition(-change, dimension);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // measure cpu time
    clock_t start, finish;
    start = clock();
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        // euqilibrate
        if (i > numberOfMetropolisSteps*m_equilibrationFraction) {
            m_sampler->sample(acceptedStep);
        }

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


