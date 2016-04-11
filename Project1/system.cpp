#include "system.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <time.h>

using std::cout;
using std::endl;
using std::setprecision;
using std::vector;

bool System::metropolisStepImportance() {
    /* With importance sampling
     * The proposed change is now altered to make the particles drift
     * towards areas where the wavefunction is larger.
     */

    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle
    //int dimension = Random::nextInt(m_numberOfDimensions);  // choose random dimension

    // compute proposed change
    vector<double> plusChange = driftForce(particle);

    // compute proposed change
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        plusChange[dim] = 0.5*plusChange[dim]*m_timeStep + Random::nextGaussian(0.0, sqrt(m_timeStep));
    }

    // compute negative of proposed change
    vector<double> minusChange;
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        minusChange.push_back(-plusChange[dim]);
    }

    // get old wavefunction
    double waveFuncOld = m_waveFunction->evaluate(m_particles);

    // get old position
    vector<double> oldPosition = m_particles[particle]->getPosition();

    // get new position and wavefunction
    m_particles[particle]->adjustPositionAllDimensions(plusChange);
    vector<double> newPosition = m_particles[particle]->getPosition();
    double waveFuncNew = m_waveFunction->evaluate(m_particles);

    // get old Green's function, adjust back to old position first
    m_particles[particle]->adjustPositionAllDimensions(minusChange);
    double GreensOld = evaluateGreensFunction(particle, newPosition, oldPosition);

    // get new Green's function
    m_particles[particle]->adjustPositionAllDimensions(plusChange);
    double GreensNew = evaluateGreensFunction(particle, oldPosition, newPosition);

    // accept/reject new position using Metropolis algorithm
    double ratio = ( GreensNew*pow(waveFuncNew, 2) ) / ( GreensOld*pow(waveFuncOld, 2) );

    if (ratio >= Random::nextDouble()) {
        //cout << m_particles[particle]->getPosition()[0] << endl;
        return true;
    }
    else {
        // correct position change
        m_particles[particle]->adjustPositionAllDimensions(minusChange);
        return false;
    }
}

vector<double> System::driftForce(int particle) {
    // return a 3d "drift-vector" for the chosen particle

    vector<double> driftVector;
    for (int j=0; j < m_numberOfDimensions; j++) {
        driftVector.push_back(2*m_waveFunction->computeGradient(m_particles)[particle*m_numberOfDimensions + j]);
    }
    return driftVector;
}

double System::evaluateGreensFunction(int particle, vector<double> newPosition, vector<double> oldPosition) {

    vector<double> greensVector;

    // make vector that needs to be dotted
    for (int j=0; j < m_numberOfDimensions; j++) {
        greensVector.push_back(newPosition[j] - oldPosition[j] - 0.5*m_timeStep*driftForce(particle)[j]);
    }

    double greensFunction = 0;
    // find length squared of vector
    for (int j=0; j < m_numberOfDimensions; j++) {
        greensFunction += greensVector[j]*greensVector[j];
    }
    greensFunction /= 2*m_timeStep;
    return exp(-greensFunction);
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change its position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle
    int dimension = Random::nextInt(m_numberOfDimensions);  // choose random dimension
    double change = (Random::nextDouble()*2-1)*m_stepLength;  // propose change of particle's position

    // get old wavefunction
    double waveFuncOld = m_waveFunction->evaluate(m_particles);

    // adjust position
    m_particles[particle]->adjustPosition(change, dimension);

    // get new wavefunction
    double waveFuncNew = m_waveFunction->evaluate(m_particles);

    // accept/reject new position using Metropolis algorithm
    double ratio = pow(waveFuncNew, 2) / pow(waveFuncOld, 2);

    if (ratio >= Random::nextDouble()) {
        //cout << m_particles[particle]->getPosition()[0] << endl;
        return true;
    }
    else {
        // correct position change
        m_particles[particle]->adjustPosition(-change, dimension);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, bool useImportanceSampling,
                                bool writeEnergiesToFile, bool writePositionsToFile)
{
    m_particles                 = m_initialState->getParticles();
    m_numberOfAcceptedSteps = 0;
    if (m_samplerSetup == false) {
        m_sampler                   = new Sampler(this);
        m_samplerSetup = true;
    } else {
        m_sampler->clean();
    }
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // measure cpu time
    clock_t start, finish;
    start = clock();

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep;
        if (useImportanceSampling) { acceptedStep = metropolisStepImportance(); }
        else                       { acceptedStep = metropolisStep(); }

        // compute acceptance rate
        if (acceptedStep) {
            m_numberOfAcceptedSteps += 1;
        }

        // euqilibrate
        if (i > numberOfMetropolisSteps*m_equilibrationFraction) {
            m_sampler->sample(acceptedStep, writeEnergiesToFile, writePositionsToFile);
        }

        // print progression to terminal
        if ( !(i % 100) ) {
            cout << setprecision(2) << 100*(i / (double) numberOfMetropolisSteps) << " % complete" << "\r";
            fflush(stdout);
        }

    }
    finish = clock();
    cout << "Time elapsed: " << std::setprecision(5) << ((double) (finish-start)/CLOCKS_PER_SEC)
         << " s" << endl;

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

void System::setTimeStep(double timeStep) {
    m_timeStep = timeStep;
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


