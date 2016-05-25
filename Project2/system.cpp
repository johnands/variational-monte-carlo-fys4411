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

bool System::metropolisStepSlater() {

    int particle = Random::nextInt(m_numberOfParticles);        // choose random particle
    int dimension = Random::nextInt(m_numberOfDimensions);      // choose random dimension
    double change = (Random::nextDouble()*2-1)*m_stepLength;    // propose change of particle's position
    /*vector<double> change(2);
    change[0] = (Random::nextDouble()*2-1)*m_stepLength;
    change[1] = (Random::nextDouble()*2-1)*m_stepLength;*/

    // store new proposed position
    m_particles[particle]->setNewPosition(change, dimension);
    //m_particles[particle]->setNewPositionAllDimensions(change);

    double ratio = m_waveFunction->computeRatio(m_particles, particle);
    //cout << setprecision(10) << ratio*ratio << endl;

    // this are the same for Slater
    if (ratio*ratio >= Random::nextDouble()) {
        //cout << "yes" << endl;
        m_waveFunction->updateSlaterInverse(m_particles, particle);
        m_particles[particle]->adjustPosition(change, dimension);
        //m_particles[particle]->adjustPositionAllDimensions(change);

        //double inverse = m_waveFunction->evaluate(m_particles);

        return true;
    }
    else {
        return false;
    }
}

bool System::metropolisStepSlaterImportance() {

    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle

    // calculate old quantum force
    vector<double> quantumForceOld = driftForce(particle);
    //cout << "old " <<quantumForceOld[0] << endl;
    //cout << "old " << quantumForceOld[1] << endl;

    // compute proposed new position
    vector<double> plusChange(m_numberOfDimensions);
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        plusChange[dim] = 0.5*quantumForceOld[dim]*m_timeStep +
                          Random::nextGaussian(0.0, sqrt(m_timeStep));
    }

    // compute negative of proposed change
    vector<double> minusChange;
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        minusChange.push_back(-plusChange[dim]);
    }

    // store new proposed position
    m_particles[particle]->setNewPositionAllDimensions(plusChange);

    // compute ratio (excluding Green's functions)
    double ratio = m_waveFunction->computeRatio(m_particles, particle);

    // get old position
    vector<double> oldPosition = m_particles[particle]->getPosition();

    // get new position
    vector<double> newPosition = m_particles[particle]->getNewPosition();

    // compute new quantum force (with updated position of chosen particle and ratio)
    // I need REAL position to be new to compute new drift force.....
    m_particles[particle]->adjustPositionAllDimensions(plusChange);
    vector<double> quantumForceNew = driftForce(particle);
    //cout << "new " << quantumForceNew[0] << endl;
    //cout << "new " << quantumForceNew[1] << endl;

    double exponent = 0;
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        double term1 = - (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[dim]) *
                         (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[dim]);
        double term2 =   (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[dim]) *
                         (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[dim]);
        exponent += term1 + term2;
    }

    /*for (int j=0; j < m_numberOfParticles; j++) {
        for (int dim=0; dim < m_numberOfDimensions; dim++) {
            double term1, term2;
            if (j != particle) {
                term1 = - (- 0.5*m_timeStep*quantumForceNew[2*j+dim]) *
                          (- 0.5*m_timeStep*quantumForceNew[2*j+dim]);
                term2 =   (- 0.5*m_timeStep*quantumForceOld[2*j+dim]) *
                          (- 0.5*m_timeStep*quantumForceOld[2*j+dim]);
            }
            else {
                term1 = - (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[2*j+dim]) *
                          (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[2*j+dim]);
                term2 =   (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[2*j+dim]) *
                          (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[2*j+dim]);
            }
            exponent += term1 + term2;

    }
    }*/
    double greensRatio = exp(exponent / 2*m_timeStep);

    // evaluate Greens' function with new and old position
    //double GreensOld = evaluateGreensFunction(newPosition, oldPosition, quantumForceOld);
    //double GreensNew = evaluateGreensFunction(oldPosition, newPosition, quantumForceNew);
    //cout << GreensOld << endl;
    //cout << GreensNew << endl;

    ratio *= ratio;
    ratio *= greensRatio;
    //cout << ratio << endl;

    if (ratio >= Random::nextDouble()) {
        m_waveFunction->updateSlaterInverse(m_particles, particle);       
        return true;
    }
    else {
        m_particles[particle]->adjustPositionAllDimensions(minusChange);
        return false;
    }
}

vector<double> System::driftForce(int particle) {
    // return a d-dimensional "drift vector" for the chosen particle

    vector<double> quantumForce = m_waveFunction->computeGradient(m_particles, particle);

    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        quantumForce[dim] *= 2;
    }

    return quantumForce;
}

double System::evaluateGreensFunction(vector<double> newPosition, vector<double> oldPosition,
                                      vector<double> quantumForce) {

    //vector<double> greensVector;
    double greensFunction = 0;
    // make vector that needs to be dotted
    /*for (int dim=0; dim < m_numberOfDimensions; dim++) {
        greensVector.push_back(newPosition[dim] - oldPosition[dim] - 0.5*m_timeStep*quantumForce[dim]);
    }*/

    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        greensFunction += - (newPosition[dim] - oldPosition[dim] - 0.5*m_timeStep*quantumForce[dim]) *
                            (newPosition[dim] - oldPosition[dim] - 0.5*m_timeStep*quantumForce[dim]);
    }

    // find length squared of vector
    /*for (int dim=0; dim < m_numberOfDimensions; dim++) {
        greensFunction += greensVector[dim]*greensVector[dim];
    }*/
    greensFunction /= 2*m_timeStep;
    return exp(greensFunction);
}

bool System::metropolisStepImportance() {
    /* With importance sampling
     * The proposed change is now altered to make the particles drift
     * towards areas where the wavefunction is larger.
     */

    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle
    //int dimension = Random::nextInt(m_numberOfDimensions);  // choose random dimension

    // calculate old quantum force using old position
    vector<double> quantumForceOld = driftForce(particle);

    // compute proposed change
    vector<double> plusChange(m_numberOfDimensions);
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        plusChange[dim] = 0.5*quantumForceOld[dim]*m_timeStep +
                Random::nextGaussian(0.0, m_timeStep);
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

    // get old Green's function using old quantum force
    double GreensOld = evaluateGreensFunction(newPosition, oldPosition, quantumForceOld);

    // get new Green's function using new quantum force and switching new and old position
    vector<double> quantumForceNew = driftForce(particle);
    double GreensNew = evaluateGreensFunction(oldPosition, newPosition, quantumForceNew);

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
        return true;
    }
    else {
        // correct position change
        m_particles[particle]->adjustPosition(-change, dimension);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps)
{
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

        if (getUseSlater()) {
            if (getUseImportanceSampling()) { acceptedStep = metropolisStepSlaterImportance(); }
            else                            { acceptedStep = metropolisStepSlater(); }
        }
        else {
            if (getUseImportanceSampling()) { acceptedStep = metropolisStepImportance(); }
            else                            { acceptedStep = metropolisStep(); }
        }

        // compute acceptance rate
        if (acceptedStep) {
            m_numberOfAcceptedSteps += 1;
        }

        // euqilibrate
        if (i > numberOfMetropolisSteps*m_equilibrationFraction) {
            m_sampler->sample(acceptedStep);
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

void System::setUseImportanceSampling(bool useImportanceSampling) {
    m_useImportanceSampling = useImportanceSampling;
}

void System::setWriteEnergiesToFile(bool writeEnergiesToFile) {
    m_writeEnergiesToFile = writeEnergiesToFile;
}

void System::setWritePositionsToFile(bool writePositionsToFile) {
    m_writePositionsToFile = writePositionsToFile;
}

void System::setTimeStep(double timeStep) {
    m_timeStep = timeStep;
}

void System::setUseSlater(bool useSlater) {
    m_useSlater = useSlater;
}

void System::setUseJastrow(bool useJastrow) {
    m_useJastrow = useJastrow;
}

void System::setOptimizeParameters(bool optimizeParameters) {
    m_optimizeParameters = optimizeParameters;
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

void System::setParticles(std::vector<Particle*> particles) {
    m_particles = particles;
}

void System::setParallel(bool parallel) {
    m_parallel = parallel;
}

void System::setRank(int rank) {
    m_rank = rank;
}

void System::setSize(int size) {
    m_size = size;
}

