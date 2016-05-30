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
#include <armadillo>
#include <mpi.h>

using std::cout;
using std::endl;
using std::setprecision;
using std::vector;

bool System::metropolisStepSlater() {

    int particle = Random::nextInt(m_numberOfParticles);        // choose random particle
    int dimension = Random::nextInt(m_numberOfDimensions);      // choose random dimension
    double change = (Random::nextDouble()*2-1)*m_stepLength;    // propose change of particle's position

    // store new proposed position
    m_particles[particle]->setNewPosition(change, dimension);

    // compute ratio of wave functions
    double ratio = m_waveFunction->computeRatio(m_particles, particle);

    // Metropolis test
    if (ratio*ratio >= Random::nextDouble()) {
        m_waveFunction->updateSlaterInverse(m_particles, particle);
        m_particles[particle]->adjustPosition(change, dimension);
        return true;
    }
    else {
        return false;
    }
}

bool System::metropolisStepSlaterImportance() {

    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle

    arma::mat quantumForceOld = arma::zeros<arma::mat>(m_numberOfParticles, m_numberOfDimensions);
    arma::mat quantumForceNew = quantumForceOld;
    arma::mat oldPositions    = quantumForceOld;

    setUseOldQuantumForce(true);
    for (int i=0; i < m_numberOfParticles; i++){
        // calculate drift force for particle i
        vector<double> force = driftForce(i);
        for (int j=0; j < m_numberOfDimensions; j++){
            // store in matrix
            quantumForceOld(i,j) = force[j];
            oldPositions(i,j) = m_particles[i]->getPosition()[j];
        }
    }

    // compute proposed new position
    vector<double> plusChange(m_numberOfDimensions);
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        plusChange[dim] = 0.5*quantumForceOld(particle,dim)*m_timeStep +
                          Random::nextGaussian(0.0, sqrt(m_timeStep));
    }

    //vector<double> quantumForceOld = driftForce(particle);

    // compute negative of proposed change
    vector<double> minusChange;
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        minusChange.push_back(-plusChange[dim]);
    }

    // store new proposed position
    m_particles[particle]->setNewPositionAllDimensions(plusChange);

    // compute ratio of wave functions
    double ratio = m_waveFunction->computeRatio(m_particles, particle);

    // get old position
    //vector<double> oldPosition = m_particles[particle]->getPosition();

    // get new position
    //vector<double> newPosition = m_particles[particle]->getNewPosition();

    // compute new quantum force with new position
    m_particles[particle]->adjustPositionAllDimensions(plusChange);
    //vector<double> quantumForceNew = driftForce(particle);

    setUseOldQuantumForce(false);
    double exponent = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        vector<double> force = driftForce(i);
        // here, I must divide by Rsd
        for (int j=0; j < m_numberOfDimensions; j++) {
            quantumForceNew(i,j) = force[j];
            double term1 = - (oldPositions(i,j) - m_particles[i]->getPosition()[j] - 0.5*m_timeStep*quantumForceNew(i,j)) *
                             (oldPositions(i,j) - m_particles[i]->getPosition()[j] - 0.5*m_timeStep*quantumForceNew(i,j));
            double term2 = (- oldPositions(i,j) + m_particles[i]->getPosition()[j] - 0.5*m_timeStep*quantumForceOld(i,j)) *
                           (- oldPositions(i,j) + m_particles[i]->getPosition()[j] - 0.5*m_timeStep*quantumForceOld(i,j));
            exponent += term1 + term2;
        }
    }

    double greensRatio = exp(exponent / (2*m_timeStep));
    // calculate old quantum force


    // compute ratio of new and old Greens' function
    //double greensRatio = evaluateGreensFunction(oldPosition, newPosition, quantumForceOld, quantumForceNew);

    // calculate final ratio
    ratio *= ratio;
    ratio *= greensRatio;

    // Metropolis test
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

double System::evaluateGreensFunction(vector<double> oldPosition, vector<double> newPosition,
                                      vector<double> quantumForceOld, vector<double> quantumForceNew) {

    double exponent = 0;
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        double term1 = - (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[dim]) *
                         (oldPosition[dim] - newPosition[dim] - 0.5*m_timeStep*quantumForceNew[dim]);
        double term2 =   (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[dim]) *
                         (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*quantumForceOld[dim]);
        exponent += term1 + term2;
    }
    return exp(exponent / 2*m_timeStep);
}

bool System::metropolisStepImportance() {
    /* With importance sampling
     * The proposed change is now altered to make the particles drift
     * towards areas where the wavefunction is larger.
     */

    int particle = Random::nextInt(m_numberOfParticles);    // choose random particle

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

    // evaluate old wavefunction
    double waveFuncOld = m_waveFunction->evaluate(m_particles);

    // get old position
    vector<double> oldPosition = m_particles[particle]->getPosition();

    // get new position and wavefunction
    m_particles[particle]->adjustPositionAllDimensions(plusChange);
    vector<double> newPosition = m_particles[particle]->getPosition();
    double waveFuncNew = m_waveFunction->evaluate(m_particles);

    // calculate new drift force with new position
    vector<double> quantumForceNew = driftForce(particle);

    // compute ratio of new and old Greens' function
    double greensRatio = evaluateGreensFunction(oldPosition, newPosition, quantumForceOld, quantumForceNew);

    // accept/reject new position using Metropolis algorithm
    double ratio = ( waveFuncNew*waveFuncNew * greensRatio ) / ( waveFuncOld*waveFuncOld );

    if (ratio >= Random::nextDouble()) {
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
    double ratio = ( waveFuncNew*waveFuncNew ) / ( waveFuncOld*waveFuncOld );

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
    }
    else {
        m_sampler->clean();
    }
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // measure cpu time
    //clock_t start, finish;
    //start = clock();
    double startTime = MPI_Wtime();

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
    //finish = clock();
    double endTime = MPI_Wtime();
    /*cout << "Time elapsed: " << std::setprecision(5) << ((double) (finish-start)/CLOCKS_PER_SEC)
         << " s" << endl;*/
    double totalTime = endTime - startTime;


    m_sampler->computeAverages();
    if (m_parallel) {
        if (m_rank == 0) {
            cout << "Time elapsed: " << std::setprecision(5) << totalTime <<
                    " on number of processors: " << m_size << endl;
            m_sampler->printOutputToTerminal();
        }
    }
    else {
        m_sampler->printOutputToTerminal();
    }
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

void System::setRatioSD(double ratioSD) {
    m_ratioSD = ratioSD;
}

void System::setUseOldQuantumForce(bool useOldQuantumForce) {
    m_useOldQuantumForce = useOldQuantumForce;
}

