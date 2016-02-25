#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorinteracting.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

using namespace std;


int main() {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 5;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency
    double alpha            = 0.6;          // Variational parameter
    double stepLength       = 0.73;         // Metropolis step length
    double equilibration    = 0.1;          // Amount of the total steps used
    bool useNumerical       = false;        // compute kinetic energy numerically
    double timeStep         = 0.01;         // used in importance sampling
    bool useImportanceSampling = true;


    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, useNumerical));
    //system->setHamiltonian              (new HarmonicOscillatorInteracting(system, omega, useNumerical));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (timeStep);
    system->runMetropolisSteps          (numberOfSteps, useImportanceSampling);
    return 0;
}
