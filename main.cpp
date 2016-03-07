#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinggaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorinteracting.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "steepestdescent.h"
#include "Math/random.h"

using namespace std;


int main() {

    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // oscillator frequency
    double alpha            = 0.7;          // variational parameter 1
    double beta             = 2.82843;      // variational parameter 2
    double stepLength       = 0.73;         // metropolis step length
    double equilibration    = 0.1;          // amount of the total steps used for equilibration
    double timeStep         = 0.01;         // importance sampling
    double a                = 0.0043;       // hard sphere radius
    double gamma            = 2.82843;      // trap potential strength z-direction

    bool useNumerical       = true;         // compute kinetic energy numerically
    bool useImportanceSampling = false;
    bool writeEnergiesToFile = true;
    bool writePositionsToFile = false;

    System* system = new System();
    //system->setHamiltonian              (new HarmonicOscillator(system, omega, useNumerical));
    system->setHamiltonian              (new HarmonicOscillatorInteracting(system, omega, a, gamma, useNumerical));
    //system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setWaveFunction             (new InteractingGaussian(system, alpha, beta, a));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (timeStep);
    system->runMetropolisSteps          (numberOfSteps, useImportanceSampling, writeEnergiesToFile,
                                         writePositionsToFile);

    // optimize alpha
    /*
    double stepLengthOptimize = 0.01;
    SteepestDescent* sd  = new SteepestDescent(system, stepLengthOptimize);
    sd->optimize(0.9);*/

    return 0;
}
