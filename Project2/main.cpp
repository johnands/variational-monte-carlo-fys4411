#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinggaussian.h"
#include "WaveFunctions/quantumdottwoelectrons.h"
#include "WaveFunctions/manybodyquantumdot.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorinteracting.h"
#include "Hamiltonians/harmonicoscillatorquantumdot2.h"
#include "Hamiltonians/homanybodyquantumdot.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "steepestdescent.h"
#include "Math/random.h"

using namespace std;

int main() {

    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // oscillator frequency
    double alpha            = 1.0;          // variational parameter 1
    double beta             = 0.5;        // variational parameter 2
    double stepLength       = 1.2;          // metropolis step length
    double equilibration    = 0.1;          // amount of the total steps used for equilibration
    double timeStep         = 0.005;        // importance sampling
    double a                = 0.0043;       // hard sphere radius
    double a2               = 0;          // two-body quantum dot (depends on spin)
    double gamma            = 2.82843;      // trap potential strength z-direction

    bool useNumerical           = false;    // compute kinetic energy numerically
    bool useImportanceSampling  = false;
    bool writeEnergiesToFile    = false;
    bool writePositionsToFile   = false;

    System* system = new System();   
    system->setUseJastrow               (true);
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    //system->setHamiltonian              (new HarmonicOscillator(system, omega, useNumerical));
    //system->setHamiltonian              (new HarmonicOscillatorInteracting(system, omega, a, gamma, useNumerical));
    //system->setHamiltonian              (new HarmonicOscillatorQuantumDot2(system, omega, useNumerical));
    system->setHamiltonian              (new HOManyBodyQuantumDot(system, omega, useNumerical));

    //system->setWaveFunction             (new SimpleGaussian(system, alpha));
    //system->setWaveFunction             (new InteractingGaussian(system, alpha, beta, a));
    //system->setWaveFunction             (new QuantumDotTwoElectrons(system, alpha, beta, omega, a2));
    system->setWaveFunction             (new ManyBodyQuantumDot(system, alpha, beta, omega));

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (timeStep);
    system->setUseSlater                (true);
    system->runMetropolisSteps          (numberOfSteps, useImportanceSampling, writeEnergiesToFile, writePositionsToFile);

    // optimize alpha
    /*double initialAlpha = 0.5;
    double initialBeta = 0.6;
    std::vector<double> parameters(2);
    parameters[0] = initialAlpha; parameters[1] = initialBeta;
    system->setOptimizeParameters(true);
    double stepLengthOptimize = 0.1;
    SteepestDescent* sd  = new SteepestDescent(system, stepLengthOptimize);
    sd->optimize(parameters);*/

    return 0;
}
