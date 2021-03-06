#pragma once
#include <vector>

class System {
public:
    bool metropolisStep                 ();
    bool metropolisStepImportance       ();
    bool metropolisStepSlater           ();
    bool metropolisStepSlaterImportance ();
    std::vector<double> driftForce  (int particle);
    double evaluateGreensFunction   (std::vector<double> oldPosition, std::vector<double> newPosition,
                                     std::vector<double> quantumForceOld, std::vector<double> quantumForceNew);
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setTimeStep                (double timeStep);
    void setUseImportanceSampling   (bool useImportanceSampling);
    void setWriteEnergiesToFile     (bool writeEnergiesToFile);
    void setWritePositionsToFile    (bool writePositionsToFile);
    void setUseSlater               (bool useSlater);
    void setUseJastrow              (bool useJastrow);
    void setOptimizeParameters      (bool optimizeParameters);
    void setParallel                (bool parallel);
    void setRank                    (int rank);
    void setSize                    (int size);
    void setRatioSD                 (double ratioSD);
    void setUseOldQuantumForce      (bool useOldQuantumForce);

    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setParticles               (std::vector<class Particle*> particles);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class InitialState*             getInitialState()   { return m_initialState; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    int getNumberOfAcceptedSteps()      { return m_numberOfAcceptedSteps; }
    double getTimeStep()                { return m_timeStep; }
    bool getUseImportanceSampling()     { return m_useImportanceSampling; }
    bool getWriteEnergiesToFile()       { return m_writeEnergiesToFile; }
    bool getWritePositionsToFile()      { return m_writePositionsToFile; }
    bool getUseSlater()                 { return m_useSlater; }
    bool getUseJastrow()                { return m_useJastrow; }
    bool getOptimizeParameters()        { return m_optimizeParameters; }
    bool getParallel()                  { return m_parallel; }
    int getRank()                       { return m_rank; }
    int getSize()                       { return m_size; }
    double getRatioSD()                 { return m_ratioSD; }
    bool getUseOldQuantumForce()        { return m_useOldQuantumForce; }

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_timeStep = 0.01;
    int                             m_numberOfAcceptedSteps = 0;
    bool                            m_useImportanceSampling = false;
    bool                            m_writeEnergiesToFile = false;
    bool                            m_writePositionsToFile = false;
    bool                            m_samplerSetup = false;
    bool                            m_useSlater = false;
    bool                            m_useJastrow = true;
    bool                            m_optimizeParameters = false;
    bool                            m_parallel = false;

    int                             m_rank;
    int                             m_size;

    double                          m_ratioSD;
    bool                            m_useOldQuantumForce = true;

    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};

