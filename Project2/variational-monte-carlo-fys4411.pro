TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    Hamiltonians/harmonicoscillatorinteracting.cpp \
    WaveFunctions/interactinggaussian.cpp \
    steepestdescent.cpp \
    WaveFunctions/quantumdottwoelectrons.cpp \
    Hamiltonians/harmonicoscillatorquantumdot2.cpp \
    WaveFunctions/manybodyquantumdot.cpp \
    Hamiltonians/homanybodyquantumdot.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    Hamiltonians/harmonicoscillatorinteracting.h \
    WaveFunctions/interactinggaussian.h \
    steepestdescent.h \
    WaveFunctions/quantumdottwoelectrons.h \
    Hamiltonians/harmonicoscillatorquantumdot2.h \
    WaveFunctions/manybodyquantumdot.h \
    Hamiltonians/homanybodyquantumdot.h

LIBS += -larmadillo -lblas -llapack

#release {
#    DEFINES += ARMA_NO_DEBUG
#    QMAKE_CXXFLAGS_RELEASE -= -O2
#    QMAKE_CXXFLAGS_RELEASE += -O3
#}

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK


