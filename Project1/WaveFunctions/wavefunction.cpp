#include "wavefunction.h"


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setAlpha(double alpha) {
    m_parameters[0] = alpha;
}
