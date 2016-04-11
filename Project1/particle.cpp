#include "particle.h"
#include <cassert>

Particle::Particle() {
}

void Particle::setPosition(const std::vector<double> &position) {
    // setter for position of each particle
    assert(position.size() == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
}

void Particle::adjustPositionAllDimensions(std::vector<double> change) {
    for (int j=0; j < m_numberOfDimensions; j++) {
        m_position.at(j) += change.at(j);
    }
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    // setter for number of dimensions for each particle (Particle class doesn't have access to System instance)
    m_numberOfDimensions = numberOfDimensions;
}
