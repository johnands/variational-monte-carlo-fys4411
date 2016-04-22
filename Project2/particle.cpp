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
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        m_position.at(dim) += change.at(dim);
    }
}

void Particle::setNewPosition(double change, int dimension) {
    m_newPosition = m_position;
    m_newPosition.at(dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
