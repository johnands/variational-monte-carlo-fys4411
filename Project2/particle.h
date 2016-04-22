#pragma once
#include <vector>

class Particle {
public:
    Particle();
    void setPosition(const std::vector<double> &position);
    void adjustPosition(double change, int dimension);
    void adjustPositionAllDimensions(std::vector<double> change);
    void setNumberOfDimensions(int numberOfDimensions);
    void setNewPosition(double change, int dimension);
    std::vector<double>& getPosition() { return m_position; }
    std::vector<double>& getNewPosition() {return m_newPosition; }

private:
    int     m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
    std::vector<double> m_newPosition = std::vector<double>();
};

