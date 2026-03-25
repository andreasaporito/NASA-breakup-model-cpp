#include "Collision.h"

void Collision::init() {
    Breakup::init();
    //The pdf for Collisions is: 0.0101914/(x^2.71)
    _lcPowerLawExponent = -2.71;
    //Equation 12 mu = 0.9 * chi + 2.9
    _deltaVelocityFactorOffset = std::make_pair(0.9, 2.9);
}

void Collision::calculateFragmentCount() {
    using util::operator-, util::euclideanNorm;
    using util::operator/;
    //Get the two satellites from the input
    Satellite &sat1 = _input.at(0);
    Satellite &sat2 = _input.at(1);

    //Sets the maximalCharacteristicLength which will be required later
    _maximalCharacteristicLength = std::max(sat1.getCharacteristicLength(), sat2.getCharacteristicLength());
    //Sets the satType attribute to the correct type (later required for the A/M)
    //The Default of this member is SPACECRAFT
    if (sat1.getSatType() == SatType::ROCKET_BODY || sat2.getSatType() == SatType::ROCKET_BODY) {
        _satType = SatType::ROCKET_BODY;
    }

    //Assume sat1 is always the bigger one
    if (sat1.getCharacteristicLength() < sat2.getCharacteristicLength()) {
        std::swap(sat1, sat2);
    }
    _MassRatio = sat2.getMass() / sat1.getMass();
    _maxCharacteristicLength1 = sat1.getCharacteristicLength();
    _maxCharacteristicLength2 = sat2.getCharacteristicLength();
    //Sets the _input mass which will be required later for mass conservation purpose (maximal upper bound)
    _inputMass = sat1.getMass() + sat2.getMass();

    //Contains the mass M (later filled with an adequate value)
    double mass = 0;

    //The Relative Collision Velocity [m/s]
    const double dv = euclideanNorm(sat1.getVelocity() - sat2.getVelocity());
    // Squared Relative Collision Velocity [m^2/s^2]
    const double dv2 = dv * dv;

    //Calculate the Catastrophic Ratio, if greater than 40 J/g then we have a catastrophic collision
    //A catastrophic collision means that both satellites are fully fragmented whereas in a non-catastrophic collision
    //only the smaller satellite is fragmented (see here Section: Collision in [johnson et al.]
    double catastrophicRatio = (sat2.getMass() * dv2) / (2.0 * sat1.getMass() * 1000.0);
    if (catastrophicRatio < 40.0) {
        _isCatastrophic = false;
        // The original work states this as product of the projectile's mass in [kg] and the collision velocity in [km/s]
        // The recent paper below states that the original publication lacked the exponent 2 on the collision velocity
        // Horstman, A. (2020). Enhancement of s/c Fragmentation and Environment Evolution Models.
        // Final Report, Contract N. 4000115973/15/D/SR,
        // Institute of Space System, Technische Universität Braunschweig, 26(08).
        mass = sat2.getMass() * dv2 / 1e6;
    } else {
        _isCatastrophic = true;
        mass = sat1.getMass() + sat2.getMass();
    }

    //The fragment Count, respectively Equation 4
    auto fragmentCount = static_cast<size_t>(0.1 * std::pow(mass, 0.75) *
                                             std::pow(_minimalCharacteristicLength, -1.71));
    auto fragmentCount1 = fragmentCount * (1.0 / (1.0 + _MassRatio));
    auto fragmentCount2 = fragmentCount * (1.0 - (1.0 / (1.0 + _MassRatio)));
    this->generateFragments(fragmentCount1, sat1.getMass(), sat1.getVelocity(), sat1.getName(), sat1.getPosition());
    this->generateFragments(fragmentCount2, sat2.getMass(), sat2.getVelocity(), sat2.getName(), sat2.getPosition());
}
