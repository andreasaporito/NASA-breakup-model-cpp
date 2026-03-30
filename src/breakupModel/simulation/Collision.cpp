#include "Collision.h"

void Collision::run() {
    // init() is called once here at the Collision level; each SubCollision
    // will also call it via Breakup::run(), which is correct since they each
    // need their own fully initialised state.
    this->init();

    // Ensure sat1 is always the bigger one (same invariant as before)
    if (_input.at(0).getCharacteristicLength() < _input.at(1).getCharacteristicLength()) {
        std::swap(_input.at(0), _input.at(1));
    }

    auto sub1 = std::make_unique<SubCollision>(
        _input, 1, _minimalCharacteristicLength, 
        _currentMaxGivenID, _enforceMassConservation
    );

    auto sub2 = std::make_unique<SubCollision>(
        _input, 2, _minimalCharacteristicLength, 
        _currentMaxGivenID, _enforceMassConservation
    );

    sub1->run();
    sub2->run();

    // Merge both fragment sets into this Collision's _output
    const size_t size1 = sub1->getResultSoA().size();
    const size_t size2 = sub2->getResultSoA().size();
    _output.resize(size1 + size2);
    _output.copyFrom(sub1->getResultSoA(), 0);
    _output.copyFrom(sub2->getResultSoA(), size1);
}

void Collision::init() {
    Breakup::init();
    //The pdf for Collisions is: 0.0101914/(x^2.71)
    _lcPowerLawExponent = -2.71;
    //Equation 12 mu = 0.9 * chi + 2.9
    _deltaVelocityFactorOffset = std::make_pair(0.9, 2.9);
}

void Collision::addFurtherFragments() {
    if (!_isCatastrophic) {
        // If non-catastrophic: add a remainder fragment
        Satellite &target = _input.at(0);

        // Prepend, so that the bigger satellite (target) is assigned as parent
        auto tuple = _output.prependElement();
        auto &[lc, areaToMassRatio, area, mass] = tuple;

        // One special fragment representing the cratered remainder of the target satellite hit by the projectile
        // However, it should not be heavier than the actual original parent!
        mass = std::min(_inputMass - _outputMass, target.getMass());
        lc = util::calculateCharacteristicLengthFromMass(mass);
        areaToMassRatio = calculateAreaMassRatio(lc);
        area = calculateArea(lc);

        // Update the output mass accordingly
        _outputMass += mass;
    }
    Breakup::addFurtherFragments();
}

void Collision::calculateFragmentCount() {}
void Collision::assignParentProperties() {}
void Collision::enforceMassConservation() {}

SubCollision::SubCollision(std::vector<Satellite> input,
                           int cardinality,
                           double minimalCharacteristicLength,
                           size_t currentMaxGivenID,
                           bool enforceMassConservation)
    : Collision(std::move(input), minimalCharacteristicLength,
                currentMaxGivenID, enforceMassConservation),
      _cardinality{cardinality}
{}

void SubCollision::run() {
    spdlog::info("Running SubCollision with cardinality {}.", _cardinality);
    Breakup::run();
}

void SubCollision::calculateFragmentCount() {
    using util::operator-, util::euclideanNorm;
    using util::operator/;
    //Get the two satellites from the input
    Satellite &sat1 = _input.at(0);
    Satellite &sat2 = _input.at(1);

    //Sets the maximalCharacteristicLength which will be required later
    _maximalCharacteristicLength =
        (_cardinality == 1) ? sat1.getCharacteristicLength()
                            : sat2.getCharacteristicLength();
    //Sets the satType attribute to the correct type (later required for the A/M)
    //The Default of this member is SPACECRAFT
    if (sat1.getSatType() == SatType::ROCKET_BODY || sat2.getSatType() == SatType::ROCKET_BODY) {
        _satType = SatType::ROCKET_BODY;
    }
    Satellite parent;
    //Sets the _input mass which will be required later for mass conservation purpose (maximal upper bound)
    _inputMass = (_cardinality == 1) ? sat1.getMass() : sat2.getMass();

    //Contains the mass M (later filled with an adequate value)
    double mass = 0;
    double catastrophicRatio = 0;
    //The Relative Collision Velocity [m/s]
    const double dv = euclideanNorm(sat1.getVelocity() - sat2.getVelocity());
    // Squared Relative Collision Velocity [m^2/s^2]
    const double dv2 = dv * dv;
    //just one branch for performance reasons, since the only difference is which satellite is the parent and which one is the projectile, but the rest of the calculations are identical
    if (_cardinality == 1) {
        _maximalCharacteristicLength = sat1.getCharacteristicLength();
        _inputMass = sat1.getMass();
        catastrophicRatio = (sat2.getMass() * dv2) / (sat1.getMass() * 2 * 1000.0);
        const Satellite &parent = sat1;
        //Calculate the Catastrophic Ratio, if greater than 40 J/g then we have a catastrophic collision
        //A catastrophic collision means that both satellites are fully fragmented whereas in a non-catastrophic collision
        //only the smaller satellite is fragmented (see here Section: Collision in [johnson et al.]
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
            mass = sat1.getMass();
        }
    } else {
        _maximalCharacteristicLength = sat2.getCharacteristicLength();
        _inputMass = sat2.getMass();
        catastrophicRatio = (sat1.getMass() * dv2) / (sat2.getMass() * 2 * 1000.0);
        const Satellite &parent = sat2;
        if (catastrophicRatio < 40.0) {
            _isCatastrophic = false;
            // The original work states this as product of the projectile's mass in [kg] and the collision velocity in [km/s]
            // The recent paper below states that the original publication lacked the exponent 2 on the collision velocity
            // Horstman, A. (2020). Enhancement of s/c Fragmentation and Environment Evolution Models.
            // Final Report, Contract N. 4000115973/15/D/SR,
            // Institute of Space System, Technische Universität Braunschweig, 26(08).
            mass = sat1.getMass() * dv2 / 1e6;
        } else {
            _isCatastrophic = true;
            mass = sat2.getMass();
        }
    }

    //The fragment Count, respectively Equation 4
    auto fragmentCount = static_cast<size_t>(0.1 * std::pow(mass, 0.75) *
                                             std::pow(_minimalCharacteristicLength, -1.71));
    this->generateFragments(fragmentCount, parent.getPosition());
}

void SubCollision::assignParentProperties() {
    // All fragments in this sub-simulation come from one satellite only
    const Satellite& parent = (_cardinality == 1) ? _input.at(0) : _input.at(1);
    auto debrisName = std::make_shared<const std::string>(
        parent.getName() + "-Collision-Fragment");

    auto tupleView = _output.getCMVNTuple();
    std::for_each(tupleView.begin(), tupleView.end(),
                  [&](auto& tuple) {
                      auto& [lc, mass, velocity, name] = tuple;
                      name     = debrisName;
                      velocity = parent.getVelocity();
                  });
}

void SubCollision::addFurtherFragments() {
    // Both satellites may need remainder fragments; Collision::addFurtherFragments()
    // handles the _isCatastrophic gate internally
    Collision::addFurtherFragments();
}

void SubCollision::enforceMassConservation() {
    //Enforce Mass Conservation if the output mass is greater than the input mass
    _outputMass = std::reduce(std::execution::par_unseq,_output.mass.begin(), _output.mass.end(), 0.0);
    spdlog::debug("The simulation got {} kg of input mass for fragments", _inputMass);
    spdlog::debug("The simulation produced {} kg of debris", _outputMass);
    size_t oldSize = _output.size();
    size_t newSize = oldSize ;
    const double epsilon = 1e-6;

    // Build a permutation index sorted by mass descending
    std::vector<size_t> idx(newSize);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return _output.mass[a] < _output.mass[b]; // ascending order
    });

    // Apply the permutation in-place to every per-element vector in the SoA
    auto applyPermutation = [&](auto& vec) {
        using T = typename std::decay_t<decltype(vec)>::value_type;
        std::vector<T> tmp(vec.size());
        for (size_t i = 0; i < idx.size(); ++i) { tmp[i] = vec[idx[i]]; }
        vec = std::move(tmp);
    };
    applyPermutation(_output.name);
    applyPermutation(_output.characteristicLength);
    applyPermutation(_output.areaToMassRatio);
    applyPermutation(_output.mass);
    applyPermutation(_output.area);
    applyPermutation(_output.ejectionVelocity);
    applyPermutation(_output.velocity);

    if (_outputMass > _inputMass + epsilon) {

        while (_outputMass > _inputMass && !_output.mass.empty()) {
            _outputMass -= _output.mass.back();
            _output.popBack();
            newSize -= 1;
        }

        double gap = _inputMass - _outputMass;
        if (gap > epsilon && _enforceMassConservation) {
            this->addSingleFragment(gap); 
        }
    } 
    else if (_inputMass > _outputMass + epsilon && _enforceMassConservation) {
        if (_enforceMassConservation) {
            this->addFurtherFragments();
            
            if (_output.size() != oldSize) {
                this->enforceMassConservation();
            }
        }
    }

    // Some helpful logging hints
    if (oldSize != newSize) {
        spdlog::warn("The simulation modified the number of fragments to enforce the mass conservation in SubCollision {}.", _cardinality);
        spdlog::warn("The fragment count was adapted from {} to {} fragments.", oldSize, newSize);
    }
    spdlog::warn("Initial mass was {} kg, final mass is {} kg.", _inputMass, _outputMass);
}
