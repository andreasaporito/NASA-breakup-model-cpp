#include "Breakup.h"
#include "breakupModel/util/UtilityFunctions.h"

void Breakup::run() {
    //0. Step: Prepare constants, etc.
    this->init();

    //1. Step: Generate the new Satellites
    this->calculateFragmentCount();

    //2. Step: Assign every new Satellite a value for L_c
    this->characteristicLengthDistribution();

    //3. Step: Calculate the A/M (area-to-mass-ratio), A (area) and M (mass) values for every Satellite
    this->areaToMassRatioDistribution();

    //4. Step: Enforce the Mass Conservation and remove (or add) fragments
    this->enforceMassConservation();

    //5. Step: Assign parent and by doing that assign each fragment a base velocity
    this->assignParentProperties();

    //6. Step: Calculate the Ejection velocity for every Satellite
    this->deltaVelocityDistribution();

    this->enforceKineticEnergyConservation();

    this->enforceMomentumConservation();

    //7. Step: As a last step set the _currentMaxGivenID to the new valid value
    _currentMaxGivenID += _output.size();

    // Let's output all the mass, kinetic energy and momentum conservation values here
    spdlog::debug("Initial mass was {} kg, final mass is {} kg.", _inputMass, _outputMass);

    double outputKineticEnergy = std::transform_reduce(
        std::execution::par_unseq,
        _output.mass.begin(), 
        _output.mass.end(),
        _output.velocity.begin(),
        0.0,
        std::plus<>(),
        [](double m, const std::array<double, 3>& v) {
            return util::calculateKineticEnergy(m, v);
        }
    );
    // Add the percentage difference in kinetic energy to the log message
    double kineticEnergyDifference = std::abs(outputKineticEnergy - _initialKineticEnergy);
    double kineticEnergyDifferencePercent = (kineticEnergyDifference / (_initialKineticEnergy + 1e-10)) * 100.0;
    spdlog::debug("Initial kinetic energy was {} J, final kinetic energy is {} J. Difference is {} J ({}%).", 
                 _initialKineticEnergy, outputKineticEnergy, kineticEnergyDifference, kineticEnergyDifferencePercent);

    std::array<double, 3> momentumDifference = {calculateCurrentMomentum()[0] - _initialMomentum[0], calculateCurrentMomentum()[1] - _initialMomentum[1], calculateCurrentMomentum()[2] - _initialMomentum[2]};
    double errorMagnitude = std::sqrt(momentumDifference[0]*momentumDifference[0] + momentumDifference[1]*momentumDifference[1] + momentumDifference[2]*momentumDifference[2]);
    double initialMomentumMagnitude = std::sqrt(_initialMomentum[0]*_initialMomentum[0] + _initialMomentum[1]*_initialMomentum[1] + _initialMomentum[2]*_initialMomentum[2]);
    double momentumDifferencePercent = (errorMagnitude / (initialMomentumMagnitude + 1e-10)) * 100.0;
    spdlog::debug("Initial momentum was [{}, {}, {}] kg*m/s, final momentum is [{}, {}, {}] kg*m/s, difference is {} kg*m/s ({}%).",
                 _initialMomentum[0], _initialMomentum[1], _initialMomentum[2],
                 calculateCurrentMomentum()[0], calculateCurrentMomentum()[1], calculateCurrentMomentum()[2],
                 errorMagnitude, momentumDifferencePercent);

}

Breakup &Breakup::setSeed(std::optional<unsigned long> seed) {
    if (seed.has_value()) {
        _fixRNG = std::mt19937 {seed.value()};
    } else {
        _fixRNG = std::nullopt;
    }
    return *this;
}

void Breakup::init() {
    _inputMass = 0;
    _outputMass = 0;
}

void Breakup::generateFragments(size_t fragmentCount, const std::array<double, 3> &position) {
    //_ialPosition = position;
    _output = Satellites{_currentMaxGivenID+1, SatType::DEBRIS, position, fragmentCount};
}

void Breakup::characteristicLengthDistribution() {
    std::for_each(std::execution::par_unseq, _output.characteristicLength.begin(), _output.characteristicLength.end(),
                  [&](double &lc) {
        lc = calculateCharacteristicLength();
    });
}

void Breakup::areaToMassRatioDistribution() {
    auto tupleView = _output.getAreaMassTuple();
    std::for_each(std::execution::par_unseq, tupleView.begin(), tupleView.end(),
                  [&](auto &tuple) {
        //Order in the tuple: 0: L_c | 1: A/M | 2: Area | 3: Mass
        auto &[lc, areaToMassRatio, area, mass] = tuple;
        //Calculate the A/M value in [m^2/kg]
        areaToMassRatio = calculateAreaMassRatio(lc);
        //Calculate the area A in [m^2]
        area = calculateArea(lc);
        //Calculate the mass m in [kg]
        mass = calculateMass(area, areaToMassRatio);
    });
}
// ----------------- Still to revise a bit not sure wheter it is the fastest way to navigate --------------------
void Breakup::enforceMassConservation() {
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
        spdlog::debug("The fragment count was adapted from {} to {} fragments.", oldSize, newSize);
    }
}

void Breakup::addSingleFragment(double mass) {
    auto tuple = _output.appendElement();
    auto& [lc, amr, area, m] = tuple;

    m = mass;
    lc = util::calculateCharacteristicLengthFromMass(m);
    area = calculateArea(lc);
    amr  = calculateAreaMassRatio(lc);

    _outputMass += m;
}

void Breakup::addFurtherFragments() {
    while (_outputMass < _inputMass) {
        //Order in the tuple: 0: L_c | 1: A/M | 2: Area | 3: Mass
        //Create new element and assign values
        auto tuple = _output.appendElement();
        auto &[lc, areaToMassRatio, area, mass] = tuple;
        lc = calculateCharacteristicLength();
        areaToMassRatio = calculateAreaMassRatio(lc);
        area = calculateArea(lc);
        mass = calculateMass(area, areaToMassRatio);

        //Calculate new mass
        _outputMass += mass;
    }
    //Remove the element which has lead to the exceeding of the mass budget
    _outputMass -= _output.mass.back();
    _output.popBack();
}

void Breakup::deltaVelocityDistribution() {
    using namespace util;
    auto tupleView = _output.getVelocityTuple();
    std::for_each(std::execution::par_unseq, tupleView.begin(), tupleView.end(),
                  [&](auto &tuple) {
        //Order in the tuple: 0: A/M | 1: Velocity | 2: Ejection Velocity
        auto &[areaToMassRatio, velocity, ejectionVelocity] = tuple;
        //Calculates the velocity as a scalar based on Equation 11/ 12
        const double chi = log10(areaToMassRatio);
        const double mu = _deltaVelocityFactorOffset.first * chi + _deltaVelocityFactorOffset.second;
        constexpr double sigma = 0.4;
        std::normal_distribution<> normalDistribution{mu, sigma};
        double velocityScalar = std::pow(10.0, getRandomNumber(normalDistribution));

        //Transform the scalar velocity into a cartesian vector
        ejectionVelocity = calculateVelocityVector(velocityScalar);
        velocity = velocity + ejectionVelocity;
    });
}

double Breakup::calculateCharacteristicLength() {
    using util::transformUniformToPowerLaw;
    static std::uniform_real_distribution<> uniformRealDistribution{0.0, 1.0};
    const double y = getRandomNumber(uniformRealDistribution);
    return transformUniformToPowerLaw(_minimalCharacteristicLength, _maximalCharacteristicLength, _lcPowerLawExponent, y);
}

double Breakup::calculateAreaMassRatio(double characteristicLength) {
    using namespace util;
    const double logLc = std::log10(characteristicLength);

    if (characteristicLength > 0.11) {
        //Case bigger than 11 cm
        std::normal_distribution<> n1{mu_1(_satType, logLc), sigma_1(_satType, logLc)};
        std::normal_distribution<> n2{mu_2(_satType, logLc), sigma_2(_satType, logLc)};

        return std::pow(10.0, alpha(_satType, logLc) * getRandomNumber(n1) +
            (1 - alpha(_satType, logLc)) * getRandomNumber(n2));
    } else if (characteristicLength < 0.08) {
        //Case smaller than 8 cm
        std::normal_distribution<> n{mu_soc(logLc), sigma_soc(logLc)};

        return std::pow(10.0, getRandomNumber(n));
    } else {
        //Case between 8 cm and 11 cm
        std::normal_distribution<> n1{mu_1(_satType, logLc), sigma_1(_satType, logLc)};
        std::normal_distribution<> n2{mu_2(_satType, logLc), sigma_2(_satType, logLc)};
        std::normal_distribution<> n{mu_soc(logLc), sigma_soc(logLc)};

        double y1 = std::pow(10.0, alpha(_satType, logLc) * getRandomNumber(n1) +
                                 (1.0 - alpha(_satType, logLc)) * getRandomNumber(n2));
        double y0 = std::pow(10.0, getRandomNumber(n));

        //beta * y1 + (1 - beta) * y0 = beta * y1 + y0 - beta * y0 = y0 + beta * (y1 - y0)
        return y0 + (characteristicLength - 0.08) * (y1 - y0) / (0.03);
    }
}

double Breakup::calculateArea(double characteristicLength) {
    constexpr double lcBound = 0.00167;
    if (characteristicLength < lcBound) {
        constexpr double factorLittle = 0.540424;
        return factorLittle * characteristicLength * characteristicLength;
    } else {
        constexpr double exponentBig = 2.0047077;
        constexpr double factorBig = 0.556945;
        return factorBig * std::pow(characteristicLength, exponentBig);
    }
}

double Breakup::calculateMass(double area, double areaMassRatio) {
    return area / areaMassRatio;
}

std::array<double, 3> Breakup::calculateVelocityVector(double velocity) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::uniform_real_distribution<double> distPhi(0.0, 2.0 * util::PI);
    double z = getRandomNumber(dist);
    double phi = getRandomNumber(distPhi);
    double xy_radius = std::sqrt(1 - z * z);

    return {
        velocity * xy_radius * std::cos(phi),
        velocity * xy_radius * std::sin(phi),
        velocity * z
    };
}

std::array<double, 3> Breakup::calculateCurrentMomentum(){
    auto output_momentum = std::transform_reduce(
        std::execution::par_unseq,
        _output.mass.begin(),
        _output.mass.end(),
        _output.velocity.begin(),
        std::array<double, 3>{0.0, 0.0, 0.0},
        [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
            return std::array<double, 3>{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        },
        [](double m, const std::array<double, 3>& v) {
            return std::array<double, 3>{m * v[0], m * v[1], m * v[2]};
        });
    return output_momentum;
}

void Breakup::enforceMomentumConservation() {
    auto currentP = calculateCurrentMomentum();
    // Momentum error vector
    std::array<double, 3> error = {
        currentP[0] - _initialMomentum[0],
        currentP[1] - _initialMomentum[1],
        currentP[2] - _initialMomentum[2]
    };
    // Tilt the velocity of all fragments in a way that the momentum error is reduced, 
    // but the velocity direction is not changed significantly. This is a heuristic approach and may not be perfect, 
    // but it should reduce the momentum error significantly.
    for (size_t j = 0; j < _output.size(); ++j){
        double mass = _output.mass[j];
        auto &v = _output.velocity[j];
        double vel_magnitude = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

        // Calculate how much this particle 'contributes' to the error
        // 0.5 is a damping factor (tilting factor) to avoid over-correction
        double factor = 0.5 / (_output.velocity.size() * mass);

        v[0] -= factor * error[0] / (vel_magnitude + 1e-10);
        v[1] -= factor * error[1] / (vel_magnitude + 1e-10);
        v[2] -= factor * error[2] / (vel_magnitude + 1e-10);

        double newRescale = vel_magnitude / std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] *= newRescale;
        v[1] *= newRescale;
        v[2] *= newRescale;
    }

    currentP = calculateCurrentMomentum();
    error = {
        currentP[0] - _initialMomentum[0],
        currentP[1] - _initialMomentum[1],
        currentP[2] - _initialMomentum[2]
    };

}

void Breakup::enforceKineticEnergyConservation() {
    
    double outputKineticEnergy = std::transform_reduce(
        std::execution::par_unseq,
        _output.mass.begin(), 
        _output.mass.end(),
        _output.velocity.begin(),
        0.0,
        std::plus<>(),
        [](double m, const std::array<double, 3>& v) {
            return util::calculateKineticEnergy(m, v);
        }
    );

    if (outputKineticEnergy > 1e-10) {
        double scalingFactor = std::sqrt(_initialKineticEnergy / outputKineticEnergy);
        
        std::for_each(std::execution::par_unseq, _output.velocity.begin(), _output.velocity.end(),
            [scalingFactor](std::array<double, 3>& v) {
                v[0] *= scalingFactor;
                v[1] *= scalingFactor;
                v[2] *= scalingFactor;
            }
        );

    } else {
        spdlog::error("Kinetic energy conservation failed: Output fragments have no velocity/mass.");
    }
}
