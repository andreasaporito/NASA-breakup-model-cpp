#include "Breakup.h"

void Breakup::run() {
    //0. Step: Prepare constants, etc.
    this->init();

    //1. Step: Generate the new Satellites
    this->calculateFragmentCount();

    for (auto & subBatch : _subOutputs) {

        //2. Step: Assign every new Satellite a value for L_c
        this->characteristicLengthDistribution(subBatch.fragments);
        
        //3. Step: Calculate the A/M (area-to-mass-ratio), A (area) and M (mass) values for every Satellite
        this->areaToMassRatioDistribution(subBatch.fragments);
        //4. Step: Enforce the Mass Conservation and remove (or add) fragments
        this->enforceMassConservation(subBatch);

        //5. Step: Assign parent and by doing that assign each fragment a base velocity
        this->assignParentProperties(subBatch);

        //6. Step: Calculate the Ejection velocity for every Satellite
        this->deltaVelocityDistribution(subBatch.fragments);
    }

    // 7. Merge outputs
    this->mergeSubOutputs();
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

void Breakup::generateFragments(size_t fragmentCount, double targetMass, const std::array<double,3> &parentVelocity, const std::string namePtr, const std::array<double, 3> &position) {
    
    SubCollision sub;
    sub.fragments = Satellites{_currentMaxGivenID + 1, SatType::DEBRIS, position, fragmentCount};
    sub.targetInputMass = targetMass;
    sub.currentOutputMass = 0.0;
    sub.parentNamePtr = std::make_shared<const std::string>(namePtr + "-Fragment");
    sub.parentVelocity = parentVelocity;
    
    // Append a new Satellites batch to our temporary vector
    _subOutputs.push_back(std::move(sub));
    // Update the ID so the next sub-collision gets correct IDs
    _currentMaxGivenID += fragmentCount;
}

void Breakup::characteristicLengthDistribution(Satellites &targetOutput) {
    std::for_each(std::execution::par_unseq, targetOutput.characteristicLength.begin(), targetOutput.characteristicLength.end(),
                  [&](double &lc) {
        lc = calculateCharacteristicLength();
    });
}

void Breakup::areaToMassRatioDistribution(Satellites &targetOutput) {
    auto tupleView = targetOutput.getAreaMassTuple();
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

void Breakup::enforceMassConservation(SubCollision& subBatch) {
    //Enforce Mass Conservation if the output mass is greater than the input mass
    subBatch.currentOutputMass = std::reduce(std::execution::par_unseq, 
                                             subBatch.fragments.mass.begin(), 
                                             subBatch.fragments.mass.end(), 0.0);
    spdlog::debug("Batch target mass: {} kg", subBatch.targetInputMass);
    spdlog::debug("Batch produced mass: {} kg", subBatch.currentOutputMass);
    size_t oldSize = subBatch.fragments.size();
    size_t newSize = subBatch.fragments.size();
    // Shrink and Remove Mass Excess
    if (subBatch.currentOutputMass > subBatch.targetInputMass) {
        size_t newSize = oldSize;
        while (newSize > 0 && subBatch.currentOutputMass > subBatch.targetInputMass) {
            subBatch.currentOutputMass -= subBatch.fragments.mass[newSize - 1];
            newSize--;
        }
        subBatch.fragments.resize(newSize);
    }

    if (_enforceMassConservation && subBatch.fragments.size() == oldSize) {
        this->addFurtherFragments(subBatch);
    }

    if (oldSize != subBatch.fragments.size()) {
        spdlog::warn("Fragment count adapted from {} to {}.", oldSize, subBatch.fragments.size());
        // FIX: Use local mass for debug
        spdlog::debug("Corrected batch mass: {} kg", subBatch.currentOutputMass);
    }
}

void Breakup::addFurtherFragments(SubCollision& sub) {
    while (sub.currentOutputMass < sub.targetInputMass) {
        //Order in the tuple: 0: L_c | 1: A/M | 2: Area | 3: Mass
        //Create new element and assign values
        auto tuple = sub.fragments.appendElement();
        auto &[lc, areaToMassRatio, area, mass] = tuple;
        lc = calculateCharacteristicLength();
        areaToMassRatio = calculateAreaMassRatio(lc);
        area = calculateArea(lc);
        mass = calculateMass(area, areaToMassRatio);

        //Calculate new mass
        sub.currentOutputMass += mass;
    }
    //Remove the element which has lead to the exceeding of the mass budget
    sub.currentOutputMass -= sub.fragments.mass.back();
    sub.fragments.popBack();
}

void Breakup::deltaVelocityDistribution(Satellites &targetOutput) {
    using namespace util;
    auto tupleView = targetOutput.getVelocityTuple();
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
    std::uniform_real_distribution<> uniformRealDistribution{0.0, 1.0};

    double u = getRandomNumber(uniformRealDistribution) * 2.0 - 1.0;
    double theta = getRandomNumber(uniformRealDistribution) * 2.0 * util::PI;
    double v = std::sqrt(1.0 - u * u);

    return std::array<double, 3>
            {{v * std::cos(theta) * velocity, v * std::sin(theta) * velocity, u * velocity}};
}


void Breakup::assignParentProperties(SubCollision& sub) {
    auto tupleView = sub.fragments.getVNTuple(); // Just velocity and name needed
    std::for_each(std::execution::par_unseq, tupleView.begin(), tupleView.end(),
                  [&](auto &tuple) {
                      auto &[velocity, name] = tuple;
                      velocity = sub.parentVelocity;
                      name = sub.parentNamePtr;
                  });
}

void Breakup::mergeSubOutputs() {
    size_t totalFragments = 0;
    for (const auto& sub : _subOutputs) totalFragments += sub.fragments.size();
    
    // Resize final _output once for performance
    _output.resize(totalFragments); 

    size_t offset = 0;
    for (auto& sub : _subOutputs) {
        // You'll need a method in your Satellites class to copy a range 
        // into another Satellites object at a specific offset.
        _output.copyFrom(sub.fragments, offset);
        offset += sub.fragments.size();
    }
}