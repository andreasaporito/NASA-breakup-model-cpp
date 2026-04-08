#pragma once

#include "Breakup.h"

/**
 * A collision Breakup of two satellites.
 * @attention There is no check if the two satellites are actually at the same position!
 */
class Collision : public Breakup {

protected:
    bool _isCatastrophic;

public:

    using Breakup::Breakup;
    void run() override;

private:

    void init() override;

    void calculateFragmentCount() override;

    void assignParentProperties() override;

protected:
    void addFurtherFragments() override;

public:

    bool isIsCatastrophic() const {
        return _isCatastrophic;
    }

};

class SubCollision : public Collision {
    int _cardinality;

public: 
    SubCollision(std::vector<Satellite> input,
                 int cardinality,
                 double minimalCharacteristicLength,
                 size_t currentMaxGivenID,
                 bool enforceMassConservation);
    // Bypasses Collision::run() -> runs the standard Breakup pipeline only
    void run() override;
    int getCardinality() const { return _cardinality; }

private:
    void calculateFragmentCount() override;
    void assignParentProperties() override;
    void addFurtherFragments() override;
};
