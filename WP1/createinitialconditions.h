#ifndef CREATEINITIALCONDITIONS_H
#define CREATEINITIALCONDITIONS_H


Eigen::MatrixXd CreateInitialConditionsSE(int orbitType, long double precision, long double energyLowerBound, long double energyHigherBound);

Eigen::MatrixXd CreateInitialConditionsEM(int orbitType, long double precision, long double energyLowerBound, long double energyHigherBound);

#endif // CREATEINITIALCONDITIONS_H
