#ifndef CREATEINITIALCONDITIONS_H
#define CREATEINITIALCONDITIONS_H


Eigen::MatrixXd CreateInitialConditionsSE(int orbitType, int sampleSize, long double energyLowerBound, long double energyHigherBound);

Eigen::MatrixXd CreateInitialConditionsEM(int orbitType, int sampleSize, long double energyLowerBound, long double energyHigherBound);

#endif // CREATEINITIALCONDITIONS_H
