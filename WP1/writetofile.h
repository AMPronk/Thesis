#ifndef WRITETOFILE_H
#define WRITETOFILE_H


void WriteToFile(std::map<double, Eigen::VectorXd> integrationResult, tudat::simulation_setup::NamedBodyMap bodyMap, std::string IDname, bool addSun, bool addEarth, bool addMoon);

#endif // WRITETOFILE_H
