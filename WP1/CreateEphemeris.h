//CreateEphemeris.h

#ifndef CREATEEPHEMERIS_H
#define CREATEEPHEMERIS_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::basic_astrodynamics;

std::tuple<AccelerationMap, std::vector<std::string>, std::vector<std::string>, NamedBodyMap> CreateEphemeris(bool sun, bool earth, bool moon, bool SCSE, bool SCEM);

Eigen::VectorXd CreateEphemerisInitialState(bool sun, bool earth, bool moon);

Eigen::VectorXd InitialStateSEL2();

Eigen::VectorXd InitialStateEML2();

#endif
