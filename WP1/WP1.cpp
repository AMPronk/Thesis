/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>
#include <tuple>
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/CreateEphemeris.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/CreateEphemeris.cpp"

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getOutputPath( )
{
    return "D:/FILES/Documents/Thesis/PropagationResults";
}

}

int main( )
{

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

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0009.tls");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////     SETTING OF NECESSARY VARIABLES       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool addSun = true;
    bool addEarth = true;
    bool addMoon = true;
    bool addOlfarSE = true;
    bool addOlfarEM = true;
    const double PropagationStartS1 = 0.0;
    const double PropagationEndsS1 = 10 * tudat::physical_constants::JULIAN_YEAR;
    double StepSize = 2000.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE EPHEMERIS AND ACCELERATIONS       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Start CreateEphemeris."<<std::endl;

    int nrOfBodies = addSun + addEarth + addMoon + addOlfarSE + addOlfarEM;

    AccelerationMap ModelMapSE;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    NamedBodyMap bodyMap;

    auto Ephemeris = CreateEphemeris(addSun, addEarth, addMoon, addOlfarSE, addOlfarEM);
    ModelMapSE = std::get<0>(Ephemeris);
    bodiesToPropagate = std::get<1>(Ephemeris);
    centralBodies = std::get<2>(Ephemeris);
    bodyMap = std::get<3>(Ephemeris);


    /// Create initial states of all bodies exept Olfar
    Eigen::VectorXd systemInitialState = Eigen::VectorXd(nrOfBodies * 6);
    int counter = nrOfBodies;
    if(addOlfarSE){
        counter = counter - 1;
    }
    if(addOlfarEM){
        counter = counter - 1;
    }

    systemInitialState.segment(0 , 6 * counter) = CreateEphemerisInitialState(addSun, addEarth, addMoon);

    /// Create initial states for Olfar
    // TODO

    Eigen::VectorXd SC1InitialState = InitialStateSEL2();

    Eigen::VectorXd SC2InitialState = InitialStateEML2();

    if(addOlfarSE){
        systemInitialState.segment(counter*6 , 6) = SC1InitialState;
        if(addOlfarEM){
            systemInitialState.segment(counter*6+6 , 6) = SC2InitialState;
        }
    }
    else if(addOlfarEM){
        systemInitialState.segment(counter*6 , 6) = SC2InitialState;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE INTEGRATION AND PROPAGATION       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1 =
            boost::make_shared< IntegratorSettings< > >
            (rungeKutta4, PropagationStartS1, StepSize );

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
            boost::make_shared< TranslationalStatePropagatorSettings< > >
            (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsS1, cowell );

    ///
    /// Create simulation object and propagate dynamics.
    ///
    std::cerr<<"Start dynamics simulator."<<std::endl;
    SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                bodyMap, integratorSettingsS1, propagatorSettingsS1, true, false, false );
    std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////     WRITE TO FILE       /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cerr<<"Write to file."<<std::endl;
    input_output::writeDataMapToTextFile( tempIntegrationResultS1,
                                          "WP1_TESTRUN.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
}


