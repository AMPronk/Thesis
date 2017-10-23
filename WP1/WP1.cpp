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
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/initiateconstants.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/initiateconstants.cpp"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/CreateEphemeris.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/CreateEphemeris.cpp"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/FrameTransformation.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/FrameTransformation.cpp"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/createinitialconditions.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/createinitialconditions.cpp"

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
    const long double PropagationStartS1 = 0.0;
    const long double PropagationEndsS1 = 5 * tudat::physical_constants::JULIAN_YEAR;
    long double StepSize = 100.0;


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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     TESTPARTTESTPARTTESTPARTTESTPART       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO REMOVE

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


    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double minimumStepSizeS1 = 1.0E-5;
    double maximumStepSizeS1 = 1000;
    double relativeErrorToleranceS1 = 1.0E-12;
    double absoluteErrorToleranceS1 = 1.0E-12;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1;
    integratorSettingsS1 = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
    (rungeKuttaVariableStepSize, PropagationStartS1, StepSize, coefficientSetS1,
    minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
            boost::make_shared< TranslationalStatePropagatorSettings< > >
            (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsS1, cowell );

    ///
    /// Create simulation object and propagate dynamics.
    ///
    std::cerr<<"Start L2 dynamics simulator."<<std::endl;
    SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                bodyMap, integratorSettingsS1, propagatorSettingsS1, true, false, false );
    std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );


    std::cerr<<"Write L2 to file."<<std::endl;
    input_output::writeDataMapToTextFile( tempIntegrationResultS1,
                                          "WP1_L2locations.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE SET OF INITIAL CONDITIONS         //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cerr<<"Start initial conditions creation."<<std::endl;
    Eigen::MatrixXd InitialConditionsSE = CreateInitialConditionsSE(0,1.0E-4,3.0,3.2);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE INTEGRATION AND PROPAGATION       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //bool run, can be switched of for testing.
    bool run = 1;

    if(run){
        std::cerr<<"Start propagation."<<std::endl;
        for(int setcount = 0; setcount <= InitialConditionsSE.rows() ; setcount++)
        {
            std::cerr<<"Start propagation " +
                      boost::lexical_cast< std::string >( setcount )
                      + " out of " +
                      boost::lexical_cast< std::string >( InitialConditionsSE.rows() )
                      + "."<<std::endl;

            SC1InitialState = FrameTransformationSE(InitialConditionsSE(setcount,0),
                                                    InitialConditionsSE(setcount,1),
                                                    InitialConditionsSE(setcount,2),
                                                    InitialConditionsSE(setcount,3));


            if(addOlfarSE){
                systemInitialState.segment(counter*6 , 6) = SC1InitialState;
                if(addOlfarEM){
                    systemInitialState.segment(counter*6+6 , 6) = SC2InitialState;
                }
            }
            else if(addOlfarEM){
                systemInitialState.segment(counter*6 , 6) = SC2InitialState;
            }

            std::cerr<<"Start create integ en propag."<<std::endl;
            ///
            /// Create integrator and propagator
            ///
            const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet =
                    numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
            double minimumStepSize = 1.0E-3;
            double maximumStepSize = 1000;
            double relativeErrorTolerance = 1.0E-12;
            double absoluteErrorTolerance = 1.0E-12;

            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            (rungeKuttaVariableStepSize, PropagationStartS1, StepSize, coefficientSet,
            minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsS1, cowell );

            std::cerr<<"Start dynamics simulator."<<std::endl;
            ///
            /// Create simulation object and propagate dynamics.
            ///
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );
            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


            ///
            /// Write to file
            ///
            std::cerr<<"Started writing to file."<<std::endl;

            input_output::writeDataMapToTextFile( tempIntegrationResult,
                                                  "WP1_Orbit" + boost::lexical_cast< std::string >( setcount )
                                                  + "_JacobiConstant" + boost::lexical_cast< std::string >( InitialConditionsSE(setcount,4) )
                                                  + ".dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
        }

    }
    else{
        std::cerr<<"Do not start propagation."<<std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////     WRITE TO FILE       /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}



