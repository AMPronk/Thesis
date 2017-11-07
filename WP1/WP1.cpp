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
    bool addMoon = false;
    bool addOlfarSE = true;
    bool addOlfarEM = false;
    const long double PropagationStart = 0.0; //1.0 * tudat::physical_constants::SIDEREAL_YEAR;
    const long double PropagationEnds = 1.0 * tudat::physical_constants::SIDEREAL_YEAR;
    long double StepSize = 0.5;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE EPHEMERIS AND ACCELERATIONS       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Start CreateEphemeris."<<std::endl;

    AccelerationMap ModelMapSE;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    NamedBodyMap bodyMap;

    auto Ephemeris = CreateEphemeris(addSun, addEarth, addMoon, addOlfarSE, addOlfarEM, PropagationStart, PropagationEnds);
    ModelMapSE = std::get<0>(Ephemeris);
    bodiesToPropagate = std::get<1>(Ephemeris);
    centralBodies = std::get<2>(Ephemeris);
    bodyMap = std::get<3>(Ephemeris);


    /// Create initial states of all bodies exept Olfar

    int counter = 0;
    if(addOlfarSE){
        counter = counter + 1;
    }
    if(addOlfarEM){
        counter = counter + 1;
    }

    Eigen::VectorXd systemInitialState = Eigen::VectorXd(counter * 6);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     TESTPARTTESTPARTTESTPARTTESTPART       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO REMOVE

    bool runtest = 1;

    Eigen::VectorXd SC1InitialState = InitialStateSEL2();

    Eigen::VectorXd SC2InitialState = InitialStateEML2();

    if(runtest){

        if(addOlfarSE){
            systemInitialState.segment(0 , 6) = SC1InitialState;
            if(addOlfarEM){
                systemInitialState.segment(6 , 6) = SC2InitialState;
            }
        }
        else if(addOlfarEM){
            systemInitialState.segment(0 , 6) = SC2InitialState;
        }

        const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
                numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
        double StepSizeS1 = 1;
        double minimumStepSizeS1 = 1E-5;
        double maximumStepSizeS1 = 1E2;
        double relativeErrorToleranceS1 = 1E-20;
        double absoluteErrorToleranceS1 = 1E-12;

        double PropagationEndsS1 = 2.0 * tudat::physical_constants::SIDEREAL_DAY;

//        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1;
//        integratorSettingsS1 = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
//                (rungeKuttaVariableStepSize, PropagationStart, StepSizeS1, coefficientSetS1,
//                 minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );

        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1 =
                boost::make_shared< IntegratorSettings< > >
                (rungeKutta4, PropagationStart, 1 );

        boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
                boost::make_shared< TranslationalStatePropagatorSettings< > >
                (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsS1, cowell);

        ///
        /// Create simulation object and propagate dynamics.
        ///
        std::cerr<<"Start L2 dynamics simulator."<<std::endl;
        SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                    bodyMap, integratorSettingsS1, propagatorSettingsS1, true, true, false );
        std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );


        std::cerr<<"Write L2-propagation to files."<<std::endl;
        //Write SC data to file
        input_output::writeDataMapToTextFile( tempIntegrationResultS1,
                                              "WP1_L2locations_Spacecraft.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        //Retrieve planet data
        std::map< double, Eigen::VectorXd > sunBarycentricStates;
        std::map< double, Eigen::VectorXd > earthBarycentricStates;
        std::map< double, Eigen::VectorXd > moonBarycentricStates;

        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = tempIntegrationResultS1.begin( );
             stateIterator != tempIntegrationResultS1.end( ); stateIterator++ )
        {
            if(addSun){
                sunBarycentricStates[ stateIterator->first ] = bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris(
                            stateIterator->first );
            }
            if(addEarth){
                earthBarycentricStates[ stateIterator->first ] = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris(
                            stateIterator->first );
            }
            if(addMoon){
                moonBarycentricStates[ stateIterator->first ] = bodyMap.at( "Moon" )->getStateInBaseFrameFromEphemeris(
                            stateIterator->first );
            }
        }

        //write planet data to file
        if(addSun){
            input_output::writeDataMapToTextFile( sunBarycentricStates,
                                                  "WP1_L2locations_Sun.dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
        }
        if(addEarth){
            input_output::writeDataMapToTextFile( earthBarycentricStates,
                                                  "WP1_L2locations_Earth.dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
        }
        if(addMoon){
            input_output::writeDataMapToTextFile( moonBarycentricStates,
                                                  "WP1_L2locations_Moon.dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
        }


        tempIntegrationResultS1.clear();
        sunBarycentricStates.clear();
        earthBarycentricStates.clear();
        moonBarycentricStates.clear();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE SET OF INITIAL CONDITIONS         //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cerr<<"Start initial conditions creation."<<std::endl;
    Eigen::MatrixXd InitialConditionsSE = CreateInitialConditionsSE(2,50,3.0,3.2);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE INTEGRATION AND PROPAGATION       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //bool run, can be switched off for testing.
    bool run = 0;

    if(run){
        std::cerr<<"Start propagation."<<std::endl;

        std::map< double, Eigen::VectorXd > sunBarycentricStates;
        std::map< double, Eigen::VectorXd > earthBarycentricStates;
        std::map< double, Eigen::VectorXd > moonBarycentricStates;

        for(int setcount = 2; setcount <= InitialConditionsSE.rows() ; setcount++)
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
                systemInitialState.segment(0 , 6) = SC1InitialState;
                if(addOlfarEM){
                    systemInitialState.segment(6 , 6) = SC2InitialState;
                }
            }
            else if(addOlfarEM){
                systemInitialState.segment(0 , 6) = SC2InitialState;
            }

            std::cerr<<"Start create integ en propag."<<std::endl;
            ///
            /// Create integrator and propagator
            ///
            const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet =
                    numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
            double minimumStepSize = 1.0E-5;
            double maximumStepSize = 1.0E2;
            double relativeErrorTolerance = 1.0E-20;
            double absoluteErrorTolerance = 1.0E-12;

//            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
//            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
//                    (rungeKuttaVariableStepSize, PropagationStart, StepSize, coefficientSet,
//                     minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                    boost::make_shared< IntegratorSettings< > >
                    (rungeKutta4, PropagationStart, -100 );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEnds, cowell );

            std::cerr<<"Start dynamics simulator."<<std::endl;
            ///
            /// Create simulation object and propagate dynamics.
            ///
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            std::cerr<<"Read planetary states."<<std::endl;

            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = tempIntegrationResult.begin( );
                 stateIterator != tempIntegrationResult.end( ); stateIterator++ )
            {
                if(addSun){
                    sunBarycentricStates[ stateIterator->first ] = bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris(
                                stateIterator->first );
                }
                if(addEarth){
                    earthBarycentricStates[ stateIterator->first ] = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris(
                                stateIterator->first );
                }
                if(addMoon){
                    moonBarycentricStates[ stateIterator->first ] = bodyMap.at( "Moon" )->getStateInBaseFrameFromEphemeris(
                                stateIterator->first );
                }
            }


            ///
            /// Write to file
            ///
            std::cerr<<"Started writing to file."<<std::endl;

            //write planet data to file
            if(addSun){
                input_output::writeDataMapToTextFile( sunBarycentricStates,
                                                      "WP1_SunOrbit" + boost::lexical_cast< std::string >( setcount )
                                                      + "_JacobiConstant" + boost::lexical_cast< std::string >( InitialConditionsSE(setcount,4) )
                                                      + ".dat",
                                                      tudat_applications::getOutputPath( ),
                                                      "",
                                                      std::numeric_limits< double >::digits10,
                                                      std::numeric_limits< double >::digits10,
                                                      "," );
            }
            if(addEarth){
                input_output::writeDataMapToTextFile( earthBarycentricStates,
                                                      "WP1_EarthOrbit" + boost::lexical_cast< std::string >( setcount )
                                                      + "_JacobiConstant" + boost::lexical_cast< std::string >( InitialConditionsSE(setcount,4) )
                                                      + ".dat",
                                                      tudat_applications::getOutputPath( ),
                                                      "",
                                                      std::numeric_limits< double >::digits10,
                                                      std::numeric_limits< double >::digits10,
                                                      "," );
            }
            if(addMoon){
                input_output::writeDataMapToTextFile( moonBarycentricStates,
                                                      "WP1_MoonOrbit" + boost::lexical_cast< std::string >( setcount )
                                                      + "_JacobiConstant" + boost::lexical_cast< std::string >( InitialConditionsSE(setcount,4) )
                                                      + ".dat",
                                                      tudat_applications::getOutputPath( ),
                                                      "",
                                                      std::numeric_limits< double >::digits10,
                                                      std::numeric_limits< double >::digits10,
                                                      "," );
            }

            //Write SC data to file
            input_output::writeDataMapToTextFile( tempIntegrationResult,
                                                  "WP1_Orbit" + boost::lexical_cast< std::string >( setcount )
                                                  + "_JacobiConstant" + boost::lexical_cast< std::string >( InitialConditionsSE(setcount,4) )
                                                  + ".dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            tempIntegrationResult.clear();
            sunBarycentricStates.clear();
            earthBarycentricStates.clear();
            moonBarycentricStates.clear();
        }

    }
    else{
        std::cerr<<"Do not start propagation."<<std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////     WRITE TO FILE       /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}



