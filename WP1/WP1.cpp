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
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/writetofile.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/writetofile.cpp"



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
    bool addOlfarEM = false; //TODO, script doesn't support this right now.

    bool orbitType1 = true;
    bool orbitType2 = true;
    bool orbitType3 = true;
    bool orbitType4 = false;
    bool orbitType5 = false;
    bool orbitType6 = false;

    const long double PropagationLength = 1.0 * tudat::physical_constants::SIDEREAL_YEAR;

    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double initStepSizeB = -1.0;
    double initStepSizeF = 1.0;

    double minimumStepSizeB = -1.0E-5;
    double maximumStepSizeB = -50.0;
    double minimumStepSizeF = 1.0E-5;
    double maximumStepSizeF = 50.0;

    double relativeErrorTolerance = -1.0E-20;
    double absoluteErrorTolerance = -1.0E-12;

    int StepSize = 50;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE EPHEMERIS AND ACCELERATIONS       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Start CreateEphemeris."<<std::endl;

    AccelerationMap ModelMapSE;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    NamedBodyMap bodyMap;

    auto Ephemeris = CreateEphemeris(addSun, addEarth, addMoon, addOlfarSE, addOlfarEM, 0.0, PropagationLength);
    ModelMapSE = std::get<0>(Ephemeris);
    bodiesToPropagate = std::get<1>(Ephemeris);
    centralBodies = std::get<2>(Ephemeris);
    bodyMap = std::get<3>(Ephemeris);


    /// Create initial states of all bodies exept Olfar

    Eigen::VectorXd systemInitialState = Eigen::VectorXd(6);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     TESTPARTTESTPARTTESTPARTTESTPART       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO REMOVE

    bool runtest = 0;

    Eigen::VectorXd SC1InitialState = InitialStateSEL2();

    if(runtest){

        systemInitialState = SC1InitialState;

        const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
                numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
        double StepSizeS1 = -20.0;
        double minimumStepSizeS1 = -1.0;//E-5;
        double maximumStepSizeS1 = -50.0;
        double relativeErrorToleranceS1 = 1.0E-20;
        double absoluteErrorToleranceS1 = 1.0E-12;

        double PropagationStartS1 = PropagationLength;
        double PropagationEndsS1 = 0.0;

//        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1;
//        integratorSettingsS1 = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
//                (rungeKuttaVariableStepSize, PropagationStart, StepSizeS1, coefficientSetS1,
//                 minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );

        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1 =
                    boost::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, PropagationStartS1, -50 );

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
    ///////////////////////     CREATE NECESSARY VARIABLE SAVINGS               ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    Eigen::MatrixXd InitialConditionsSE;
    Eigen::MatrixXd InitialConditionsEM;
    long double PropagationStart;
    long double PropagationEnds;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     RUN ORBIT TYPE 1: ASYMPTOTIC, TOWARDS L2        ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(orbitType1)
    {
        std::cerr<<"Start orbit type 1."<<std::endl;

        InitialConditionsSE = CreateInitialConditionsSE(1,StepSize,3.0,3.2);

        for(int setcount = 2; setcount <= InitialConditionsSE.rows() ; setcount++)
        {
            std::cerr<<"Start OT1 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditionsSE.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            SC1InitialState = FrameTransformationSE(InitialConditionsSE(setcount,0),
                                                    InitialConditionsSE(setcount,1),
                                                    InitialConditionsSE(setcount,2),
                                                    InitialConditionsSE(setcount,3));

            systemInitialState = SC1InitialState;



            ///
            /// Create Integrator and Propagator
            ///

            PropagationStart = PropagationLength;
            PropagationEnds = 0.0;


            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

//            boost::shared_ptr< IntegratorSettings< > > integratorSettings =
//                        boost::make_shared< IntegratorSettings< > >
//                        ( rungeKutta4, PropagationStart, -50 );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEnds, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditionsSE(setcount,4);
            double alpha1 = InitialConditionsSE(setcount,5);
            double alpha2 = InitialConditionsSE(setcount,6);
            double beta = InitialConditionsSE(setcount,7);
            double time = InitialConditionsSE(setcount,8);

            std::string filename = "orbittype1_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filename, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType2)
    {
        std::cerr<<"Start orbit type 2."<<std::endl;

        InitialConditionsSE = CreateInitialConditionsSE(2,StepSize,3.0,3.2);

        for(int setcount = 2; setcount <= InitialConditionsSE.rows() ; setcount++)
        {
            std::cerr<<"Start OT2 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditionsSE.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            SC1InitialState = FrameTransformationSE(InitialConditionsSE(setcount,0),
                                                    InitialConditionsSE(setcount,1),
                                                    InitialConditionsSE(setcount,2),
                                                    InitialConditionsSE(setcount,3));

            systemInitialState = SC1InitialState;



            ///
            /// Create Integrator and Propagator
            ///

            PropagationStart = 0.0;
            PropagationEnds = PropagationLength;


            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeF, coefficientSet,
                    minimumStepSizeF, maximumStepSizeF, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEnds, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditionsSE(setcount,4);
            double alpha1 = InitialConditionsSE(setcount,5);
            double alpha2 = InitialConditionsSE(setcount,6);
            double beta = InitialConditionsSE(setcount,7);
            double time = InitialConditionsSE(setcount,8);

            std::string filename = "orbittype2_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filename, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }





}



