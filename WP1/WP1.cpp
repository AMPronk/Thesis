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
    bool orbitType4 = true;
    bool orbitType5 = false;
    bool orbitType6 = false;

    //Propagation constants used
    const long double PropagationLength = 0.9 * tudat::physical_constants::SIDEREAL_YEAR; //cannot be more than one year!

    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double initStepSizeB = -1.0;
    double initStepSizeF = 1.0;

    double minimumStepSizeB = -1.0E-5;
    double maximumStepSizeB = -500.0;
    double minimumStepSizeF = 1.0E-5;
    double maximumStepSizeF = 500.0;

    double relativeErrorTolerance = -1.0E-20;
    double absoluteErrorTolerance = -1.0E-12;

    //Initial Conditions Creation constants used
    int InitCondSampleSize = 50;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE EPHEMERIS AND ACCELERATIONS       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Start CreateEphemeris."<<std::endl;

    AccelerationMap ModelMapSE;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    NamedBodyMap bodyMap;

    auto Ephemeris = CreateEphemeris(addSun, addEarth, addMoon, addOlfarSE, addOlfarEM, 0.0, 2.0 * tudat::physical_constants::SIDEREAL_YEAR);
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

    bool runtest = 1;

    Eigen::VectorXd SC1InitialState = InitialStateSEL2();

    if(runtest){

        systemInitialState = SC1InitialState;
        double PropagationStart = 0.0;
        double PropagationEnds = 2.0 * tudat::physical_constants::SIDEREAL_YEAR;

        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1;
        integratorSettingsS1 = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                (rungeKuttaVariableStepSize, PropagationStart, initStepSizeF, coefficientSet,
                 minimumStepSizeF, maximumStepSizeF, relativeErrorTolerance, absoluteErrorTolerance );

        boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
                boost::make_shared< TranslationalStatePropagatorSettings< > >
                (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEnds, cowell);

        ///
        /// Create simulation object and propagate dynamics.
        ///
        std::cerr<<"Start L2 dynamics simulator."<<std::endl;
        SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                    bodyMap, integratorSettingsS1, propagatorSettingsS1, true, true, false );
        std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );


        WriteToFile(tempIntegrationResultS1, bodyMap, "L2trajectory", addSun, addEarth, addMoon);

        tempIntegrationResultS1.clear();

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

        InitialConditionsSE = CreateInitialConditionsSE(1,InitCondSampleSize,3.0,3.2);

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

            PropagationStart = tudat::physical_constants::SIDEREAL_YEAR;
            PropagationEnds = tudat::physical_constants::SIDEREAL_YEAR - PropagationLength;


            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

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

        InitialConditionsSE = CreateInitialConditionsSE(2,InitCondSampleSize,3.0,3.2);

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

            PropagationStart = tudat::physical_constants::SIDEREAL_YEAR;
            PropagationEnds = tudat::physical_constants::SIDEREAL_YEAR + PropagationLength;


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

    if(orbitType3)
    {
        std::cerr<<"Start orbit type 3."<<std::endl;

        InitialConditionsSE = CreateInitialConditionsSE(3,0.5*InitCondSampleSize,3.0,3.2);

        for(int setcount = 2; setcount <= InitialConditionsSE.rows() ; setcount++)
        {
            std::cerr<<"Start OT3 propagation " +
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

            PropagationStart = tudat::physical_constants::SIDEREAL_YEAR;
            double PropagationEndsB = tudat::physical_constants::SIDEREAL_YEAR - 0.5 * PropagationLength;
            double PropagationEndsF = tudat::physical_constants::SIDEREAL_YEAR + 0.5 * PropagationLength;

            // START FORWARD INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsF;
            integratorSettingsF = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeF, coefficientSet,
                     minimumStepSizeF, maximumStepSizeF, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsF =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditionsSE(setcount,4);
            double alpha1 = InitialConditionsSE(setcount,5);
            double alpha2 = InitialConditionsSE(setcount,6);
            double beta = InitialConditionsSE(setcount,7);
            double time = InitialConditionsSE(setcount,8);

            std::string filenameF = "orbittype3_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsB, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorB(
                        bodyMap, integratorSettingsB, propagatorSettingsB, true, true, false );

            tempIntegrationResult = dynamicsSimulatorB.getEquationsOfMotionNumericalSolution( );

            std::string filenameB = "orbittype3_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "B_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType4)
    {
        std::cerr<<"Start orbit type 4."<<std::endl;

        InitialConditionsSE = CreateInitialConditionsSE(4,0.5*InitCondSampleSize,3.0,3.2);

        for(int setcount = 2; setcount <= InitialConditionsSE.rows() ; setcount++)
        {
            std::cerr<<"Start OT4 propagation " +
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

            PropagationStart = tudat::physical_constants::SIDEREAL_YEAR;
            double PropagationEndsB = tudat::physical_constants::SIDEREAL_YEAR - 0.5 * PropagationLength;
            double PropagationEndsF = tudat::physical_constants::SIDEREAL_YEAR + 0.5 * PropagationLength;

            // START FORWARD INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsF;
            integratorSettingsF = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeF, coefficientSet,
                     minimumStepSizeF, maximumStepSizeF, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsF =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditionsSE(setcount,4);
            double alpha1 = InitialConditionsSE(setcount,5);
            double alpha2 = InitialConditionsSE(setcount,6);
            double beta = InitialConditionsSE(setcount,7);
            double time = InitialConditionsSE(setcount,8);

            std::string filenameF = "orbittype4_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMapSE, bodiesToPropagate, systemInitialState, PropagationEndsB, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorB(
                        bodyMap, integratorSettingsB, propagatorSettingsB, true, true, false );

            tempIntegrationResult = dynamicsSimulatorB.getEquationsOfMotionNumericalSolution( );

            std::string filenameB = "orbittype4_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "B_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, addSun, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }




}



