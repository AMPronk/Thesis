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
    bool addOlfarSE = addSun;
    bool addOlfarEM = addMoon; //Select EITHER addOlfarSE OR addOlfarEM -- NEVER BOTH.

    bool runL2 = true;

    bool orbitType1 = false;
    bool orbitType2 = false;
    bool orbitType3 = false;
    bool orbitType4 = false;
    bool orbitType5 = false;
    bool orbitType6 = false;

    //Propagation constants used
    double InitialStateTime = tudat::physical_constants::SIDEREAL_YEAR;
    long double PropagationLength;

    if(addOlfarSE){
        PropagationLength = 1.0 * tudat::physical_constants::SIDEREAL_YEAR; //cannot be more than one year!
    }
    else if(addOlfarEM){
        PropagationLength = 0.3 * tudat::physical_constants::SIDEREAL_YEAR; //cannot be more than one year!
    }

    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double initStepSizeB = -1.0;
    double initStepSizeF = 1.0;

    double minimumStepSizeB;
    double maximumStepSizeB;
    double minimumStepSizeF;
    double maximumStepSizeF;

    if(addOlfarSE){
        minimumStepSizeB = -1.0;
        maximumStepSizeB = -500.0;
        minimumStepSizeF = 1.0;
        maximumStepSizeF = 500.0;
    }
    else if(addOlfarEM){
        minimumStepSizeB = -0.1;
        maximumStepSizeB = -50.0;
        minimumStepSizeF = 0.1;
        maximumStepSizeF = 50.0;
    }

    double relativeErrorTolerance = -1.0E-20;
    double absoluteErrorTolerance = -1.0E-12;

    //Initial Conditions Creation constants used
    int InitCondSampleSize = 17;

    double ClowSE = 3.0;
    double ChighSE = 3.002;

    double ClowEM = 3.19;//3.195;
    double ChighEM = 3.2;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE EPHEMERIS AND ACCELERATIONS       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Start CreateEphemeris."<<std::endl;

    AccelerationMap ModelMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    NamedBodyMap bodyMap;

    auto Ephemeris = CreateEphemeris(addSun, addEarth, addMoon, addOlfarSE, addOlfarEM, InitialStateTime, PropagationLength);
    ModelMap = std::get<0>(Ephemeris);
    bodiesToPropagate = std::get<1>(Ephemeris);
    centralBodies = std::get<2>(Ephemeris);
    bodyMap = std::get<3>(Ephemeris);


    /// Create empty initial state vectors for spacecraft

    Eigen::VectorXd SCInitialState = Eigen::VectorXd(6);



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     TESTPARTTESTPARTTESTPARTTESTPART       //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO REMOVE


    if(addOlfarSE)
    {
        std::cerr<<"Start SE L2 calculation."<<std::endl;
        SCInitialState = InitialStateSEL2();
    }
    else
    {
        std::cerr<<"Start EM L2 calculation."<<std::endl;
        SCInitialState = InitialStateEML2();
    }


    if(runL2){

        double PropagationStart = InitialStateTime;
        double PropagationEnds = InitialStateTime + PropagationLength; //tudat::physical_constants::SIDEREAL_YEAR;

        boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1;
        integratorSettingsS1 = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                (rungeKuttaVariableStepSize, PropagationStart, initStepSizeF, coefficientSet,
                 minimumStepSizeF, maximumStepSizeF, relativeErrorTolerance, absoluteErrorTolerance );

        boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
                boost::make_shared< TranslationalStatePropagatorSettings< > >
                (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEnds, cowell);

        ///
        /// Create simulation object and propagate dynamics.
        ///
        SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                    bodyMap, integratorSettingsS1, propagatorSettingsS1, true, true, false );
        std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );


        WriteToFile(tempIntegrationResultS1, bodyMap, "L2trajectory", addSun, addEarth, addMoon);

        tempIntegrationResultS1.clear();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE NECESSARY VARIABLE SAVINGS               ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    Eigen::MatrixXd InitialConditions;
    long double PropagationStart;
    long double PropagationEnds;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     RUN ORBIT TYPE 1: ASYMPTOTIC, TOWARDS L2        ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(orbitType1)
    {
        std::cerr<<"Start orbit type 1."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(1,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(1,InitCondSampleSize,ClowEM,ChighEM);
        }


        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT1 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }


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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEnds, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filename = "orbittype1_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filename, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType2)
    {
        std::cerr<<"Start orbit type 2."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(2,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(2,InitCondSampleSize,ClowEM,ChighEM);
        }

        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT2 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }


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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEnds, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filename = "orbittype2_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filename, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType3)
    {
        std::cerr<<"Start orbit type 3."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(3,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(3,InitCondSampleSize,ClowEM,ChighEM);
        }

        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT3 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }

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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filenameF = "orbittype3_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, false, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsB, cowell );

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
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType4)
    {
        std::cerr<<"Start orbit type 4."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(4,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(4,InitCondSampleSize,ClowEM,ChighEM);
        }

        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT4 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }

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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filenameF = "orbittype4_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, false, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsB, cowell );

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
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType5)
    {
        std::cerr<<"Start orbit type 5."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(5,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(5,InitCondSampleSize,ClowEM,ChighEM);
        }

        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT5 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }

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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filenameF = "orbittype5_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, false, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsB, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorB(
                        bodyMap, integratorSettingsB, propagatorSettingsB, true, true, false );

            tempIntegrationResult = dynamicsSimulatorB.getEquationsOfMotionNumericalSolution( );

            std::string filenameB = "orbittype5_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "B_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }

    if(orbitType6)
    {
        std::cerr<<"Start orbit type 6."<<std::endl;

        if(addOlfarSE)
        {
            InitialConditions = CreateInitialConditionsSE(6,InitCondSampleSize,ClowSE,ChighSE);
        }
        if(addOlfarEM)
        {
            InitialConditions = CreateInitialConditionsEM(6,InitCondSampleSize,ClowEM,ChighEM);
        }

        for(int setcount = 2; setcount <= InitialConditions.rows() ; setcount++)
        {
            std::cerr<<"Start OT6 propagation " +
                       boost::lexical_cast< std::string >( setcount )
                       + " out of " +
                       boost::lexical_cast< std::string >( InitialConditions.rows() )
                       + "."<<std::endl;

            ///
            /// Create initial conditions
            ///

            if(addOlfarSE)
            {
                SCInitialState = FrameTransformationSE(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }
            if(addOlfarEM)
            {
                SCInitialState = FrameTransformationEM(InitialConditions(setcount,0),
                                                       InitialConditions(setcount,1),
                                                       InitialConditions(setcount,2),
                                                       InitialConditions(setcount,3));
            }

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
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsF, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorF(
                        bodyMap, integratorSettingsF, propagatorSettingsF, true, true, false );

            std::map< double, Eigen::VectorXd > tempIntegrationResult = dynamicsSimulatorF.getEquationsOfMotionNumericalSolution( );

            double C = InitialConditions(setcount,4);
            double alpha1 = InitialConditions(setcount,5);
            double alpha2 = InitialConditions(setcount,6);
            double beta = InitialConditions(setcount,7);
            double time = InitialConditions(setcount,8);

            std::string filenameF = "orbittype6_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "F_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameF, false, addEarth, addMoon);

            tempIntegrationResult.clear();

            // START BACKWARDS INTEGRATION
            boost::shared_ptr< IntegratorSettings< > > integratorSettingsB;
            integratorSettingsB = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    (rungeKuttaVariableStepSize, PropagationStart, initStepSizeB, coefficientSet,
                     minimumStepSizeB, maximumStepSizeB, relativeErrorTolerance, absoluteErrorTolerance );

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsB =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, ModelMap, bodiesToPropagate, SCInitialState, PropagationEndsB, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///

            SingleArcDynamicsSimulator< > dynamicsSimulatorB(
                        bodyMap, integratorSettingsB, propagatorSettingsB, true, true, false );

            tempIntegrationResult = dynamicsSimulatorB.getEquationsOfMotionNumericalSolution( );

            std::string filenameB = "orbittype6_nr=" + boost::lexical_cast< std::string >( setcount )
                    + "B_C=" + boost::lexical_cast< std::string >( C )
                    + "_a1=" + boost::lexical_cast< std::string >( alpha1 )
                    + "_a2=" + boost::lexical_cast< std::string >( alpha2 )
                    + "_b=" + boost::lexical_cast< std::string >( beta )
                    + "_t=" + boost::lexical_cast< std::string >( time );

            WriteToFile(tempIntegrationResult, bodyMap, filenameB, false, addEarth, addMoon);

            tempIntegrationResult.clear();

        }
    }





}



