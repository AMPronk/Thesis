//CreateEphemeris.cpp

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>
#include <tuple>
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"

#include "functions/librationPointLocationFunction.h"
#include "functions/librationPointLocationFunction1.h"
#include "functions/librationPointLocationFunction2.h"
#include "C:/tudatBundle/tudatApplications/Thesis/WP1/CreateEphemeris.h"

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

std::tuple<AccelerationMap, std::vector< std::string >, std::vector< std::string >, NamedBodyMap> CreateEphemeris(bool sun, bool earth, bool moon, bool SCSE, bool SCEM, double InitTime, double PropagationLength)
{
    ///
    /// The Create Ephemeris script creates the ephemeris model based on the input variables.
    /// Booleans sun, earth, moon, SCSE and SCEM determine which bodies are used in the model.
    /// Doubles InitTime and PropagationLength determine the timespan for which the ephemeris model is created.
    ///

    // Create empty bodymap and bodysettings variables to fill up.
    NamedBodyMap bodyMap;
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;

    // Determine which planetary bodies should be created in the ephemeris model, based on the input booleans.
    std::vector< std::string > bodiesToCreate;

    if(sun)
    {
        bodiesToCreate.push_back( "Sun" );
    }
    if(earth)
    {
        bodiesToCreate.push_back( "Earth" );
    }
    if(moon)
    {
        bodiesToCreate.push_back( "Moon" );
    }

    // Retrieve the standard settings for each of the bodies, over a timespan determined by the input doubles
    bodySettings = getDefaultBodySettings( bodiesToCreate, InitTime - PropagationLength - 3000.0, InitTime + PropagationLength + 3000.0 );


    // Set the initial conditions of the planetary bodies
    if(sun)
    {
        // The Sun remains at the center of the ephemeris model
        Eigen::Vector6d SunInitialState;
        SunInitialState( 0 ) = 0.0;
        SunInitialState( 1 ) = 0.0;
        SunInitialState( 2 ) = 0.0;
        SunInitialState( 3 ) = 0.0;
        SunInitialState( 4 ) = 0.0;
        SunInitialState( 5  ) = 0.0;

        bodySettings [ "Sun" ]-> ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    SunInitialState, "SSB", "J2000" );
        bodySettings [ "Sun" ]-> gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( SunGravitationalParameter );

        bodySettings [ "Sun" ]-> atmosphereSettings = NULL;
        bodySettings [ "Sun" ]-> rotationModelSettings = NULL;
        bodySettings [ "Sun" ]-> shapeModelSettings = NULL;
    }

    if(earth)
    {
        // The Earth has a circular orbit around the Sun
        Eigen::Vector6d EarthInitialStateKeplerian;
        EarthInitialStateKeplerian( semiMajorAxisIndex) = AU;
        EarthInitialStateKeplerian( eccentricityIndex) = 0.0;
        EarthInitialStateKeplerian( inclinationIndex) = 0.0;
        EarthInitialStateKeplerian( argumentOfPeriapsisIndex) = 0.0;
        EarthInitialStateKeplerian( longitudeOfAscendingNodeIndex) = 0.0;
        EarthInitialStateKeplerian( trueAnomalyIndex ) = 0.0;

        bodySettings [ "Earth" ]-> ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                    EarthInitialStateKeplerian, InitTime, SunGravitationalParameter, "SSB", "J2000" );
        bodySettings [ "Earth" ]-> gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( EarthGravitationalParameter );
        bodySettings [ "Earth" ]-> atmosphereSettings = NULL;
        bodySettings [ "Earth" ]-> rotationModelSettings = NULL;
        bodySettings [ "Earth" ]-> shapeModelSettings = NULL;
    }

    if(moon)
    {
        // The Moon has a circular orbit around the Earth
        Eigen::Vector6d MoonInitialStateKeplerian;
        MoonInitialStateKeplerian( semiMajorAxisIndex) = MoonEarthDistance;
        MoonInitialStateKeplerian( eccentricityIndex) = 0.0;
        MoonInitialStateKeplerian( inclinationIndex) = 0.0;
        MoonInitialStateKeplerian( argumentOfPeriapsisIndex) = 0.0;
        MoonInitialStateKeplerian( longitudeOfAscendingNodeIndex) = 0.0;
        MoonInitialStateKeplerian( trueAnomalyIndex ) = 0.0;

        bodySettings [ "Moon" ]-> ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                    MoonInitialStateKeplerian, InitTime, EarthGravitationalParameter, "Earth", "J2000" );
        bodySettings [ "Moon" ]-> gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( MoonGravitationalParameter );
        bodySettings [ "Moon" ]-> atmosphereSettings = NULL;
        bodySettings [ "Moon" ]-> rotationModelSettings = NULL;
        bodySettings [ "Moon" ]-> shapeModelSettings = NULL;
    }

    // Enter the settings implemented above into the ephemeris model
    bodyMap = createBodies(bodySettings);

    if(SCSE)
    {
        // Create the spacecraft point mass that is attracted by the sun and the earth
        bodyMap[ "Spacecraft1" ] = boost::make_shared< simulation_setup::Body >( );
    }

    if(SCEM)
    {
        // Create the spacecraft point mass that is attracted by the earth and the moon
        bodyMap[ "Spacecraft2" ] = boost::make_shared< simulation_setup::Body >( );
    }

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSCSE;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSCEM;

    // Attractions in the Sun Earth system
    if(SCSE){
        if(sun){
            accelerationsOfSCSE[ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        if(earth){
            accelerationsOfSCSE[ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Spacecraft1" ] = accelerationsOfSCSE;
        bodiesToPropagate.push_back( "Spacecraft1" );
        centralBodies.push_back( "Earth" );
    }

    // Attractions in the Earth Moon system
    if(SCEM){
        if(earth){
            accelerationsOfSCEM[ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        if(moon){
            accelerationsOfSCEM[ "Moon" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Spacecraft2" ] = accelerationsOfSCEM;
        bodiesToPropagate.push_back( "Spacecraft2" );
        centralBodies.push_back( "Moon" );
    }

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    std::cerr<<"End CreateEphemeris."<<std::endl;

    return std::make_tuple(accelerationModelMap, bodiesToPropagate, centralBodies, bodyMap);
}


Eigen::VectorXd InitialStateSEL2()
{
    ///
    /// InitialStateSEL2 sets the initial conditions for the second Lagrange point of the Sun Earth system
    ///

    // begin with a null vector
    Eigen::Vector6d L2Initial;
    L2Initial(0) = 0.0;
    L2Initial(1) = 0.0;
    L2Initial(2) = 0.0;
    L2Initial(3) = 0.0;
    L2Initial(4) = 0.0;
    L2Initial(5) = 0.0;

    // Create object containing the functions.
    boost::shared_ptr< LibrationPointLocationFunction2SE > LibrationPointLocationFunction = boost::make_shared< LibrationPointLocationFunction2SE >( 1 );

    // The termination condition.
    tudat::root_finders::NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double > >(
                             LibrationPointLocationFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );
    // Test Newton-Raphson object.
    tudat::root_finders::NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    long double gammaL = newtonRaphson.execute( LibrationPointLocationFunction, LibrationPointLocationFunction->getTrueRootLocation() );

    // Apply manual corrections, based on trial and error, to ensure the stability of the point
    L2Initial [xCartesianPositionIndex] = gammaL * AU + 507.776127;
    L2Initial [yCartesianVelocityIndex] = L2Initial[xCartesianPositionIndex] * 2.0 * pi / secInYear - 0.000448950212 ;

    return L2Initial;
}

Eigen::VectorXd InitialStateEML2()
{
    ///
    /// InitialStateEML2 sets the initial conditions for the second Lagrange point of the Earth Moon system
    ///

    // begin with a null vector
    Eigen::Vector6d L2Initial;
    L2Initial(0) = 0.0;
    L2Initial(1) = 0.0;
    L2Initial(2) = 0.0;
    L2Initial(3) = 0.0;
    L2Initial(4) = 0.0;
    L2Initial(5) = 0.0;

    // Create object containing the functions.
    boost::shared_ptr< LibrationPointLocationFunction2EM > LibrationPointLocationFunction = boost::make_shared< LibrationPointLocationFunction2EM >( 1 );

    // The termination condition.
    tudat::root_finders::NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double > >(
                             LibrationPointLocationFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );
    // Test Newton-Raphson object.
    tudat::root_finders::NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    double gammaL = newtonRaphson.execute( LibrationPointLocationFunction, LibrationPointLocationFunction->getTrueRootLocation( ) );

    // Apply manual corrections, based on trial and error, to ensure the stability of the point
    L2Initial [xCartesianPositionIndex] = MoonEarthDistance * gammaL + 9.9997183293e4;
    L2Initial [yCartesianVelocityIndex] = L2Initial[xCartesianPositionIndex] * 2.0 * pi / secInMonth - 0.95575;

    return L2Initial;

}
