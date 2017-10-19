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

std::tuple<AccelerationMap, std::vector< std::string >, std::vector< std::string >, NamedBodyMap> CreateEphemeris(bool sun, bool earth, bool moon, bool SCSE, bool SCEM)
{

    unsigned int NrOfBodies = 0;


    NamedBodyMap bodyMap;

    if(sun)
    {
        NrOfBodies = NrOfBodies + 1;

        // Create sun point mass
        bodyMap [ "Sun" ] = boost::make_shared< simulation_setup::Body >();
        bodyMap [ "Sun" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( SunGravitationalParameter ) );
    }

    if(earth)
    {
        NrOfBodies = NrOfBodies + 1;

        // Create earth point mass
        bodyMap [ "Earth" ] = boost::make_shared< simulation_setup::Body >();
        bodyMap [ "Earth" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( EarthGravitationalParameter ) );
    }

    if(moon)
    {
        NrOfBodies = NrOfBodies + 1;

        // Create moon point mass
        bodyMap [ "Moon" ] = boost::make_shared< simulation_setup::Body >();
        bodyMap [ "Moon" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( MoonGravitationalParameter ) );
    }

    if(SCSE)
    {
        NrOfBodies = NrOfBodies + 1;

        // Create the spacecraft point mass that is attracted by the sun and the earth
        bodyMap[ "Spacecraft1" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Spacecraft1" ]->setConstantBodyMass( 5.0 );
    }
    if(SCEM)
    {
        NrOfBodies = NrOfBodies + 1;

        // Create the spacecraft point mass that is attracted by the earth and the moon
        bodyMap[ "Spacecraft2" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Spacecraft2" ]->setConstantBodyMass( 5.0 );
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

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSun;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSCSE;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSCEM;

    if(sun){
        if(earth){
            accelerationsOfSun[ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Sun" ] = accelerationsOfSun;
        bodiesToPropagate.push_back( "Sun" );
        centralBodies.push_back( "SSB" );
    }

    if(earth){
        if(sun){
            accelerationsOfEarth[ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Earth" ] = accelerationsOfEarth;
        bodiesToPropagate.push_back( "Earth" );
        centralBodies.push_back( "SSB" );
    }

    if(moon){
        if(earth){
            accelerationsOfMoon[ "Earth" ].push_back(
                    boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Moon" ] = accelerationsOfMoon;
        bodiesToPropagate.push_back( "Moon" );
        centralBodies.push_back( "Earth" );
    }

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
        centralBodies.push_back( "SSB" );
    }

    if(SCEM){
        if(moon){
            accelerationsOfSCEM[ "Moon" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        if(earth){
            accelerationsOfSCEM[ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        accelerationMap[ "Spacecraft2" ] = accelerationsOfSCEM;
        bodiesToPropagate.push_back( "Spacecraft2" );
        centralBodies.push_back( "Earth" );
    }

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    std::cerr<<"End CreateEphemeris."<<std::endl;

    return std::make_tuple(accelerationModelMap, bodiesToPropagate, centralBodies, bodyMap);
}


Eigen::VectorXd CreateEphemerisInitialState(bool sun, bool earth, bool moon)
{

    Eigen::Vector6d SunInitialState;
    SunInitialState [xCartesianPositionIndex] = - muSE * AU;
    SunInitialState [yCartesianPositionIndex] = 0.0;
    SunInitialState [zCartesianPositionIndex] = 0.0;
    SunInitialState [xCartesianVelocityIndex] = 0.0;
    SunInitialState [yCartesianVelocityIndex] = SunInitialState[xCartesianPositionIndex] * 2.0 * pi / (secInYear);
    SunInitialState [zCartesianVelocityIndex] = 0.0;

    Eigen::Vector6d EarthInitialState;
    EarthInitialState [xCartesianPositionIndex] = (1.0 - muSE) * AU;
    EarthInitialState [yCartesianPositionIndex] = 0.0;
    EarthInitialState [zCartesianPositionIndex] = 0.0;
    EarthInitialState [xCartesianVelocityIndex] = 0.0;
    EarthInitialState [yCartesianVelocityIndex] = EarthInitialState[xCartesianPositionIndex] * 2.0 * pi / (secInYear);
    EarthInitialState [zCartesianVelocityIndex] = 0.0;

    Eigen::Vector6d MoonInitialState;
    MoonInitialState [xCartesianPositionIndex] = - MoonEarthDistance;
    MoonInitialState [yCartesianPositionIndex] = 0.0;
    MoonInitialState [zCartesianPositionIndex] = 0.0;
    MoonInitialState [xCartesianVelocityIndex] = 0.0;
    MoonInitialState [yCartesianVelocityIndex] = MoonInitialState[xCartesianPositionIndex] * 2.0 * pi / secInMonth;
    MoonInitialState [zCartesianVelocityIndex] = 0.0;

    int nrOfBodies = sun + earth + moon;
    Eigen::VectorXd systemInitialState = Eigen::VectorXd(nrOfBodies * 6);

    //begin constructing system initial state
    if(sun){
        if(earth){
            if(moon){
                systemInitialState.segment(0,6) = SunInitialState;
                systemInitialState.segment(6,6) = EarthInitialState;
                systemInitialState.segment(12,6) = MoonInitialState;
            }
            else{
                systemInitialState.segment(0,6) = SunInitialState;
                systemInitialState.segment(6,6) = EarthInitialState;
            }
        }
        else{
            if(moon){
                systemInitialState.segment(0,6) = SunInitialState;
                systemInitialState.segment(6,6) = MoonInitialState;
            }
            else{
                systemInitialState.segment(0,6) = SunInitialState;
            }
        }
    }
    else{
        if(earth){
            if(moon){
                systemInitialState.segment(0,6) = EarthInitialState;
                systemInitialState.segment(6,6) = MoonInitialState;
            }
            else{
                systemInitialState.segment(0,6) = EarthInitialState;
            }
        }
        else{
            if(moon){
                systemInitialState.segment(0,6) = MoonInitialState;
            }
            else{

            }
        }
    }

    return systemInitialState;
}

Eigen::VectorXd InitialStateSEL2()
{
    //using equations from Astrodynamics reader, p66

    long double alpha = muSE/(1.0-muSE);
    long double beta = pow(alpha/3.0 , 1.0/3.0);
    long double gamma = beta + pow(beta,2.0) / 3.0 - pow(beta,3.0) / 9.0 - pow(beta,4.0) * 31.0/81.0;

    Eigen::VectorXd earthInitial = CreateEphemerisInitialState(false,true,false);

    Eigen::VectorXd L2Initial = earthInitial;

    long double z = pow( muSE/3.0 , 1.0/3.0 );
    long double unscaledX = z + pow(z,2.0)/3.0 - pow(z,3.0)/9.0 + pow(z,4.0)*50.0/81.0 + 1.0 - muSE;

    // Create object containing the functions.
    boost::shared_ptr< LibrationPointLocationFunction2 > LibrationPointLocationFunction = boost::make_shared< LibrationPointLocationFunction2 >( 1 );

    // The termination condition.
    tudat::root_finders::NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double > >(
                                 LibrationPointLocationFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );
    // Test Newton-Raphson object.
    tudat::root_finders::NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    double gammaL = newtonRaphson.execute( LibrationPointLocationFunction, LibrationPointLocationFunction->getInitialGuess( ) );

    L2Initial [xCartesianPositionIndex] = (gammaL + 1 - muSE) * AU; //unscaledX * AU; // (1.0-muSE+gamma)*AU;
    L2Initial [yCartesianVelocityIndex] = L2Initial[xCartesianPositionIndex] * 2.0 * pi / (secInYear);

    return L2Initial;
}

Eigen::VectorXd InitialStateEML2()
{

    //using equations from Astrodynamics reader, p66
    long double alpha = muEM/(1.0-muEM);
    long double beta = pow(alpha/3.0 , 1.0/3.0);
    long double gamma = beta + pow(beta,2.0) / 3.0 - pow(beta,3.0) / 9.0 - pow(beta,4.0) * 31.0/81.0;

    Eigen::VectorXd moonInitial = CreateEphemerisInitialState(false,false,true);

    Eigen::VectorXd L2Initial = moonInitial;

    L2Initial [xCartesianPositionIndex] = moonInitial[xCartesianPositionIndex] * (1.0 + gamma);
    L2Initial [yCartesianVelocityIndex] = - 1023.0 / moonInitial[xCartesianPositionIndex] * L2Initial[xCartesianPositionIndex];

    return L2Initial;

}
