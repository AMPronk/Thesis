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

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getOutputPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return ( filePath_.substr( 0, filePath_.length( ) -
                               std::string( "WP1.cpp" ).length( ) ) ) + "PropagationResults";
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
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::ephemerides;

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started environment and vehicle creation."<<std::endl;

    // Set simulation time settings.
    const double simulationStartEpoch = physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = physical_constants::JULIAN_YEAR + 50.0 * physical_constants::JULIAN_DAY;
    double InitialStepSize = 1;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Olfar" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Olfar" ]->setConstantBodyMass( 3.5 );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started acceleration creation."<<std::endl;

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfOlfar;


    accelerationsOfOlfar[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );

    accelerationsOfOlfar[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                 basic_astrodynamics::central_gravity ) );

    accelerationsOfOlfar[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                  basic_astrodynamics::central_gravity ) );

    accelerationMap[ "Olfar" ] = accelerationsOfOlfar;
    bodiesToPropagate.push_back( "Olfar" );

    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////     SEGMENT 1: EARTH to SE-L2       ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment 1 calculations."<<std::endl;

    ///
    /// Set the moment at which the segment propagation starts. Note that backward propagation is used, thus the propagation start time is the end time for the trajectory
    /// Use propagation start to find the location of the sun at this point. This is later used to switch between coordinate systems.
    ///

    const double PropagationStartS1 = physical_constants::JULIAN_YEAR + 0.001 * physical_constants::JULIAN_DAY;
    InitialStepSize = -1;

    Eigen::VectorXd sunBarycentricStatesS1 = bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris(PropagationStartS1);
    double sunDistance = pow( pow( sunBarycentricStatesS1[1] , 2 ) + pow( sunBarycentricStatesS1[2] , 2 ) , 0.5 );
    double sunAngle;
    if (sunBarycentricStatesS1[1] <= 0 )
    {
        sunAngle = std::atan(sunBarycentricStatesS1[2]/sunBarycentricStatesS1[1]) + 2 * std::acos(0.0) ;
    }
    else
    {
        sunAngle = std::atan(sunBarycentricStatesS1[2]/sunBarycentricStatesS1[1]);
    }
    double rotationAngle = 2 * std::acos(0.0) - sunAngle;
    double rotationSpeed = 4 * std::acos(0.0) / (365.2422 * 24 * 60 * 60);

    ///
    /// Set the integrator settings, which are the same for every set of initial conditions
    ///
    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double minimumStepSizeS1 = -1.0E-5;
    double maximumStepSizeS1 = -10;
    double relativeErrorToleranceS1 = 1.0E-8;
    double absoluteErrorToleranceS1 = 1.0E-8;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1 =
            boost::make_shared< IntegratorSettings< > >
            (rungeKutta4, PropagationStartS1, InitialStepSize );
            //boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            //(rungeKuttaVariableStepSize, PropagationStartS1, InitialStepSize, coefficientSetS1,
             //minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );


    ///
    /// Calculation of the required constants to calculate the possible initial conditions -- performed in the CR3BP ///
    ///
    const double sunMass = 1.989E30;
    const double earthMass = 5.972E24;
    const double earthSunMu = earthMass/(earthMass+sunMass);
    const double L2Distance = 1.0 + pow( earthMass / (3.0 * sunMass) , 1.0/3.0);

    const double rho = earthSunMu * pow(L2Distance - 1.0 + earthSunMu,-3) + (1.0 - earthSunMu) * pow(L2Distance + earthSunMu,-3);
    const double a = 2.0 * rho + 1.0;
    const double b = rho - 1.0;

    const double labda = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
    const double nu = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

    Eigen::Vector4d u1;
    u1[0] = 1;
    u1[1] = u1[0] * ( pow( labda , 2 ) - a ) / ( 2*labda );
    u1[2] = u1[0] * labda;
    u1[3] = u1[1] * labda;

    Eigen::Vector4d u2;
    u2[0] = u1[0];
    u2[1] = - u1[1];
    u2[2] = - u1[2];
    u2[3] = u1[3];

    Eigen::Vector4cd w1;
    w1[0] = 1.0;
    w1[1] = w1[0] * ( -pow( nu , 2.0 ) - a ) / (2 * nu * std::complex<double>(0,-1));
    w1[2] = w1[0] * nu * std::complex<double>(0,-1);
    w1[3] = w1[1] * nu * std::complex<double>(0,-1);

    //std::cerr<< "w1 =  " + boost::lexical_cast< std::string >( w1 ) <<std::endl;


    /// Calculation of the possible initial vectors in the CR3BP, using the equation for ASYMPTOTIC orbits:
    /// [x, y, x', y'] = alpha1 * e ^ (labda*t) * u1 + 2 * beta * e ^ (nu*t) * w1
    /// Alpha1 and Beta are the sliding variables, which create the set of initial conditions.

    // Set the number of samples on each dimension
    int sampleSize = 2;

    // Create the vectors of each of the sliding variables. For segment 1, alpha must be smaller than 0
    Eigen::VectorXd alpha1 = Eigen::VectorXd::LinSpaced(sampleSize, 0.0, 0.0001);
    //double alpha2 = 0.0;
    Eigen::VectorXd beta = Eigen::VectorXd::LinSpaced(sampleSize,0.0,0.0001);
    double time = 0.0;
    int conditionsCounter;

    // Create the variables that are used to save the inital conditions in the CR3BP

    double x;
    double y;
    double xdot;
    double ydot;
    double x2;
    double y2;
    double xdot2;
    double ydot2;

    // Create the two matrix that will hold the initial conditions in the frame used for propagation.
    // initialConditionsS1 saves only x, y, xdot and ydot, along with the used alpha and beta. This is used for later reference.
    // initialCartesianElementsS1 holds the full set of Cartesian Elements and is used for the propagation. It is overwritten each sequence.
    Eigen::MatrixXd initialConditionsS1 ( 6 , sampleSize * sampleSize );
    Eigen::Vector6d initialCartesianElementsS1;

    // Create the matrix that will hold the output trajectories for each compination of alpha and beta.
    Eigen::MatrixXd cartesianIntegrationResultS1;

//    for( int alpha1Counter = 0; alpha1Counter <= sampleSize; alpha1Counter++ )
//    {
//        std::cerr<<"Next alpha."<<std::endl;
//        for( int betaCounter = 0; betaCounter <= sampleSize; betaCounter++ )
//        {
//            std::cerr<<"Next beta."<<std::endl;
//            conditionsCounter = sampleSize * alpha1Counter + betaCounter;

//            ///
//            /// for the given alpha and beta, using t=0, set the initial conditions
//            ///
//std::cerr<< "CR3BP coordinates" <<std::endl;
//            // initial conditions in the CR3BP
//            x = alpha1[alpha1Counter] * exp( labda * time ) * u1[0] + 2 * std::real( beta[betaCounter] * exp( nu * sqrt(std::complex<double>(-1,0)) * time) * w1[0] );
//            y = alpha1[alpha1Counter] * exp( labda * time ) * u1[1] + 2 * std::real( beta[betaCounter] * exp( nu * sqrt(std::complex<double>(-1,0)) * time) * w1[1] );
//            xdot = alpha1[alpha1Counter] * exp( labda * time ) * u1[2] + 2 * std::real( beta[betaCounter] * exp( nu * sqrt(std::complex<double>(-1,0)) * time) * w1[2] );
//            ydot = alpha1[alpha1Counter] * exp( labda * time ) * u1[3] + 2 * std::real( beta[betaCounter] * exp( nu * sqrt(std::complex<double>(-1,0)) * time) * w1[3]);


//std::cerr<< "Corrected coordinates" <<std::endl;
//            // initial conditions in the earth centered frame
//            initialConditionsS1(conditionsCounter, 0) = alpha1[alpha1Counter];
//            initialConditionsS1(conditionsCounter, 1) = beta[betaCounter];

//            initialConditionsS1(conditionsCounter, 2) = sunDistance * ( ( x + L2Distance - 1.0 + earthSunMu ) * std::cos(rotationAngle) + y * std::sin(rotationAngle) );
//            initialConditionsS1(conditionsCounter, 3) = sunDistance * ( y * std::cos(rotationAngle) - ( x + L2Distance - 1.0 + earthSunMu ) * std::sin(rotationAngle) );

//            initialConditionsS1(conditionsCounter, 4) = sunDistance * ( xdot * std::cos(rotationAngle) + ydot * std::sin(rotationAngle) )
//                    - rotationSpeed * sqrt( pow (initialConditionsS1(conditionsCounter,2) , 2.0 ) + pow( initialConditionsS1(conditionsCounter, 3) , 2.0 )) * sin(-rotationAngle);

//            initialConditionsS1(conditionsCounter, 5) = sunDistance * ( ydot * std::cos(rotationAngle) - xdot * std::sin(rotationAngle) )
//                    + rotationSpeed * sqrt( pow (initialConditionsS1(conditionsCounter,2) , 2.0 ) + pow( initialConditionsS1(conditionsCounter, 3) , 2.0 )) * cos(-rotationAngle);

//            //initial conditions in the earth centered frame, as used for the propagation
//            initialCartesianElementsS1[ xCartesianPositionIndex ] = initialConditionsS1(conditionsCounter, 2);
//            initialCartesianElementsS1[ yCartesianPositionIndex ] = initialConditionsS1(conditionsCounter, 3);
//            initialCartesianElementsS1[ zCartesianPositionIndex ] = 0.0;
//            initialCartesianElementsS1[ xCartesianVelocityIndex ] = initialConditionsS1(conditionsCounter, 4);
//            initialCartesianElementsS1[ yCartesianVelocityIndex ] = initialConditionsS1(conditionsCounter, 5);
//            initialCartesianElementsS1[ zCartesianVelocityIndex ] = 0.0;

//            ///
//            /// Set the propagator settings
//            ///

//            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
//                    boost::make_shared< TranslationalStatePropagatorSettings< > >
//                    (centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianElementsS1, simulationStartEpoch, cowell );

//            ///
//            /// Create simulation object and propagate dynamics.
//            ///
//            std::cerr<< "Start Propagator" <<std::endl;
//            SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
//                        bodyMap, integratorSettingsS1, propagatorSettingsS1, true, false, false );
//            std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );

//            std::cerr<< "Start Writing to file" <<std::endl;

//            // Write segment 1 propagation to file for later access.
//            input_output::writeDataMapToTextFile( tempIntegrationResultS1,
//                                                  "WP1_Segment1_alpha1" +
//                                                  boost::lexical_cast< std::string >( alpha1[alpha1Counter] ) + "_beta" +
//                                                  boost::lexical_cast< std::string >( beta[betaCounter] ) + ".dat",
//                                                  tudat_applications::getOutputPath( ),
//                                                  "",
//                                                  std::numeric_limits< double >::digits10,
//                                                  std::numeric_limits< double >::digits10,
//                                                  "," );

//        }
//    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////     SEGMENT 2: SE-L2 to Earth-Cut       ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment 2 calculations."<<std::endl;

    ///
    /// Set the moment at which the segment propagation starts and ends.
    /// Use propagation start to find the location of the sun at this point. This is later used to switch between coordinate systems.
    ///

    const double PropagationStartS2 = physical_constants::JULIAN_YEAR + 2 * physical_constants::JULIAN_DAY;
    const double PropagationEndsS2 = physical_constants::JULIAN_YEAR + 2.01 * physical_constants::JULIAN_DAY;
    InitialStepSize = 1E-5;

    Eigen::VectorXd sunBarycentricStatesS2 = bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris(PropagationStartS2);
    sunDistance = pow( pow( sunBarycentricStatesS2[1] , 2 ) + pow( sunBarycentricStatesS2[2] , 2 ) , 0.5 );

    if (sunBarycentricStatesS2[1] <= 0 )
    {
        sunAngle = std::atan(sunBarycentricStatesS2[2]/sunBarycentricStatesS2[1]) + 2 * std::acos(0.0) ;
    }
    else
    {
        sunAngle = std::atan(sunBarycentricStatesS2[2]/sunBarycentricStatesS2[1]);
    }
    rotationAngle = 2 * std::acos(0.0) - sunAngle;
    rotationSpeed = 4 * std::acos(0.0) / (365.2422 * 24 * 60 * 60);

    ///
    /// Set the integrator settings, which are the same for every set of initial conditions
    ///
    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS2 =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double minimumStepSizeS2 = 1.0E-10;
    double maximumStepSizeS2 = 1.0E-2;
    //double relativeErrorToleranceS2 = 1.0;
    //double absoluteErrorToleranceS2 = 1.0;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS2 =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            (rungeKuttaVariableStepSize, PropagationStartS2, InitialStepSize, coefficientSetS2,
             minimumStepSizeS2, maximumStepSizeS2);//, relativeErrorToleranceS2, absoluteErrorToleranceS2 );


    ///
    /// Calculation of the required constants to calculate the possible initial conditions -- performed in the CR3BP ///
    /// TAKEN OUT: SAME AS FOR SEGMENT 1
    ///
//    const double sunMass = 1.989E30;
//    const double earthMass = 5.972E24;
//    const double earthSunMu = earthMass/(earthMass+sunMass);
//    const double L2Distance = 1.0 + pow( earthMass / (3.0 * sunMass) , 1.0/3.0);

//    const double rho = earthSunMu * pow(L2Distance - 1.0 + earthSunMu,-3) + (1.0 - earthSunMu) * pow(L2Distance + earthSunMu,-3);
//    const double a = 2.0 * rho + 1.0;
//    const double b = rho - 1.0;

//    const double labda = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
//    const double nu = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

//    Eigen::Vector4d u1;
//    u1[0] = 1;
//    u1[1] = u1[0] * ( pow( labda , 2 ) - a ) / ( 2*labda );
//    u1[2] = u1[0] * labda;
//    u1[3] = u1[1] * labda;

//    Eigen::Vector4d u2;
//    u2[0] = u1[0];
//    u2[1] = - u1[1];
//    u2[2] = - u1[2];
//    u2[3] = u1[3];

//    Eigen::Vector4cd w1;
//    w1[0] = 1.0;
//    w1[1] = w1[0] * ( -pow( nu , 2.0 ) - a ) / (2 * nu * sqrt(std::complex<double>(-1,0)));
//    w1[2] = w1[0] * nu * sqrt(std::complex<double>(-1,0));
//    w1[3] = w1[1] * nu * sqrt(std::complex<double>(-1,0));

//    //std::cerr<< "w1 =  " + boost::lexical_cast< std::string >( w1 ) <<std::endl;


    /// Calculation of the possible initial vectors in the CR3BP, using the equation for ASYMPTOTIC orbits:
    /// [x, y, x', y'] = alpha1 * e ^ (labda*t) * u1 + 2 * beta * e ^ (nu*t) * w1
    /// Alpha1 and Beta are the sliding variables, which create the set of initial conditions.

    // Set the number of samples on each dimension
    sampleSize = 10;

    // Create the vectors of each of the sliding variables. For segment 1, alpha must be smaller than 0
    Eigen::VectorXd alpha2 = Eigen::VectorXd::LinSpaced(sampleSize, -1E-16, 1E-16);
    //double alpha1 = 0.0;
    beta = Eigen::VectorXd::LinSpaced(sampleSize, -1E-4, 1E-4);
    time = 0.0;

    // Create the two matrix that will hold the initial conditions in the frame used for propagation.
    // initialConditionsS2 saves only x, y, xdot and ydot, along with the used alpha and beta. This is used for later reference.
    // initialCartesianElementsS2 holds the full set of Cartesian Elements and is used for the propagation. It is overwritten each sequence.
    Eigen::MatrixXd initialConditionsS2 ( 6 , sampleSize * sampleSize );
    Eigen::Vector6d initialCartesianElementsS2;

    // Create the matrix that will hold the output trajectories for each compination of alpha and beta.
    Eigen::MatrixXd cartesianIntegrationResultS2;

    // Create a variable that will hold an angle;
    double gamma;

    for( int alpha2Counter = 0; alpha2Counter <= sampleSize; alpha2Counter++ )
    {
        std::cerr<<"Next alpha."<<std::endl;
        for( int betaCounter = 0; betaCounter <= sampleSize; betaCounter++ )
        {
            std::cerr<<"Next beta."<<std::endl;
            conditionsCounter = sampleSize * alpha2Counter + betaCounter;

            ///
            /// for the given alpha and beta, using t=0, set the initial conditions
            ///
//std::cerr<< "CR3BP coordinates" <<std::endl;
            // initial conditions in the CR3BP
            x = alpha2[alpha2Counter] * exp( - labda * time ) * u1[0] + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time) * w1[0] );
            y = alpha2[alpha2Counter] * exp( - labda * time ) * u1[1] + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time) * w1[1] );
            xdot = alpha2[alpha2Counter] * exp( - labda * time ) * u1[2] + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time) * w1[2] );
            ydot = alpha2[alpha2Counter] * exp( - labda * time ) * u1[3] + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time) * w1[3]);

            // scale the coordinates
            x2 = sunDistance * (x + L2Distance -1.0 + earthSunMu );
            y2 = sunDistance * y;
            xdot2 = sunDistance * xdot;
            ydot2 = sunDistance * ydot;

            // add rotation and rotational speed to get initial conditions in the earth centered frame
            initialConditionsS2(conditionsCounter, 0) = alpha2[alpha2Counter];
            initialConditionsS2(conditionsCounter, 1) = beta[betaCounter];

            initialConditionsS2(conditionsCounter, 2) = x2 * std::cos(rotationAngle) + y2 * std::sin(rotationAngle) ;
            initialConditionsS2(conditionsCounter, 3) = y2 * std::cos(rotationAngle) - x2 * std::sin(rotationAngle) ;

            initialConditionsS2(conditionsCounter, 4) = xdot2 * std::cos(rotationAngle) + ydot2 * std::sin(rotationAngle) + rotationSpeed * (- std::sin(rotationAngle) * x2 + std::cos(rotationAngle) * y2);
            initialConditionsS2(conditionsCounter, 5) = ydot2 * std::cos(rotationAngle) - xdot2 * std::sin(rotationAngle) + rotationSpeed * (- std::cos(rotationAngle) * x2 - std::sin(rotationAngle) * y2) ;

            //initial conditions in the earth centered frame, as used for the propagation
            initialCartesianElementsS2[ xCartesianPositionIndex ] = initialConditionsS2(conditionsCounter, 2);
            initialCartesianElementsS2[ yCartesianPositionIndex ] = initialConditionsS2(conditionsCounter, 3);
            initialCartesianElementsS2[ zCartesianPositionIndex ] = 0.0;
            initialCartesianElementsS2[ xCartesianVelocityIndex ] = initialConditionsS2(conditionsCounter, 4);
            initialCartesianElementsS2[ yCartesianVelocityIndex ] = initialConditionsS2(conditionsCounter, 5);
            initialCartesianElementsS2[ zCartesianVelocityIndex ] = 0.0;

            ///
            /// Set the propagator settings
            ///

            boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS2 =
                    boost::make_shared< TranslationalStatePropagatorSettings< > >
                    (centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianElementsS2, PropagationEndsS2, cowell );

            ///
            /// Create simulation object and propagate dynamics.
            ///
            //std::cerr<< "Start Propagator" <<std::endl;
            SingleArcDynamicsSimulator< > dynamicsSimulatorS2(
                        bodyMap, integratorSettingsS2, propagatorSettingsS2, true, false, false );
            std::map< double, Eigen::VectorXd > tempIntegrationResultS2 = dynamicsSimulatorS2.getEquationsOfMotionNumericalSolution( );

            //std::cerr<< "Start Writing to file" <<std::endl;

            // Write segment 2 propagation to file for later access.
            input_output::writeDataMapToTextFile( tempIntegrationResultS2,
                                                  "WP1_Segment2_alpha2" +
                                                  boost::lexical_cast< std::string >( alpha2[alpha2Counter] ) + "_beta" +
                                                  boost::lexical_cast< std::string >( beta[betaCounter] ) + ".dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

        }
    }



    // Set initial Cartesian elements for Olfar.
    Eigen::Vector6d initialCartesianElementsS3;

        /// FINAL CONDITIONS OF TRAJECTORY SEGMENT 3: Moon - Earth L2 elliptic orbit, asymptotic trajectory energy ///
    initialCartesianElementsS3[ xCartesianPositionIndex ] = 6378.0E3 + 550.0E3;
    initialCartesianElementsS3[ yCartesianPositionIndex ] = 0.05; // NOTE: This value is larger than for the real satellite
    initialCartesianElementsS3[ zCartesianPositionIndex ] = 0.0;
    initialCartesianElementsS3[ xCartesianVelocityIndex ] = 1000.0;
    initialCartesianElementsS3[ yCartesianVelocityIndex ] = 8000.0;
    initialCartesianElementsS3[ zCartesianVelocityIndex ] = 0.0;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////       PROPAGATE SEGMENT 3          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment 3 propagation."<<std::endl;

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS3 =
            boost::make_shared< TranslationalStatePropagatorSettings< > >
            (centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianElementsS3, simulationEndEpoch, cowell );

    const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS3 =
            numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    double minimumStepSizeS3 = 1.0E-9;
    double maximumStepSizeS3 = 1.0E6;
    double relativeErrorToleranceS3 = 1.0E-14;
    double absoluteErrorToleranceS3 = 1.0E14;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS3 =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            (rungeKuttaVariableStepSize, simulationStartEpoch, InitialStepSize, coefficientSetS3,
             minimumStepSizeS3, maximumStepSizeS3, relativeErrorToleranceS3, absoluteErrorToleranceS3 );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulatorS3(
                bodyMap, integratorSettingsS3, propagatorSettingsS3, true, false, false );
    std::map< double, Eigen::VectorXd > cartesianIntegrationResultS3 = dynamicsSimulatorS3.getEquationsOfMotionNumericalSolution( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        LINK SEGMENTS INTO TRAjECTORIES               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment linking."<<std::endl;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Writing output to file."<<std::endl;

    //    // Write Olfar propagation history to file.
    //    input_output::writeDataMapToTextFile( cartesianIntegrationResultS1,
    //                                          "OlfarOrbitCartesian.dat",
    //                                          tudat_applications::getOutputPath( ),
    //                                          "",
    //                                          std::numeric_limits< double >::digits10,
    //                                          std::numeric_limits< double >::digits10,
    //                                          "," );

    //    // Write Olfar propagation history to file.
    //    input_output::writeDataMapToTextFile( keplerianIntegrationResult,
    //                                          "OlfarOrbitKeplerian.dat",
    //                                          tudat_applications::getOutputPath( ),
    //                                          "",
    //                                          std::numeric_limits< double >::digits10,
    //                                          std::numeric_limits< double >::digits10,
    //                                          "," );

    //    std::map< double, Eigen::VectorXd > sunBarycentricStates;
    //    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResultS1.begin( );
    //         stateIterator != cartesianIntegrationResultS1.end( ); stateIterator++ )
    //    {
    //        sunBarycentricStates[ stateIterator->first ] = bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris(
    //                    stateIterator->first );
    //    }

    //    // Write Earth barycentric state history to file.
    //    input_output::writeDataMapToTextFile( sunBarycentricStates,
    //                                          "SunBarycentricStatesPart2.dat",
    //                                          tudat_applications::getOutputPath( ),
    //                                          "",
    //                                          std::numeric_limits< double >::digits10,
    //                                          std::numeric_limits< double >::digits10,
    //                                          "," );

    std::cerr<<"COMPLETED."<<std::endl;

}


