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
    //    // Declare file path string assigned to filePath.
    //    // __FILE__ only gives the absolute path of the header file!
    //    std::string filePath_( __FILE__ );

    //    // Strip filename from temporary string and return root-path string.
    //    return ( filePath_.substr( 0, filePath_.length( ) -
    //                               std::string( "WP1.cpp" ).length( ) ) ) + "PropagationResults";
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

    const double pi = 2 * std::acos(0.0) ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started environment and vehicle creation."<<std::endl;

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 4;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Sun";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Moon";
    bodyNames[ 3 ] = "Olfar";

    // Gravitational Parameters from SMAD, p955
    double SunGravitationalParameter = 1.32712440041E20;
    double EarthGravitationalParameter = 3.986004356E14;
    double MoonGravitationalParameter = 4.90280015E12;

    //    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
    //            getDefaultBodySettings(bodyNames);
    NamedBodyMap bodyMap; // = createBodies( bodySettings );

    // Create sun point mass
    bodyMap [ "Sun" ] = boost::make_shared< simulation_setup::Body >();
    //bodyMap [ "Sun" ]-> setConstantBodyMass( 1.989E30 );
    bodyMap [ "Sun" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( SunGravitationalParameter ) );

    // Create earth point mass
    bodyMap [ "Earth" ] = boost::make_shared< simulation_setup::Body >();
    //bodyMap [ "Earth" ]-> setConstantBodyMass( 5.972E24 );
    bodyMap [ "Earth" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( EarthGravitationalParameter ) );

    // Create moon point mass
    bodyMap [ "Moon" ] = boost::make_shared< simulation_setup::Body >();
    //bodyMap [ "Moon" ]-> setConstantBodyMass( 7.3477E22 );
    bodyMap [ "Moon" ]-> setGravityFieldModel( boost::make_shared< GravityFieldModel >( MoonGravitationalParameter ) );

    // Create Olfar point mass
    bodyMap[ "Olfar" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Olfar" ]->setConstantBodyMass( 5.0 );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started acceleration creation."<<std::endl;

    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;

    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > currentAccelerations;
        for( unsigned int j = 0; j < bodyNames.size( ) - 1 ; j++ )
        {
            // Create central gravity acceleration between each 2 bodies.
            if( i != j )
            {
                currentAccelerations[ bodyNames.at( j ) ].push_back(
                            boost::make_shared< AccelerationSettings >( central_gravity ) );\
            }
        }
        accelerationMap[ bodyNames.at( i ) ] = currentAccelerations;
    }

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToPropagate = bodyNames;
    unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );

    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.resize( numberOfNumericalBodies );

    // Set central bodies with hyrarchical setting
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        // Set SSB as central body for Sun
        if( i == 0 )
        {
            centralBodies[ i ] = "SSB";
        }
        // Set Sun as central body for Earth
        else if( i == 1 )
        {
            centralBodies[ i ] = "Sun";
        }
        // Set Earth as central body for all objects
        else
        {
            centralBodies[ i ] = "Earth";
        }
    }

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////     SET INITIAL CONDITIONS OF BODIES       /////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::Vector6d SunInitialState;
    SunInitialState [xCartesianPositionIndex] = 0.0;
    SunInitialState [yCartesianPositionIndex] = 0.0;
    SunInitialState [zCartesianPositionIndex] = 0.0;
    SunInitialState [xCartesianVelocityIndex] = 0.0;
    SunInitialState [yCartesianVelocityIndex] = 0.0;
    SunInitialState [zCartesianVelocityIndex] = 0.0;

    //Astronomical Unit in meters
    const double AU = 149597870700;

    Eigen::Vector6d EarthInitialStateKeplerian;
    EarthInitialStateKeplerian( semiMajorAxisIndex) = 1.00000011 * AU;
    EarthInitialStateKeplerian( eccentricityIndex) = 0.0;
    EarthInitialStateKeplerian( inclinationIndex) = 0.0;
    EarthInitialStateKeplerian( argumentOfPeriapsisIndex) = 0.0;
    EarthInitialStateKeplerian( longitudeOfAscendingNodeIndex) = 0.0;
    EarthInitialStateKeplerian( trueAnomalyIndex ) = 0.0;

    Eigen::Vector6d MoonInitialStateKeplerian;
    MoonInitialStateKeplerian( semiMajorAxisIndex ) = 0.002570 * AU;
    MoonInitialStateKeplerian( eccentricityIndex ) = 0.0;
    MoonInitialStateKeplerian( inclinationIndex ) = 0.0;
    MoonInitialStateKeplerian( argumentOfPeriapsisIndex ) = 0.0;
    MoonInitialStateKeplerian( longitudeOfAscendingNodeIndex ) = 0.0;
    MoonInitialStateKeplerian( trueAnomalyIndex ) = 0.0;


    const Eigen::Vector6d EarthInitialState = convertKeplerianToCartesianElements(
                EarthInitialStateKeplerian,
                SunGravitationalParameter );

    const Eigen::Vector6d MoonInitialState = convertKeplerianToCartesianElements(
                MoonInitialStateKeplerian,
                EarthGravitationalParameter );

    //begin constructing system initial state
    Eigen::VectorXd systemInitialState = Eigen::VectorXd(24);
    systemInitialState.segment(0,6) = SunInitialState;
    systemInitialState.segment(6,6) = EarthInitialState;
    systemInitialState.segment(12,6) = MoonInitialState;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////     SEGMENT 1: EARTH to SE-L2       ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment 1 calculations."<<std::endl;

    ///
    /// Set the moment at which the segment propagation starts. Note that backward propagation is used, thus the propagation start time is the end time for the trajectory
    ///

    const double PropagationStartS1 = 0.0;
    const double PropagationEndsS1 = -0.5E7;
    double StepSize = -1000.0;

    /// Set the location of the Sun at PropagationStartS1
    double sunDistance = EarthInitialStateKeplerian(semiMajorAxisIndex);
    double sunAngle = pi;
    double rotationAngle = pi - sunAngle;
    double rotationSpeed = 2 * pi / (365.2422 * 24 * 60 * 60);

    ///
    /// Set the integrator settings, which are the same for every set of initial conditions
    ///
    //const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
    //        numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    //double minimumStepSizeS1 = -1.0E-5;
    //double maximumStepSizeS1 = -10;
    //double relativeErrorToleranceS1 = 1.0E-8;
    //double absoluteErrorToleranceS1 = 1.0E-8;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS1 =
            boost::make_shared< IntegratorSettings< > >
            (rungeKutta4, PropagationStartS1, StepSize );
    //boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
    //(rungeKuttaVariableStepSize, PropagationStartS1, InitialStepSize, coefficientSetS1,
    //minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );


    ///
    /// Calculation of the required constants to calculate the possible initial conditions -- performed in the CR3BP ///
    ///
    const double sunMass = 1.989E30;
    const double earthMass = 5.972E24;
    const double earthSunMu = earthMass/(earthMass+sunMass);
    double L2Distance = 1.0 + pow( earthMass / (3.0 * sunMass) , 1.0/3.0);

    double rho = earthSunMu * pow(L2Distance - 1.0 + earthSunMu,-3) + (1.0 - earthSunMu) * pow(L2Distance + earthSunMu,-3);
    double a = 2.0 * rho + 1.0;
    double b = rho - 1.0;

    double labda = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
    double nu = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

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


    /// Calculation of the possible initial vectors in the CR3BP, using the equation:
    /// [x, y, x', y'] = alpha1 * e ^ (labda*t) * u1 + alpha2 * e ^ (-labda*t) * u2 + 2 * beta * e ^ (nu*t) * w1
    /// Alpha1, Alpha2, Beta and Time are the sliding variables, which create the set of initial conditions.

    // Set the number of samples on each dimension
    int sampleSizeA1 = 2;
    int sampleSizeA2 = 1; //set to 0??
    int sampleSizeB = 2;
    int sampleSizeT = 1;

    // Create the vectors of each of the sliding variables. For segment 1, alpha1 must be zero and alpha2 must be larger than 0
    Eigen::VectorXd alpha1 = Eigen::VectorXd::LinSpaced(sampleSizeA1, -5E-10, 0.0);
    Eigen::VectorXd alpha2 = Eigen::VectorXd::LinSpaced(sampleSizeA2, -1E-300, 1E-300);
    Eigen::VectorXd beta = Eigen::VectorXd::LinSpaced(sampleSizeB, -2.0E-10, -0.65E-10);
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(sampleSizeT, 0.0, 10000);

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
    Eigen::Matrix4d initialConditionsS1;
    Eigen::Vector6d initialCartesianElementsS1;

    // Create the matrix that will hold the output trajectories for each compination of alpha and beta.
    Eigen::MatrixXd cartesianIntegrationResultS1;

    for(int alpha1Counter = 0; alpha1Counter <= alpha1.size(); alpha1Counter++)
    {
        std::cerr<<"S1, Alpha1 counter " + boost::lexical_cast< std::string >( alpha1Counter )<<std::endl;
        for(int alpha2Counter = 0; alpha2Counter <= alpha2.size(); alpha2Counter++)
        {
            //std::cerr<<"Alpha2+."<<std::endl;
            std::cerr<<"S1, Alpha2 counter " + boost::lexical_cast< std::string >( alpha2Counter )<<std::endl;
            for(int betaCounter = 0; betaCounter <= beta.size(); betaCounter++)
            {
                std::cerr<<"S1, Beta counter " + boost::lexical_cast< std::string >( betaCounter )<<std::endl;
                //std::cerr<<"Beta+."<<std::endl;
                for(int timeCounter = 0; timeCounter <= time.size(); timeCounter++)
                {

                    //std::cerr<<"Time+."<<std::endl;
                    ///
                    /// for the given alpha1, alpha2, beta and time, set the initial conditions
                    ///

                    // initial conditions in the CR3BP
                    x =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[0]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[0]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[0] );
                    y =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[1]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[1]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[1] );
                    xdot =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[2]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[2]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[2] );
                    ydot =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[3]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[3]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[3] );

                    // scale the coordinates
                    x2 = sunDistance * (x + L2Distance -1.0 + earthSunMu );
                    y2 = sunDistance * y;
                    xdot2 = sunDistance * xdot;
                    ydot2 = sunDistance * ydot;

                    // add rotation and rotational speed to get initial conditions in the earth centered frame

                    initialConditionsS1(0) = x2 * std::cos(rotationAngle) + y2 * std::sin(rotationAngle) ;
                    initialConditionsS1(1) = y2 * std::cos(rotationAngle) - x2 * std::sin(rotationAngle) ;

                    initialConditionsS1(2) = xdot2 * std::cos(rotationAngle) + ydot2 * std::sin(rotationAngle) + rotationSpeed * (- std::sin(rotationAngle) * x2 + std::cos(rotationAngle) * y2) ;
                    initialConditionsS1(3) = ydot2 * std::cos(rotationAngle) - xdot2 * std::sin(rotationAngle) + rotationSpeed * (- std::cos(rotationAngle) * x2 - std::sin(rotationAngle) * y2) ;

                    //initial conditions in the earth centered frame, as used for the propagation
                    initialCartesianElementsS1[ xCartesianPositionIndex ] = initialConditionsS1(0);
                    initialCartesianElementsS1[ yCartesianPositionIndex ] = initialConditionsS1(1);
                    initialCartesianElementsS1[ zCartesianPositionIndex ] = 0.0;
                    initialCartesianElementsS1[ xCartesianVelocityIndex ] = initialConditionsS1(2);
                    initialCartesianElementsS1[ yCartesianVelocityIndex ] = initialConditionsS1(3);
                    initialCartesianElementsS1[ zCartesianVelocityIndex ] = 0.0;

                    systemInitialState.segment(18,6) = initialCartesianElementsS1;

                    ///
                    /// Set the propagator settings
                    ///

                    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS1 =
                            boost::make_shared< TranslationalStatePropagatorSettings< > >
                            (centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, PropagationEndsS1, cowell );

                    ///
                    /// Create simulation object and propagate dynamics.
                    ///
                    SingleArcDynamicsSimulator< > dynamicsSimulatorS1(
                                bodyMap, integratorSettingsS1, propagatorSettingsS1, true, false, false );
                    std::map< double, Eigen::VectorXd > tempIntegrationResultS1 = dynamicsSimulatorS1.getEquationsOfMotionNumericalSolution( );

                    // Write segment 1 propagation to file for later access.
                    input_output::writeDataMapToTextFile( tempIntegrationResultS1,
                                                          "WP1_Segment1_alpha1=" + boost::lexical_cast< std::string >( alpha1[alpha1Counter] )
                                                          + "_alpha2=" + boost::lexical_cast< std::string >( alpha2[alpha2Counter] )
                                                          + "_beta=" + boost::lexical_cast< std::string >( beta[betaCounter] )
                                                          + "_time=" + boost::lexical_cast< std::string >( time[timeCounter] )
                                                          + ".dat",
                                                          tudat_applications::getOutputPath( ),
                                                          "",
                                                          std::numeric_limits< double >::digits10,
                                                          std::numeric_limits< double >::digits10,
                                                          "," );

                }
            }
        }
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////     SEGMENT 2: SE-L2 to Earth-Cut       ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cerr<<"Started segment 2 calculations."<<std::endl;

    ///
    /// Set the moment at which the segment propagation starts.
    ///

    const double PropagationStartS2 = 0.0;
    const double PropagationEndsS2 = 0.5E7;
    StepSize = 1000.0;

    /// if PropagationStartS1 = PropagationStartS2, sun variables remain the same
    //double sunDistance = EarthInitialStateKeplerian(semiMajorAxisIndex);
    //double sunAngle = pi;
    //double rotationAngle = pi - sunAngle;
    //double rotationSpeed = 2 * pi / (365.2422 * 24 * 60 * 60);

    ///
    /// Set the integrator settings, which are the same for every set of initial conditions
    ///
    //const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
    //        numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    //double minimumStepSizeS1 = -1.0E-5;
    //double maximumStepSizeS1 = -10;
    //double relativeErrorToleranceS1 = 1.0E-8;
    //double absoluteErrorToleranceS1 = 1.0E-8;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS2 =
            boost::make_shared< IntegratorSettings< > >
            (rungeKutta4, PropagationStartS2, StepSize );
    //boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
    //(rungeKuttaVariableStepSize, PropagationStartS1, InitialStepSize, coefficientSetS1,
    //minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );


    ///
    /// Calculation of the required constants to calculate the possible initial conditions -- performed in the CR3BP ///
    /// Not done for section two; identical to section one
    ///


    /// Calculation of the possible initial vectors in the CR3BP, using the equation:
    /// [x, y, x', y'] = alpha1 * e ^ (labda*t) * u1 + alpha2 * e ^ (-labda*t) * u2 + 2 * beta * e ^ (nu*t) * w1
    /// Alpha1, Alpha2, Beta and Time are the sliding variables, which create the set of initial conditions.

    // Set the number of samples on each dimension
    sampleSizeA1 = 1; //set to 0??
    sampleSizeA2 = 2;
    sampleSizeB = 2;
    sampleSizeT = 1;

    // Create the vectors of each of the sliding variables. For segment 2, alpha1 must be larger than zero and alpha2 must be zer0
    alpha1 = Eigen::VectorXd::LinSpaced(sampleSizeA1, -1E-300, 1E-300);
    alpha2 = Eigen::VectorXd::LinSpaced(sampleSizeA2, 0.0, 5E-10);
    beta = Eigen::VectorXd::LinSpaced(sampleSizeB, -2.5E-10, -1.0E-10);
    time = Eigen::VectorXd::LinSpaced(sampleSizeT, 0.0, 10000);

    // Create the two matrix that will hold the initial conditions in the frame used for propagation.
    // initialConditionsS1 saves only x, y, xdot and ydot, along with the used alpha and beta. This is used for later reference.
    // initialCartesianElementsS1 holds the full set of Cartesian Elements and is used for the propagation. It is overwritten each sequence.
    Eigen::Matrix4d initialConditionsS2;
    Eigen::Vector6d initialCartesianElementsS2;

    // Create the matrix that will hold the output trajectories for each compination of alpha and beta.
    Eigen::MatrixXd cartesianIntegrationResultS2;

    for(int alpha1Counter = 0; alpha1Counter <= alpha1.size(); alpha1Counter++)
    {
        std::cerr<<"S2, Alpha1 counter " + boost::lexical_cast< std::string >( alpha1Counter )<<std::endl;
        for(int alpha2Counter = 0; alpha2Counter <= alpha2.size(); alpha2Counter++)
        {
            std::cerr<<"S2, alpha2 counter " + boost::lexical_cast< std::string >( alpha2Counter )<<std::endl;
            //std::cerr<<"Alpha2+."<<std::endl;
            for(int betaCounter = 0; betaCounter <= beta.size(); betaCounter++)
            {
                std::cerr<<"S2, Beta counter " + boost::lexical_cast< std::string >( betaCounter )<<std::endl;
                //std::cerr<<"Beta+."<<std::endl;
                for(int timeCounter = 0; timeCounter <= time.size(); timeCounter++)
                {

                    //std::cerr<<"Time+."<<std::endl;
                    ///
                    /// for the given alpha1, alpha2, beta and time, set the initial conditions
                    ///

                    // initial conditions in the CR3BP
                    x =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[0]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[0]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[0] );
                    y =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[1]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[1]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[1] );
                    xdot =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[2]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[2]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[2] );
                    ydot =
                            alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[3]
                            + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[3]
                            + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[3] );

                    // scale the coordinates
                    x2 = sunDistance * (x + L2Distance -1.0 + earthSunMu );
                    y2 = sunDistance * y;
                    xdot2 = sunDistance * xdot;
                    ydot2 = sunDistance * ydot;

                    // add rotation and rotational speed to get initial conditions in the earth centered frame

                    initialConditionsS2(0) = x2 * std::cos(rotationAngle) + y2 * std::sin(rotationAngle) ;
                    initialConditionsS2(1) = y2 * std::cos(rotationAngle) - x2 * std::sin(rotationAngle) ;

                    initialConditionsS2(2) = xdot2 * std::cos(rotationAngle) + ydot2 * std::sin(rotationAngle) + rotationSpeed * (- std::sin(rotationAngle) * x2 + std::cos(rotationAngle) * y2) ;
                    initialConditionsS2(3) = ydot2 * std::cos(rotationAngle) - xdot2 * std::sin(rotationAngle) + rotationSpeed * (- std::cos(rotationAngle) * x2 - std::sin(rotationAngle) * y2) ;

                    //initial conditions in the earth centered frame, as used for the propagation
                    initialCartesianElementsS2[ xCartesianPositionIndex ] = initialConditionsS2(0);
                    initialCartesianElementsS2[ yCartesianPositionIndex ] = initialConditionsS2(1);
                    initialCartesianElementsS2[ zCartesianPositionIndex ] = 0.0;
                    initialCartesianElementsS2[ xCartesianVelocityIndex ] = initialConditionsS2(2);
                    initialCartesianElementsS2[ yCartesianVelocityIndex ] = initialConditionsS2(3);
                    initialCartesianElementsS2[ zCartesianVelocityIndex ] = 0.0;

                    systemInitialState.segment(18,6) = initialCartesianElementsS2;

                    ///
                    /// Set the propagator settings
                    ///

                    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS2 =
                            boost::make_shared< TranslationalStatePropagatorSettings< > >
                            (centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, PropagationEndsS2, cowell );

                    ///
                    /// Create simulation object and propagate dynamics.
                    ///
                    SingleArcDynamicsSimulator< > dynamicsSimulatorS2(
                                bodyMap, integratorSettingsS2, propagatorSettingsS2, true, false, false );
                    std::map< double, Eigen::VectorXd > tempIntegrationResultS2 = dynamicsSimulatorS2.getEquationsOfMotionNumericalSolution( );

                    // Write segment 2 propagation to file for later access.
                    input_output::writeDataMapToTextFile( tempIntegrationResultS2,
                                                          "WP1_Segment2_alpha1=" + boost::lexical_cast< std::string >( alpha1[alpha1Counter] )
                                                          + "_alpha2=" + boost::lexical_cast< std::string >( alpha2[alpha2Counter] )
                                                          + "_beta=" + boost::lexical_cast< std::string >( beta[betaCounter] )
                                                          + "_time=" + boost::lexical_cast< std::string >( time[timeCounter] )
                                                          + ".dat",
                                                          tudat_applications::getOutputPath( ),
                                                          "",
                                                          std::numeric_limits< double >::digits10,
                                                          std::numeric_limits< double >::digits10,
                                                          "," );

                }
            }
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////       PROPAGATE SEGMENT 3          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment 3 calculations."<<std::endl;

    ///
    /// Set the moment at which the segment propagation starts. Note that backward propagation is used, thus the propagation start time is the end time for the trajectory
    ///

    const double PropagationStartS3 = 0.0;
    const double PropagationEndsS3 = -0.1E7;
    StepSize = -500.0;

    /// Sun location remains unchanged

    ///
    /// Set the integrator settings, which are the same for every set of initial conditions
    ///
    //const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSetS1 =
    //        numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
    //double minimumStepSizeS1 = -1.0E-5;
    //double maximumStepSizeS1 = -10;
    //double relativeErrorToleranceS1 = 1.0E-8;
    //double absoluteErrorToleranceS1 = 1.0E-8;

    boost::shared_ptr< IntegratorSettings< > > integratorSettingsS3 =
            boost::make_shared< IntegratorSettings< > >
            (rungeKutta4, PropagationStartS3, StepSize );
    //boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
    //(rungeKuttaVariableStepSize, PropagationStartS1, InitialStepSize, coefficientSetS1,
    //minimumStepSizeS1, maximumStepSizeS1, relativeErrorToleranceS1, absoluteErrorToleranceS1 );


    ///
    /// Calculation of the required constants to calculate the possible initial conditions -- performed in the CR3BP ///
    ///
    //earthMass = 5.972E24;
    const double moonMass = 7.342E22;
    const double moonEarthMu = moonMass/(earthMass+moonMass);
    L2Distance = 1.0 + pow( moonMass / (3.0 * earthMass) , 1.0/3.0);

    rho = moonEarthMu * pow(L2Distance - 1.0 + moonEarthMu,-3) + (1.0 - moonEarthMu) * pow(L2Distance + moonEarthMu,-3);
    a = 2.0 * rho + 1.0;
    b = rho - 1.0;

    labda = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
    nu = pow( pow( pow(a,2) + 2.0 * a * b - 8.0 * a + pow(b,2) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

    //Eigen::Vector4d u1;
    u1[0] = 1;
    u1[1] = u1[0] * ( pow( labda , 2 ) - a ) / ( 2*labda );
    u1[2] = u1[0] * labda;
    u1[3] = u1[1] * labda;

    //Eigen::Vector4d u2;
    u2[0] = u1[0];
    u2[1] = - u1[1];
    u2[2] = - u1[2];
    u2[3] = u1[3];

    //Eigen::Vector4cd w1;
    w1[0] = 1.0;
    w1[1] = w1[0] * ( -pow( nu , 2.0 ) - a ) / (2 * nu * std::complex<double>(0,-1));
    w1[2] = w1[0] * nu * std::complex<double>(0,-1);
    w1[3] = w1[1] * nu * std::complex<double>(0,-1);

    //std::cerr<< "w1 =  " + boost::lexical_cast< std::string >( w1 ) <<std::endl;


    /// Calculation of the possible initial vectors in the CR3BP, using the equation:
    /// [x, y, x', y'] = alpha1 * e ^ (labda*t) * u1 + alpha2 * e ^ (-labda*t) * u2 + 2 * beta * e ^ (nu*t) * w1
    /// Alpha1, Alpha2, Beta and Time are the sliding variables, which create the set of initial conditions.

    // Set the number of samples on each dimension
    sampleSizeA1 = 2;
    sampleSizeA2 = 1; //set to 0??
    sampleSizeB = 2;
    sampleSizeT = 1;
    int sampleSizeM = 1;

    // Create the vectors of each of the sliding variables. For segment 1, alpha1 must be zero and alpha2 must be larger than 0
    alpha1 = Eigen::VectorXd::LinSpaced(sampleSizeA1, -1E-6, 1E-6);
    alpha2 = Eigen::VectorXd::LinSpaced(sampleSizeA2, -1E-300, 1E-300);
    beta = Eigen::VectorXd::LinSpaced(sampleSizeB, -1E-5, 0);
    time = Eigen::VectorXd::LinSpaced(sampleSizeT, 0.0, 10000);
    Eigen::VectorXd moonAngle = Eigen::VectorXd::LinSpaced(sampleSizeM, -1.5*pi, 1.0*pi);


    // Create the two matrix that will hold the initial conditions in the frame used for propagation.
    // initialConditionsS3 saves only x, y, xdot and ydot, along with the used alpha and beta. This is used for later reference.
    // initialCartesianElementsS3 holds the full set of Cartesian Elements and is used for the propagation. It is overwritten each sequence.
    Eigen::Matrix4d initialConditionsS3;
    Eigen::Vector6d initialCartesianElementsS3;

    // Create the matrix that will hold the output trajectories for each compination of alpha and beta.
    Eigen::MatrixXd cartesianIntegrationResultS3;

    for(int alpha1Counter = 0; alpha1Counter <= alpha1.size(); alpha1Counter++)
    {
        std::cerr<<"S3, Alpha1 counter " + boost::lexical_cast< std::string >( alpha1Counter )<<std::endl;
        for(int alpha2Counter = 0; alpha2Counter <= alpha2.size(); alpha2Counter++)
        {
            std::cerr<<"S3, Alpha2 counter " + boost::lexical_cast< std::string >( alpha2Counter )<<std::endl;
            //std::cerr<<"Alpha2+."<<std::endl;
            for(int betaCounter = 0; betaCounter <= beta.size(); betaCounter++)
            {
                std::cerr<<"S3, Beta counter " + boost::lexical_cast< std::string >( betaCounter )<<std::endl;
                //std::cerr<<"Beta+."<<std::endl;
                for(int timeCounter = 0; timeCounter <= time.size(); timeCounter++)
                {

                    for(int moonCounter = 0; moonCounter <= moonAngle.size(); moonCounter++)
                    {

                        /// Set the location of the Moon relative to the Earth at PropagationStartS3
                        MoonInitialStateKeplerian( trueAnomalyIndex ) = moonAngle[moonCounter];

                        double MoonDistance = MoonInitialStateKeplerian(semiMajorAxisIndex);
                        double MoonAngle = moonAngle[moonCounter];
                        double MoonRotationAngle = pi - MoonAngle;
                        double MoonRotationSpeed = 2 * pi / (27.321611 * 24 * 60 * 60);

                        ///
                        /// for the given alpha1, alpha2, beta and time, set the initial conditions
                        ///

                        // initial conditions in the Moon-Earth CR3BP
                        x =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[0]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[0]
                                + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[0] );
                        y =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[1]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[1]
                                + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[1] );
                        xdot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[2]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[2]
                                + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[2] );
                        ydot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[3]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[3]
                                + 2 * std::real( beta[betaCounter] * exp( nu * std::complex<double>(0,-1) * time[timeCounter]) * w1[3] );

                        // scale the coordinates
                        x2 = MoonDistance * (x + L2Distance + moonEarthMu );
                        y2 = MoonDistance * y;
                        xdot2 = MoonDistance * xdot;
                        ydot2 = MoonDistance * ydot;

                        // add rotation and rotational speed to get initial conditions in the earth centered frame

                        initialConditionsS3(0) = x2 * std::cos(MoonRotationAngle) + y2 * std::sin(MoonRotationAngle) ;
                        initialConditionsS3(1) = y2 * std::cos(MoonRotationAngle) - x2 * std::sin(MoonRotationAngle) ;

                        initialConditionsS3(2) = xdot2 * std::cos(MoonRotationAngle) + ydot2 * std::sin(MoonRotationAngle)
                                + MoonRotationSpeed * (- std::sin(MoonRotationAngle) * x2 + std::cos(MoonRotationAngle) * y2)
                                + rotationSpeed * (- std::sin(rotationAngle) * x2 + std::cos(rotationAngle) * y2) ;
                        initialConditionsS3(3) = ydot2 * std::cos(MoonRotationAngle) - xdot2 * std::sin(MoonRotationAngle)
                                + MoonRotationSpeed * (- std::cos(MoonRotationAngle) * x2 - std::sin(MoonRotationAngle) * y2)
                                + rotationSpeed * (- std::cos(rotationAngle) * x2 - std::sin(rotationAngle) * y2) ;

                        //initial conditions in the earth centered frame, as used for the propagation
                        initialCartesianElementsS3[ xCartesianPositionIndex ] = initialConditionsS3(0);
                        initialCartesianElementsS3[ yCartesianPositionIndex ] = initialConditionsS3(1);
                        initialCartesianElementsS3[ zCartesianPositionIndex ] = 0.0;
                        initialCartesianElementsS3[ xCartesianVelocityIndex ] = initialConditionsS3(2);
                        initialCartesianElementsS3[ yCartesianVelocityIndex ] = initialConditionsS3(3);
                        initialCartesianElementsS3[ zCartesianVelocityIndex ] = 0.0;

                        systemInitialState.segment(18,6) = initialCartesianElementsS3;

                        ///
                        /// Set the propagator settings
                        ///

                        boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsS3 =
                                boost::make_shared< TranslationalStatePropagatorSettings< > >
                                (centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, PropagationEndsS3, cowell );

                        ///
                        /// Create simulation object and propagate dynamics.
                        ///
                        SingleArcDynamicsSimulator< > dynamicsSimulatorS3(
                                    bodyMap, integratorSettingsS3, propagatorSettingsS3, true, false, false );
                        std::map< double, Eigen::VectorXd > tempIntegrationResultS3 = dynamicsSimulatorS3.getEquationsOfMotionNumericalSolution( );

                        // Write segment 1 propagation to file for later access.
                        input_output::writeDataMapToTextFile( tempIntegrationResultS3,
                                                              "WP1_Segment3_alpha1=" + boost::lexical_cast< std::string >( alpha1[alpha1Counter] )
                                                              + "_alpha2=" + boost::lexical_cast< std::string >( alpha2[alpha2Counter] )
                                                              + "_beta=" + boost::lexical_cast< std::string >( beta[betaCounter] )
                                                              + "_time=" + boost::lexical_cast< std::string >( time[timeCounter] )
                                                              + "_moonStart=" + boost::lexical_cast< std::string >( moonAngle[moonCounter] )
                                                              + ".dat",
                                                              tudat_applications::getOutputPath( ),
                                                              "",
                                                              std::numeric_limits< double >::digits10,
                                                              std::numeric_limits< double >::digits10,
                                                              "," );

                    }
                }
            }
        }
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        LINK SEGMENTS INTO TRAjECTORIES               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cerr<<"Started segment linking."<<std::endl;

    std::vector< boost::filesystem::path > FileNames = input_output::listAllFilesInDirectory(tudat_applications::getOutputPath( ));
    std::vector< std::string > S1Files;
    std::vector< std::string > S2Files;
    std::vector< std::string > S3Files;

    int i=0;
    int j=0;
    int k=0;

    for(int fileCounter = 0; fileCounter <= FileNames.size(); fileCounter++)
    {
        std::string currentFile = FileNames[fileCounter].string();

        if(std::stoi( currentFile.substr(11,1) ) == 1)
        {
            S1Files.push_back(FileNames[fileCounter].string());
            i=i+1;
        }
        else if(std::stoi( currentFile.substr(11,1) ) == 2)
        {
            S2Files.push_back(FileNames[fileCounter].string());
            j=j+1;
        }
        else if(std::stoi( currentFile.substr(11,1) ) == 3)
        {
            S3Files.push_back(FileNames[fileCounter].string());
            k=k+1;
        }

    }

    std::cerr<<"Sectors divided"<<std::endl;


    //Note to self: x@19, y@20
    int stop;
    int S1endtime;
    double S1endX;
    int S2starttime;
    double S2startX;
    int S2endtime;
    double S2endY;
    int S3starttime;
    double S3startY;

    for(int S1Counter = 0; S1Counter <= S1Files.size(); S1Counter++)
        // For each of the Segment1 trajectories, find the passage through y=0 near the Lagrange point if there is any.
    {
        std::string path = S1Files.at(S1Counter);
        Eigen::MatrixXd S1segment = input_output::readMatrixFromFile( path , ",");

        stop = 0;
        for(int S1time = 0; S1time <= S1segment.row(0).size(); S1time++)
        {
            if(stop == 0)
            {
                if(std::abs(S1segment(20,S1time)) <= 1E-5)
                {
                    if(S1segment(19,S1time) > 1E5 ){
                        stop = 1;
                        S1endtime = S1time;
                        S1endX = S1segment(19,S1time);
                    }
                }
            }
        }

        if(stop == 1)
            // if there is no passage through y=0 near the Lagrange point, this segment is not used for further calculations.
        {
            for(int S2Counter = 0; S2Counter <= S2Files.size(); S2Counter++)
                //For each of the Segment 2 trajectories, find the passage through y=0 near the Lagrange point if there is any. If there is, there is a S1-S2 link.

            {
                std::string path = S2Files.at(S2Counter);
                Eigen::MatrixXd S2segment = input_output::readMatrixFromFile( path, ",");

                stop = 0;
                int match12 = 0;
                for(int S2time = 0; S2time <= S2segment.row(0).size(); S2time++)
                {
                    if(stop == 0)
                    {
                        if(std::abs(S2segment(20,S2time)) <= 1E-5)
                        {
                            if(S2segment(19,S2time) > 1E5 ){
                                stop = 1;
                                S2starttime = S2time;
                                S2startX = S2segment(19,S2time);

                                if(std::abs(S1endX-S2startX) <= 100 )
                                {
                                    match12 = 1;
                                }
                            }
                        }
                    }
                }

                if(match12 == 1)
                    // if there is no match between segments 1 and 2, there is no use working on with this segment
                {

                    stop = 0;

                    for(int S2time = S2starttime; S2time<= S2segment.row(0).size(); S2time++)
                    {
                        if(stop == 0)
                        {
                            if(std::abs(S2segment(19,S2time)) <= 1E-5)
                            {
                                if(S2segment(20,S2time) > 0)
                                {
                                    stop = 1;
                                    S2endtime = S2time;
                                    S2endY = S2segment(20,S2time);
                                }
                            }
                        }
                    }

                    if(stop == 1)
                        // if there is no segment 2 ending in the desired area, there is no use working on with this segment
                    {
                        for(int S3Counter = 0; S3Counter <= S3Files.size(); S3Counter++)
                        {
                            std::string path = S3Files.at(S3Counter);
                            Eigen::MatrixXd S3segment = input_output::readMatrixFromFile( path , ",");

                            stop = 0;
                            int match23 = 0;

                            for(int S3time = 0; S3time <= S3segment.row(0).size(); S3time++)
                            {
                                if(stop == 0)
                                {
                                    if(std::abs(S3segment(19,S3time)) <= 1E-5)
                                    {
                                        if(S3segment(20,S3time) > 1E5 ){
                                            stop = 1;
                                            S3starttime = S3time;
                                            S3startY = S3segment(20,S3time);

                                            if(std::abs(S2endY-S3startY) <= 100 )
                                            {
                                                match23 = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


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


