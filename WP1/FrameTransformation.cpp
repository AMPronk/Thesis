#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>
#include <tuple>
#include "FrameTransformation.h"

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


Eigen::Vector6d FrameTransformationSE(long double x, long double y, long double xdot, long double ydot)
{
    //define the distance of the L2 point from the Earth
    Eigen::Vector6d L2location = InitialStateSEL2();
    long double L2Distance = L2location[0];

    //define the speed at which the system rotates: One full rotation per year
    long double rotationSpeed = 2.0 * pi / secInYear;
    //define the angle at which the L2 point is currently placed, compared to the x-axis.
    //TODO, shouldn't be hardcode
    long double rotationAngle = 0.0;

    // scale the coordinates for distance
    long double x2 = (AU * x) + L2Distance;
    long double y2 = AU * y;
    // scale the speeds for distance and time
    long double xdot2 = AU * rotationSpeed * xdot ;
    long double ydot2 = AU * rotationSpeed * ydot;


    // add rotation and rotational speed to get initial conditions in the earth centered frame

    long double x3 = x2 * std::cos(rotationAngle) + y2 * std::sin(rotationAngle) ;
    long double y3 = y2 * std::cos(rotationAngle) - x2 * std::sin(rotationAngle) ;

    long double xdot3 = xdot2 * std::cos(rotationAngle) + ydot2 * std::sin(rotationAngle) -
            rotationSpeed * (- std::sin(rotationAngle) * x2 + std::cos(rotationAngle) * y2) ;
    long double ydot3 = ydot2 * std::cos(rotationAngle) - xdot2 * std::sin(rotationAngle) -
            rotationSpeed * (- std::cos(rotationAngle) * x2 - std::sin(rotationAngle) * y2) ;


    //initial conditions in the earth centered frame, as used for the propagation
    Eigen::Vector6d CartesianElements;
    CartesianElements[ xCartesianPositionIndex ] = x3;
    CartesianElements[ yCartesianPositionIndex ] = y3;
    CartesianElements[ zCartesianPositionIndex ] = 0.0;
    CartesianElements[ xCartesianVelocityIndex ] = xdot3;
    CartesianElements[ yCartesianVelocityIndex ] = ydot3;
    CartesianElements[ zCartesianVelocityIndex ] = 0.0;

    return CartesianElements;

}

Eigen::Vector6d FrameTransformationEM(long double x, long double y, long double xdot, long double ydot)
{
    //define the distance of the L2 point from the SSB
    Eigen::Vector6d L2location = InitialStateEML2();
    long double L2Distance = abs(L2location[0]);
    long double L2Speed = L2location[4];

    //define the speed at which the system rotates: set rotational time for the moon
    long double rotationSpeed = L2Speed / L2Distance; //2.0 * pi / secInMonth;
    //define the angle at which the L2 point is currently placed, compared to the x-axis.
    //TODO, shouldn't be hardcode
    long double rotationAngle = 0.0;

//    std::cerr<<"x = " + boost::lexical_cast< std::string >(x) + ". y = " + boost::lexical_cast< std::string >( y )
//               + ". xdot = " + boost::lexical_cast< std::string >( xdot ) + ". ydot = " + boost::lexical_cast< std::string >( ydot )<<std::endl;


    // scale the coordinates for distance
    long double x2 = (MoonEarthDistance * x) + L2Distance;
    long double y2 = MoonEarthDistance * y;
    // scale the speeds for distance and time
    long double xdot2 = MoonEarthDistance * rotationSpeed * xdot;
    long double ydot2 = MoonEarthDistance * rotationSpeed * ydot;

//    std::cerr<<"x2 = " + boost::lexical_cast< std::string >(x2) + ". y2 = " + boost::lexical_cast< std::string >( y2 )
//               + ". xdot2 = " + boost::lexical_cast< std::string >( xdot2 ) + ". ydot2 = " + boost::lexical_cast< std::string >( ydot2 )<<std::endl;


    // add rotation and rotational speed to get initial conditions in the earth centered frame

    long double x3 = x2;
    long double y3 = y2;

    long double xdot3 = xdot2 - rotationSpeed * y2 ;
    long double ydot3 = ydot2 + rotationSpeed * x2 ;

//    std::cerr<<"x3 = " + boost::lexical_cast< std::string >(x3) + ". y3 = " + boost::lexical_cast< std::string >( y3 )
//               + ". xdot3 = " + boost::lexical_cast< std::string >( xdot3 ) + ". ydot3 = " + boost::lexical_cast< std::string >( ydot3 )<<std::endl;


    //initial conditions in the earth centered frame, as used for the propagation
    Eigen::Vector6d CartesianElements;
    CartesianElements[ xCartesianPositionIndex ] = x3;
    CartesianElements[ yCartesianPositionIndex ] = y3;
    CartesianElements[ zCartesianPositionIndex ] = 0.0;
    CartesianElements[ xCartesianVelocityIndex ] = xdot3;
    CartesianElements[ yCartesianVelocityIndex ] = ydot3;
    CartesianElements[ zCartesianVelocityIndex ] = 0.0;

    return CartesianElements;

}
