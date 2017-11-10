#include "writetofile.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>
#include <tuple>

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

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getOutputPath( )
{
    return "D:/FILES/Documents/Thesis/PropagationResults";
}

}

void WriteToFile(std::map< double, Eigen::VectorXd > integrationResult, NamedBodyMap bodyMap, std::string IDname, bool addSun, bool addEarth, bool addMoon)
{
    std::map< double, Eigen::VectorXd > sunBarycentricStates;
    std::map< double, Eigen::VectorXd > earthBarycentricStates;
    std::map< double, Eigen::VectorXd > moonBarycentricStates;

    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
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

    //Write SC data to file
    input_output::writeDataMapToTextFile( integrationResult,
                                          "WP1_Orbit_" + IDname
                                          + ".dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    //write planet data to file
    if(addSun){
        input_output::writeDataMapToTextFile( sunBarycentricStates,
                                              "WP1_SunOrbit_" + IDname
                                              + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }
    if(addEarth){
        input_output::writeDataMapToTextFile( earthBarycentricStates,
                                              "WP1_EarthOrbit_" + IDname
                                              + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }
    if(addMoon){
        input_output::writeDataMapToTextFile( moonBarycentricStates,
                                              "WP1_MoonOrbit_" + IDname
                                              + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }

}
