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

    int nrofpoints = 1000;

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

    std::string SCFile = tudat_applications::getOutputPath() + "/WP1_Orbit_" + IDname + ".dat";

    std::ofstream out(SCFile);
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to SCFile.

    int counterlimit = integrationResult.size()/nrofpoints;
    int SCcounter = counterlimit;


    for(typename std::map< double, Eigen::VectorXd >::iterator it = integrationResult.begin();
        it != integrationResult.end();
        it++)
    {
        if(SCcounter == counterlimit)
        {
            Eigen::VectorXd State = it->second;
            std::cout << it->first << " , "
                      << State(0) << " , "
                      << State(1) << " , "
                      << State(2) << " , "
                      << State(3) << " , "
                      << State(4) << " , "
                      << State(5) << "\n";
            SCcounter = 0;
        }
        SCcounter++;
    }


    //write planet data to file
    if(addSun){
        std::string SunFile = tudat_applications::getOutputPath() + "/WP1_SunOrbit_" + IDname + ".dat";

        std::ofstream out(SunFile);
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to SCFile.

        int Counter = counterlimit;

        for(typename std::map< double, Eigen::VectorXd >::iterator it = sunBarycentricStates.begin();
            it != sunBarycentricStates.end();
            it++)
        {
            if(Counter == counterlimit)
            {
                Eigen::VectorXd State = it->second;
                std::cout << it->first << " , "
                          << State(0) << " , "
                          << State(1) << " , "
                          << State(2) << " , "
                          << State(3) << " , "
                          << State(4) << " , "
                          << State(5) << "\n";
                Counter = 0;
            }
            Counter++;
        }
    }

    if(addEarth){
        std::string EarthFile = tudat_applications::getOutputPath() + "/WP1_EarthOrbit_" + IDname + ".dat";

        std::ofstream out(EarthFile);
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to SCFile.

        int Counter = counterlimit;

        for(typename std::map< double, Eigen::VectorXd >::iterator it = earthBarycentricStates.begin();
            it != earthBarycentricStates.end();
            it++)
        {
            if(Counter == counterlimit)
            {
                Eigen::VectorXd State = it->second;
                std::cout << it->first << " , "
                          << State(0) << " , "
                          << State(1) << " , "
                          << State(2) << " , "
                          << State(3) << " , "
                          << State(4) << " , "
                          << State(5) << "\n";
                Counter = 0;
            }
            Counter++;
        }
    }
    if(addMoon){
        std::string MoonFile = tudat_applications::getOutputPath() + "/WP1_MoonOrbit_" + IDname + ".dat";

        std::ofstream out(MoonFile);
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to SCFile.

        int Counter = counterlimit;

        for(typename std::map< double, Eigen::VectorXd >::iterator it = moonBarycentricStates.begin();
            it != moonBarycentricStates.end();
            it++)
        {
            if(Counter == counterlimit)
            {
                Eigen::VectorXd State = it->second;
                std::cout << it->first << " , "
                          << State(0) << " , "
                          << State(1) << " , "
                          << State(2) << " , "
                          << State(3) << " , "
                          << State(4) << " , "
                          << State(5) << "\n";
                Counter = 0;
            }
            Counter++;
        }
    }

    std::cout.rdbuf(coutbuf); //reset to standard output again

}
