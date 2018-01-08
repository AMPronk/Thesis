#include "initiateconstants.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <complex>
#include <iostream>
#include <tuple>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////     SETTING OF NECESSARY CONSTANTS       //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const long double pi = tudat::mathematical_constants::LONG_PI ;

// Gravitational Parameters from tudat
const long double SunGravitationalParameter = 1.32712440018E20;
const long double EarthGravitationalParameter = 3.986004418E14;
const long double MoonGravitationalParameter = 4.9048695E12;

const long double SSBGravitationalParameter = SunGravitationalParameter + EarthGravitationalParameter;


const long double secInYear = tudat::physical_constants::SIDEREAL_YEAR;
const long double secInMonth = 27.322 * tudat::physical_constants::SIDEREAL_DAY;

const long double muSE = EarthGravitationalParameter / (SunGravitationalParameter + EarthGravitationalParameter);
const long double muEM = MoonGravitationalParameter / (EarthGravitationalParameter + MoonGravitationalParameter);

//Astronomical Unit in meters
const long double AU = pow( (SunGravitationalParameter + EarthGravitationalParameter ) * pow( secInYear/(2.0*pi) , 2.0 ), 1.0/3.0 );//Re + Rs;
const long double MoonEarthDistance = pow( (EarthGravitationalParameter + MoonGravitationalParameter ) * pow( secInMonth/(2.0*pi) , 2.0), 1.0/3.0 ); //0.002570 * AU; //SMAD p962;



bool InitiateConstants()
{
    return 1;
}


