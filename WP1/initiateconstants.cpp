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
const long double MoonGravitationalParameter = SunGravitationalParameter / ( 328900.56 * ( 1.0 + 81.30059 ) );
const long double SSBGravitationalParameter = SunGravitationalParameter + EarthGravitationalParameter;


const long double secInYear = tudat::physical_constants::SIDEREAL_YEAR;
const long double secInMonth = 2372380.0; //this is not the right number, but I made manual changes everywhere so I can't take it out anymore :(

const long double muSE = EarthGravitationalParameter / (SunGravitationalParameter + EarthGravitationalParameter);
const long double muEM = MoonGravitationalParameter / (EarthGravitationalParameter + MoonGravitationalParameter);

//Astronomical Unit in meters
//const long double Re = pow( SunGravitationalParameter * pow( secInYear/(2.0*pi) , 2.0 ), 1.0/3.0 );//tudat::physical_constants::ASTRONOMICAL_UNIT;
//const long double Rs = pow( EarthGravitationalParameter * pow( secInYear/(2.0*pi) , 2.0 ), 1.0/3.0 );//tudat::physical_constants::ASTRONOMICAL_UNIT;
const long double AU = pow( (SunGravitationalParameter + EarthGravitationalParameter ) * pow( secInYear/(2.0*pi) , 2.0 ), 1.0/3.0 );//Re + Rs;
const long double MoonEarthDistance = pow( (EarthGravitationalParameter + MoonGravitationalParameter ) * pow( secInMonth/(2.0*pi) , 2.0), 1.0/3.0 ); //0.002570 * AU; //SMAD p962;

 //MoonEarthDistance * 2.0 * pi / 1023.0;


bool InitiateConstants()
{
    return 1;
}


