#include "initiateconstants.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////     SETTING OF NECESSARY CONSTANTS       //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const long double pi = tudat::mathematical_constants::LONG_PI ;
// Gravitational Parameters from O&CD&M, p835
const long double SunGravitationalParameter = 1.32712438E20;
const long double EarthGravitationalParameter = 3.98600441E14;
const long double MoonGravitationalParameter = 4.90279898E12;

//Astronomical Unit in meters
const long double AU = tudat::physical_constants::ASTRONOMICAL_UNIT;

const long double MoonEarthDistance = 0.002570 * AU; //SMAD p962;

const long double secInYear = tudat::physical_constants::JULIAN_YEAR;
const long double secInMonth = MoonEarthDistance * 2.0 * pi / 1023.0;

const long double muSE = EarthGravitationalParameter / (SunGravitationalParameter + EarthGravitationalParameter);
const long double muEM = MoonGravitationalParameter / (EarthGravitationalParameter + MoonGravitationalParameter);

bool InitiateConstants()
{
    return 1;
}
