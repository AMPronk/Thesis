#include "createinitialconditions.h"
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

Eigen::MatrixXd CreateInitialConditionsSE(int orbitType, int sampleSize, long double energyLowerBound, long double energyHigherBound)
{
    ///
    /// CreateInitialConditionsSE creates a set of initial conditions in the Sun Earth system, determined by the input variables
    /// int orbitType determines which type of orbit is used (definition below)
    /// int sampleSize determines how large the initial sample of initial conditions is
    /// doubles energyLowerBound and energyHigherBound determine the allowed regions for the jacobi's constant of the initial conditions
    ///

    ///orbit type definition
    // orbitType = 0 --> periodic orbit.                                                        alpha1 = 0 & alpha2 = 0.
    // orbitType = 1 --> asymptotic orbit, towards L2.                                          alpha1 = 0 & alpha2 =/ 0.
    // orbitType = 2 --> asymptotic orbit, away from L2.                                        alpha1 =/ 0 & alpha2 = 0.
    // orbitType = 3 --> non-transit orbit, inside of hillsphere.                               alpha1 < 0 & alpha2 < 0.
    // orbitType = 4 --> non-transit orbit, outside of hillsphere.                              alpha1 > 0 & alpha2 > 0.
    // orbitType = 5 --> transit orbit, into the hillsphere                                     alpha1 < 0 & alpha2 > 0.
    // orbitType = 6 --> transit orbit, out of the hillsphere                                   alpha1 > 0 & alpha2 < 0.

    /// create possible sets for alpha1, alpha2, beta, time

    //Create required empty vectors
    Eigen::MatrixXd InitialConditions(1,9);
    Eigen::VectorXd alpha1(1);
    alpha1.setZero();
    Eigen::VectorXd alpha2(1);
    alpha2.setZero();
    Eigen::VectorXd beta;
    Eigen::VectorXd time;

    //Set the limit value for the variables
    double limitvalue = 5.0E-2;

    //Create vector with the initial values of beta
    long double LowerBeta = -limitvalue;
    long double HigherBeta = limitvalue;
    beta = Eigen::VectorXd::LinSpaced(4*sampleSize, LowerBeta, HigherBeta);

    //Create vector with the initial values of t
    long double LowerTime = 0.0;
    long double HigherTime = 3.0;
    int TimeCount = sampleSize;
    time = Eigen::VectorXd::LinSpaced(TimeCount, LowerTime, HigherTime);

    //Per orbit type, create vectors with the initial values of alpha1 and alpha2
    if(orbitType == 0)
    {
        std::cerr<<"Periodic Orbit."<<std::endl;
        alpha1(0) = 0.0;
        alpha2(0) = 0.0;
    }
    else if(orbitType == 1)
    {
        std::cerr<<"Asymptotic, towards L2, backwards propagation. "<<std::endl;
        alpha1[0] = 0.0;

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 2)
    {
        std::cerr<<"Asymptotic, away from L2, positive propagation. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        alpha2[0] = 0.0;
    }
    else if(orbitType == 3)
    {
        std::cerr<<"non-transit orbit, inside hillsphere. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = 0.0;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = 0.0;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 4)
    {
        std::cerr<<"non-transit orbit, outside hillsphere. "<<std::endl;
        long double LowerAlpha1 = 0.0;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = 0.0;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 5)
    {
        std::cerr<<"transit orbit, inside hillsphere. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = 0.0;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = 0.0;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 6)
    {
        std::cerr<<"transit orbit, outside hillsphere. "<<std::endl;
        long double LowerAlpha1 = 0.0;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = 0.0;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else
    {
        std::cerr<<"Invalid orbit type. empty matrix returned."<<std::endl;
        return InitialConditions;
    }

    /// Create the eigenvectors required for transformation of axis

    Eigen::Vector6d L2location = InitialStateSEL2();
    long double L2Distance = L2location[0]/AU + 1.0 - muSE;

    long double rho = muSE * pow(L2Distance - 1.0 + muSE,-3.0) + (1.0 - muSE) * pow(L2Distance + muSE,-3.0);
    long double a = 2.0 * rho + 1.0;
    long double b = rho - 1.0;

    long double labda = pow( pow( pow(a,2.0) + 2.0 * a * b - 8.0 * a + pow(b,2.0) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
    long double nu = pow( pow( pow(a,2.0) + 2.0 * a * b - 8.0 * a + pow(b,2.0) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

    Eigen::Vector4d u1;
    u1[0] = 1.0;
    u1[1] = u1[0] * ( pow( labda , 2.0 ) - a ) / ( 2.0 * labda );
    u1[2] = u1[0] * labda;
    u1[3] = u1[1] * labda;

    Eigen::Vector4d u2;
    u2[0] = u1[0];
    u2[1] = - u1[1];
    u2[2] = - u1[2];
    u2[3] = u1[3];

    Eigen::Vector4cd w1;
    w1[0] = 1.0;
    w1[1] = w1[0] * ( -pow( static_cast<double>(nu) , 2.0 ) - static_cast<double>(a) ) / (2.0 * static_cast<double>(nu) * std::complex<double>(0.0,-1.0));
    w1[2] = w1[0] * static_cast<double>(nu) * std::complex<double>(0.0,-1.0);
    w1[3] = w1[1] * static_cast<double>(nu) * std::complex<double>(0.0,-1.0);

    long double x;
    long double y;
    long double xdot;
    long double ydot;
    long double C;
    int resultcounter = 0;

    /// transfor axis to get initial conditions in the xy plane

    for(int alpha1Counter = 0; alpha1Counter < alpha1.size(); alpha1Counter++)
        {
            for(int alpha2Counter = 0; alpha2Counter < alpha2.size(); alpha2Counter++)
            {
                for(int betaCounter = 0; betaCounter < beta.size(); betaCounter++)
                {
                    for(int timeCounter = 0; timeCounter < time.size(); timeCounter++)
                    {
                        // for the given alpha1, alpha2, beta and time, set the initial conditions

                        x =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[0]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[0]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[0] );
                        y =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[1]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[1]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[1] );

                        xdot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[2]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[2]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[2] );
                        ydot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[3]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[3]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[3] );

                        C = pow( x + L2Distance , 2.0 ) + pow( y , 2.0 ) +
                                2.0 * (1.0 - muSE)/ pow( pow( x + L2Distance + muSE , 2.0 ) + pow( y , 2.0 ) , 0.5 ) +
                                2.0 * muSE / pow( pow( x + L2Distance - 1.0 + muSE , 2.0 ) + pow( y , 2.0 ) , 0.5 ) -
                                pow( pow(xdot,2.0) + pow(ydot,2.0) , 0.5 );

                        // Check if the found initial conditions have acceptable energy level. If yes, save the initial conditions
                        if((energyLowerBound <= C) && (C <= energyHigherBound))
                        {
                            InitialConditions.conservativeResize(InitialConditions.rows()+1,InitialConditions.cols());

                            InitialConditions(InitialConditions.rows()-1,0) = x;
                            InitialConditions(InitialConditions.rows()-1,1) = y;
                            InitialConditions(InitialConditions.rows()-1,2) = xdot;
                            InitialConditions(InitialConditions.rows()-1,3) = ydot;
                            InitialConditions(InitialConditions.rows()-1,4) = C;
                            InitialConditions(InitialConditions.rows()-1,5) = alpha1[alpha1Counter];
                            InitialConditions(InitialConditions.rows()-1,6) = alpha2[alpha2Counter];
                            InitialConditions(InitialConditions.rows()-1,7) = beta[betaCounter];
                            InitialConditions(InitialConditions.rows()-1,8) = time[timeCounter];
                            resultcounter++;
                        }

                    }
                }
            }
    }

    std::cerr<<"Finished. Found " +
               boost::lexical_cast< std::string >( resultcounter )
               + " accepted SE initial conditions. "<<std::endl;

    return InitialConditions;
}

Eigen::MatrixXd CreateInitialConditionsEM(int orbitType, int sampleSize, long double energyLowerBound, long double energyHigherBound)
{
    ///
    /// CreateInitialConditionsEM creates a set of initial conditions in the Earth Moon system, determined by the input variables
    /// int orbitType determines which type of orbit is used (definition below)
    /// int sampleSize determines how large the initial sample of initial conditions is
    /// doubles energyLowerBound and energyHigherBound determine the allowed regions for the jacobi's constant of the initial conditions
    ///

    ///orbit type definition
    // orbitType = 0 --> periodic orbit.                                                        alpha1 = 0 & alpha2 = 0.
    // orbitType = 1 --> asymptotic orbit, towards L2.                                          alpha1 = 0 & alpha2 =/ 0.
    // orbitType = 2 --> asymptotic orbit, away from L2.                                        alpha1 =/ 0 & alpha2 = 0.
    // orbitType = 3 --> non-transit orbit, inside of hillsphere.                               alpha1 < 0 & alpha2 < 0.
    // orbitType = 4 --> non-transit orbit, outside of hillsphere.                              alpha1 > 0 & alpha2 > 0.
    // orbitType = 5 --> transit orbit, into the hillsphere                                     alpha1 < 0 & alpha2 > 0.
    // orbitType = 6 --> transit orbit, out of the hillsphere                                   alpha1 > 0 & alpha2 < 0.

    /// create possible sets for alpha1, alpha2, beta, time

    //Create required empty vectors
    Eigen::MatrixXd InitialConditions(1,9);
    Eigen::VectorXd alpha1(1);
    alpha1.setZero();
    Eigen::VectorXd alpha2(1);
    alpha2.setZero();
    Eigen::VectorXd beta;
    Eigen::VectorXd time;

    //Set the limit value for the variables
    double limitvalue = 1.0E-2;

    //Create vector with the initial values of beta
    long double LowerBeta = -limitvalue;
    long double HigherBeta = limitvalue;
    beta = Eigen::VectorXd::LinSpaced(4*sampleSize, LowerBeta, HigherBeta);

    //Create vector with the initial values of t
    long double LowerTime = 0.0;
    long double HigherTime = 3.0;
    int TimeCount = sampleSize; //todo, that's just ugly man.
    time = Eigen::VectorXd::LinSpaced(TimeCount, LowerTime, HigherTime);

    //Per orbit type, create vectors with the initial values of alpha1 and alpha2
    if(orbitType == 0)
    {
        std::cerr<<"Periodic Orbit."<<std::endl;
        alpha1(0) = 0.0;
        alpha2(0) = 0.0;
    }
    else if(orbitType == 1)
    {
        std::cerr<<"Asymptotic, towards L2, backwards propagation. "<<std::endl;
        alpha1[0] = 0.0;

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 2)
    {
        std::cerr<<"Asymptotic, away from L2, positive propagation. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        alpha2[0] = 0.0;
    }
    else if(orbitType == 3)
    {
        std::cerr<<"non-transit orbit, inside hillsphere. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = 0.0;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = 0.0;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 4)
    {
        std::cerr<<"non-transit orbit, outside hillsphere. "<<std::endl;
        long double LowerAlpha1 = 0.0;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = 0.0;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 5)
    {
        std::cerr<<"transit orbit, inside hillsphere. "<<std::endl;
        long double LowerAlpha1 = -limitvalue;
        long double HigherAlpha1 = 0.0;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = 0.0;
        long double HigherAlpha2 = limitvalue;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else if(orbitType == 6)
    {
        std::cerr<<"transit orbit, outside hillsphere. "<<std::endl;
        long double LowerAlpha1 = 0.0;
        long double HigherAlpha1 = limitvalue;
        int Alpha1Count = sampleSize;
        alpha1.resize(Alpha1Count);
        alpha1 = Eigen::VectorXd::LinSpaced(Alpha1Count, LowerAlpha1, HigherAlpha1);

        long double LowerAlpha2 = -limitvalue;
        long double HigherAlpha2 = 0.0;
        int Alpha2Count = sampleSize;
        alpha2.resize(Alpha2Count);
        alpha2 = Eigen::VectorXd::LinSpaced(Alpha2Count, LowerAlpha2, HigherAlpha2);
    }
    else
    {
        std::cerr<<"Invalid orbit type. empty matrix returned."<<std::endl;
        return InitialConditions;
    }

    /// Create the eigenvectors required for the transformation of axis

    Eigen::Vector6d L2location = InitialStateEML2();
    long double L2Distance = abs(L2location[0]) / MoonEarthDistance + 1.0 - muEM;

    long double rho = muEM * pow(L2Distance - 1.0 + muEM, -3.0) + (1.0 - muEM) * pow(L2Distance + muEM,-3.0);
    long double a = 2.0 * rho + 1.0;
    long double b = rho - 1.0;

    long double labda = pow( pow( pow(a,2.0) + 2.0 * a * b - 8.0 * a + pow(b,2.0) + 8.0 * b + 16.0,0.5) + a - b - 4.0,0.5) / pow( 2.0, 0.5 );
    long double nu = pow( pow( pow(a,2.0) + 2.0 * a * b - 8.0 * a + pow(b,2.0) + 8.0 * b + 16.0,0.5) - a + b + 4.0,0.5) / pow( 2.0, 0.5 );

    Eigen::Vector4d u1;
    u1[0] = 1.0;
    u1[1] = u1[0] * ( pow( labda , 2.0 ) - a ) / ( 2.0 * labda );
    u1[2] = u1[0] * labda;
    u1[3] = u1[1] * labda;

    Eigen::Vector4d u2;
    u2[0] = u1[0];
    u2[1] = - u1[1];
    u2[2] = - u1[2];
    u2[3] = u1[3];

    Eigen::Vector4cd w1;
    w1[0] = 1.0;
    w1[1] = w1[0] * ( -pow( static_cast<double>(nu) , 2.0 ) - static_cast<double>(a) ) / (2.0 * static_cast<double>(nu) * std::complex<double>(0.0,-1.0));
    w1[2] = w1[0] * static_cast<double>(nu) * std::complex<double>(0.0,-1.0);
    w1[3] = w1[1] * static_cast<double>(nu) * std::complex<double>(0.0,-1.0);

    long double x;
    long double y;
    long double xdot;
    long double ydot;
    long double C;
    int resultcounter = 0;

    /// transfor axis to get initial conditions in the xy plane

    for(int alpha1Counter = 0; alpha1Counter < alpha1.size(); alpha1Counter++)
        {
            for(int alpha2Counter = 0; alpha2Counter < alpha2.size(); alpha2Counter++)
            {
                for(int betaCounter = 0; betaCounter < beta.size(); betaCounter++)
                {
                    for(int timeCounter = 0; timeCounter < time.size(); timeCounter++)
                    {
                        // for the given alpha1, alpha2, beta and time, set the initial conditions

                        x =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[0]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[0]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[0] );
                        y =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[1]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[1]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[1] );

                        xdot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[2]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[2]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[2] );
                        ydot =
                                alpha1[alpha1Counter] * exp( - labda * time[timeCounter] ) * u1[3]
                                + alpha2[alpha2Counter] * exp( - labda * time[timeCounter] ) * u2[3]
                                + 2.0 * std::real( beta[betaCounter] * exp( (double)nu * std::complex<double>(0.0,-1.0) * time[timeCounter]) * w1[3] );

                        C = pow( x + L2Distance , 2.0 ) + pow( y , 2.0 ) +
                                2.0 * (1.0 - muEM)/ pow( pow( x + L2Distance + muEM , 2.0 ) + pow( y , 2.0 ) , 0.5 ) +
                                2.0 * muEM / pow( pow( x + L2Distance - 1.0 + muEM , 2.0 ) + pow( y , 2.0 ) , 0.5 ) -
                                pow( pow(xdot,2.0) + pow(ydot,2.0) , 0.5 );

                        // Check if the found initial conditions have acceptable energy level. If yes, save the initial conditions
                        if((energyLowerBound <= C) && (C <= energyHigherBound))
                        {
                            InitialConditions.conservativeResize(InitialConditions.rows()+1,InitialConditions.cols());

                            InitialConditions(InitialConditions.rows()-1,0) = x;
                            InitialConditions(InitialConditions.rows()-1,1) = y;
                            InitialConditions(InitialConditions.rows()-1,2) = xdot;
                            InitialConditions(InitialConditions.rows()-1,3) = ydot;
                            InitialConditions(InitialConditions.rows()-1,4) = C;
                            InitialConditions(InitialConditions.rows()-1,5) = alpha1[alpha1Counter];
                            InitialConditions(InitialConditions.rows()-1,6) = alpha2[alpha2Counter];
                            InitialConditions(InitialConditions.rows()-1,7) = beta[betaCounter];
                            InitialConditions(InitialConditions.rows()-1,8) = time[timeCounter];
                            resultcounter++;
                        }

                    }
                }
            }
    }

    std::cerr<<"Finished. Found " +
               boost::lexical_cast< std::string >( resultcounter )
               + " accepted EM initial conditions. "<<std::endl;

    return InitialConditions;
}
