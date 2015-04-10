#include <allvars.h>

/**************************************************************************************************
 NAME:	     pressure
 FUNCTION:   Calculate pressure of a gas particle
 INPUTS:     density [ 10e10 h-1 Msun/(h-1 kpc)^3 ], internal energy [ (km/sec)^2 ]
 RETURN:     pressure [ Pa/(h-1)^2 ]
**************************************************************************************************/
float pressure( float rho,
		float energy )
{
    float P;
    
    //Pressure of an ideal gas
    P = (GAMMA - 1)*(rho*U_RHO)*(energy*U_ENE);
    
    return P;
}


/**************************************************************************************************
 NAME:	     temperature
 FUNCTION:   Calculate pressure of a gas particle
 INPUTS:     internal energy [ (km/sec)^2 ]
 RETURN:     temperature [ K ]
**************************************************************************************************/
float temperature( float energy )
{
    float T;
    
    //Temperature of an ideal gas
    T = ((GAMMA - 1)/K_B)*(MU*M_ATO)*(energy*U_ENE);
  
    return T;
}