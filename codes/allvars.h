/**************************************************************************************************
			      STRUCTURES
**************************************************************************************************/
//Gadget header structure
struct gadget_head{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
    };

//Particle structure
struct part{
    //General properties
    float pos[3];
    float vel[3];
    int id;
    //Halo properties
    int id_halo;
    //Gas properties
    float mass;
    float energy;
    float rho;
    float pressure;
    float temperature;
    //Simulation properties
    float z;
    float t;
    };

    
/**************************************************************************************************
			      MACROS AND GLOBAL VARIABLES
**************************************************************************************************/
//Numerical and internal macros
#define NMAX1		1000
#define X		0
#define Y		1
#define Z		2

//Physical macros
#define GAMMA		5/3.0
#define K_B		1.3806488e-23		// [ J K^{-1} ]
#define M_ATO		1.660538921e-27		// [ Kg ]
#define MU		1			// Hydrogen gas
#define U_RHO		6.73928e-19		// [ 10e10 h-1 Msun/(h-1 kpc)^3 ] -> [ Kg/m^3 ]
#define U_ENE		1e6			// [ (km/sec)^2 ] -> [ (m/sec)^2 ]

//Particle arrays
struct part *Part;
//Gadget head
struct gadget_head Gheader;
//Number of particles in current snap
int Npart_snap;
//Number of particles in current snap with variable masses
int Npart_snap_mass;


/**************************************************************************************************
			      HEADERS
**************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <proto.h>