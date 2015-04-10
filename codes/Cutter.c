#include <allvars.h>
/*Usage:	
 	   ./Cutter.out 
	   *	<snapbase>
 	   *	<number of files per snap>
 	   *	<axis> 
 	   *	<slide>				[Mpc h^-1]
 	   *	<dx> 				[Mpc h^-1]
 	   *	<density of data sampling>
 	   *	<output filename> 
 	   *	<0-all or 1-gas or 2-DM>
*/

int main( int argc, char *argv[] )
{
    int i, file_snap, sampling, type;
    long int n_slide;

    int axis;
    float slide, dx;
    char snapbase[NMAX1], output[NMAX1];

    //Snapbase filename
    sprintf( snapbase, "%s", argv[1] );
    //Number of files per snap
    file_snap = atoi( argv[2] );
    //Axis
    axis = atoi( argv[3] );
    //Slide
    slide = atof( argv[4] );
    //Dx
    dx = atof( argv[5] );
    //Density of data sampling
    sampling = atoi( argv[6] );
    //Output filename
    sprintf( output, "%s", argv[7] );
    //Type of data
    type = atoi( argv[8] );

    //Reading data from Gadget file
    read_snap( snapbase, file_snap, type );

    //Writing data
    ascii_data_slide( Part, output, sampling, type, axis, slide, dx );

    return 0;
}