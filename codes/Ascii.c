#include <allvars.h>
/*Usage:	
 	   ./Ascii.out 
	   *	<snapbase>
 	   *	<number of files per snap>
 	   *	<density of data sampling>
 	   *	<output filename> 
 	   *	<0-all or 1-gas or 2-DM>
*/

int main( int argc, char *argv[] )
{
    int i, file_snap, sampling, type;
    
    char snapbase[NMAX1], output[NMAX1];

    //Snapbase filename
    sprintf( snapbase, "%s", argv[1] );
    //Number of files per snap
    file_snap = atoi( argv[2] );
    //Density of data sampling
    sampling = atoi( argv[3] );
    //Output filename
    sprintf( output, "%s", argv[4] );
    //Type of data
    type = atoi( argv[5] );
        
    //Reading data from Gadget file
    read_snap( snapbase, file_snap, type );

    //Writing data
    if( type == 1 )	//Gas
	ascii_data_gas( Part, output, sampling );
    else		//All or DM
	ascii_data_all( Part, output, sampling, type );
    
    
    return 0;
}