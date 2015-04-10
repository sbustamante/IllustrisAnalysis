#include <allvars.h>
/*Usage:	
 	   ./Halos.out 
	   *	<snapbase>
 	   *	<number of files per snap>
 	   *	<output filename> 
 	   *	<0-all or 1-gas or 2-DM>
 	   *	<linking lenght>
 	   *	<minim number of particles>
 	   *	<shift gas particles>
*/

int main( int argc, char *argv[] )
{
    int i, file_snap, type;
    long int n_slide;
    float linking, linkingtot;
    int min_parts, shift;

    char snapbase[NMAX1], output[NMAX1], cmd[NMAX1];

    //Snapbase filename
    sprintf( snapbase, "%s", argv[1] );
    //Number of files per snap
    file_snap = atoi( argv[2] );
    //Output filename
    sprintf( output, "%s", argv[3] );
    //Type of data
    type = atoi( argv[4] );
    //Linking lenght for FOF scheme
    linking = atof( argv[5] );
    //Minim number of particles for constructing a halo
    min_parts = atoi( argv[6] );
    //Shift gas particles
    shift = atoi( argv[7] );
    
    //Reading data from Gadget file
    read_snap( snapbase, file_snap, type );

    //Writing temporal data with positions
    ascii_data_pos( Part, "pos.tmp", 1, type );
    
    //Runing FOF scheme
    printf("\n==========================================================================\n");
    printf(" Identifying Halos\n");
    printf("==========================================================================\n");
    if( type == 0 )
	linkingtot = linking*pow(pow(Gheader.BoxSize,3)/Npart_snap,1/3.0);
    if( type == 1 )
	linkingtot = linking*pow(pow(Gheader.BoxSize,3)/Gheader.npartTotal[0],1/3.0);
    if( type == 2 )
	linkingtot = linking*pow(pow(Gheader.BoxSize,3)/Gheader.npartTotal[1],1/3.0);
    sprintf( cmd, "./fofscr/fof -e %f -m %d < ./pos.tmp", linkingtot, min_parts );
    printf( "%s\n", cmd );
    system( cmd );
    system( "rm ./pos.tmp" );
    if( shift == 1 ){
	sprintf( cmd, "bash fof_shifter.sh %s %d %d %d", 
		 "fof.grp", Gheader.npartTotal[0], Gheader.npartTotal[1], Gheader.npartTotal[4] );
	printf( "%s\n", cmd );
	system( cmd );}
      
    

    return 0;
}