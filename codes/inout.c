#include <allvars.h>

/**************************************************************************************************
 NAME:	     read_head
 FUNCTION:   This function reads a header of a Gadget snapshot file
 INPUTS:     Filename
 RETURN:     0
 COMMENTS:   Based on the version by Juan Carlos Munoz-Cuartas
**************************************************************************************************/
int read_head( char *filename )
{
    int i, Npart_Total = 0, dummi;
    char label[4];
    FILE *fp_inp;

    //Opening input file
    if( (fp_inp=fopen(filename,"r"))==NULL ){
	printf( "It is not possible to open %s (not existing file)", filename );
	return 1;
	exit(0);}
    
    //Reading header
    fread(&dummi,sizeof(dummi),1,fp_inp);
    fread(&label,sizeof(char),4,fp_inp);
    fread(&dummi,sizeof(dummi),1,fp_inp);
    fread(&dummi,sizeof(dummi),1,fp_inp);
    fread(&dummi,sizeof(dummi),1,fp_inp);
    fread(&Gheader,sizeof(Gheader),1,fp_inp);
    fread(&dummi,sizeof(dummi),1,fp_inp);
  
    //Printing number of particles of each defined class
    for( i=0;i<6;i++ )
	printf(" * %10d Particles of class %d\n",Gheader.npart[i],i); 
    printf("\n");
        
    //Printing information about particles masses
    for(i=0; i<6; i++){
	//Individual masses
	if((Gheader.npart[i] != 0) && (Gheader.mass[i] != 0.0))
	    printf(" * The mass of each particle is %d es %g\n",i,Gheader.mass[i]);
	//Individual masses for each set of particles
	if((Gheader.npart[i] != 0) && (Gheader.mass[i] == 0.0))
	    printf(" * There are individual mases for this particle set %d\n",i);}
    printf("\n");
    
    //Information of this snapshot
    printf(" * Frame's Time... %g\n",Gheader.time); 
    printf(" * Redshift... %g\n",Gheader.redshift);
    printf(" * Flagsfr... %d\n",Gheader.flag_sfr);
    printf(" * Flagfed... %d\n",Gheader.flag_feedback);
    printf(" * Flagcool... %d\n",Gheader.flag_cooling);
    printf(" * numfiles... %d\n",Gheader.num_files);
    printf(" * Boxsize... %g\n",Gheader.BoxSize);
    printf(" * Omega0... %g\n",Gheader.Omega0);
    printf(" * OmageLa... %g\n",Gheader.OmegaLambda);
    printf(" * Hubbleparam... %g\n",Gheader.HubbleParam);
    printf(" * Npart_Total=%d\n\n",Npart_Total);
    
    //Counting number of particles (nall for number>2^32 particles)
    for(i=0; i<6; i++){
	Npart_Total += Gheader.npartTotal[i];}

    
    //Allocating memory for the particle array, (including all single files of the snapshot)
    if( !(Part = (struct part *) malloc((size_t) Npart_Total*sizeof(struct part))) ){
	printf("No memory available for load the particles (%s)\n",filename);
	return 1;
	exit(0);}
   
    //Used memory
    printf(" * Using %f Mb in memory\n\n",Npart_Total*sizeof(struct part)/(1020*1024.0));

    return 0;
}


/**************************************************************************************************
 NAME:	     read_snap
 FUNCTION:   This function reads a complete Gadget snap (only gas particles) file in format 2.
 INPUTS:     Snapbase of the current snapshot, number of single component files,
	     type of particle to be read:   0-ALL      1-GAS
 RETURN:     number of particles in the snap (0 if some error)
 COMMENTS:   Based on the version of Juan Carlos Munoz-Cuartas
**************************************************************************************************/
int read_snap( char snapbase[NMAX1], 
	       int filenumber,
	       int type )
{
    int i, s, dummi, global_acum = 0;
    
    int id;
    int Npart_snap_fin;
    float pos[3], vel[3], mass, energy, rho;
        
    char label[4], filename[NMAX1];
    long int blksize_all,
	     blksize_pos, blksize_vel, blksize_ids, 
	     blksize_mss, blksize_ene, blksize_rho;
        
    FILE *fp_pos;
    FILE *fp_vel;
    FILE *fp_ids;
    FILE *fp_mss;
    FILE *fp_ene;
    FILE *fp_rho;
    
    FILE *fp_head;
    FILE *pf_out;
        
    printf("==========================================================================\n");
    printf(" Reading from gadget file %s \n",filename);
    printf("==========================================================================\n");
    
    //Sweeping all component files of this snapshot
    for( s=0;s<filenumber;s++ ){

	//Filename of current file
	sprintf( filename,"%s.%d",snapbase,s );
	if( filenumber == 1 )
	    sprintf( filename,"%s",snapbase );
	  
	//Reading first file of the sequence for this snapshot
	if( s==0 )
	    read_head( filename );

	//Openning current file
	if( !(fp_head=fopen(filename,"r")) ){
	    printf("read_gadget cannot open %s\n",filename);
	    exit(0);}

	//Different blocks of the Gadget snapshot file
	fp_pos=fopen( filename,"r" );
	fp_vel=fopen( filename,"r" );
	fp_ids=fopen( filename,"r" );
	fp_mss=fopen( filename,"r" );
	fp_ene=fopen( filename,"r" );
	fp_rho=fopen( filename,"r" );
	
	//Reading header of the s-th file of the current snapshot
	fread(&dummi,sizeof(dummi),1,fp_head);
	fread(&label,sizeof(char),4,fp_head);
	fread(&dummi,sizeof(dummi),1,fp_head);
	fread(&dummi,sizeof(dummi),1,fp_head);
	fread(&dummi,sizeof(dummi),1,fp_head);
	fread(&Gheader,sizeof(Gheader),1,fp_head);
	fread(&dummi,sizeof(dummi),1,fp_head);
	fclose(fp_head);
	
	//Counting number of particles (with and without variable masses)
	Npart_snap = 0;
	Npart_snap_mass = 0;
	for(i=0; i<6; i++){
	    Npart_snap += Gheader.npart[i];
	    if( Gheader.mass[i] == 0.0 )
		Npart_snap_mass += Gheader.npart[i];}
				
	//Seeking different blocks in the snapshot file in order to read them
	blksize_all = 6*sizeof(int)  + 1*4*sizeof(char) + sizeof(Gheader);
	//Position
	blksize_pos = blksize_all;
	fseek(fp_pos,blksize_pos,SEEK_SET);
	fread(&label,sizeof(char),4,fp_pos);
	printf( "Block %s\n", label );
	blksize_pos += 3*sizeof(int)  + 1*4*sizeof(char);
	blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + 3*Npart_snap*sizeof(float);
	
	//Velocities
	blksize_vel = blksize_all;
	fseek(fp_vel,blksize_vel,SEEK_SET);
	fread(&label,sizeof(char),4,fp_vel);
	printf( "Block %s\n", label );
	blksize_vel += 3*sizeof(int)  + 1*4*sizeof(char);
	blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + 3*Npart_snap*sizeof(float);
	
	//Identifiers
	blksize_ids = blksize_all;
	fseek(fp_ids,blksize_ids,SEEK_SET);
	fread(&label,sizeof(char),4,fp_ids);
	printf( "Block %s\n", label );
	blksize_ids += 3*sizeof(int)  + 1*4*sizeof(char);
	blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Npart_snap*sizeof(int);
	
	//Seeking in files
	fseek(fp_pos,blksize_pos,SEEK_SET);
	fseek(fp_vel,blksize_vel,SEEK_SET);
	fseek(fp_ids,blksize_ids,SEEK_SET);

	//Only if GAS particles
	if( type == 1 ){
	    //Mass
	    blksize_mss = blksize_all;
	    fseek(fp_mss,blksize_mss,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_mss);
	    printf( "Block %s\n", label );
	    blksize_mss += 3*sizeof(int)  + 1*4*sizeof(char);
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Npart_snap_mass*sizeof(float);

	    //Energy
	    blksize_ene = blksize_all;
	    fseek(fp_ene,blksize_ene,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_ene);
	    printf( "Block %s\n", label );
	    blksize_ene += 3*sizeof(int)  + 1*4*sizeof(char);
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);
	    
	    //Density
	    blksize_rho = blksize_all;
	    fseek(fp_rho,blksize_rho,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_rho += 3*sizeof(int)  + 1*4*sizeof(char);
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);

	    //Other blocks=============================================================================
	    fseek(fp_rho,blksize_all,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);
	    
	    fseek(fp_rho,blksize_all,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);

	    fseek(fp_rho,blksize_all,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);
	    
	    fseek(fp_rho,blksize_all,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);
	    
	    fseek(fp_rho,blksize_all,SEEK_SET);
	    fread(&label,sizeof(char),4,fp_rho);
	    printf( "Block %s\n", label );
	    blksize_all += 5*sizeof(int)  + 1*4*sizeof(char) + Gheader.npart[0]*sizeof(float);
	    //=========================================================================================
	    
	    //Seeking in files
	    fseek(fp_mss,blksize_mss,SEEK_SET);
	    fseek(fp_ene,blksize_ene,SEEK_SET);
	    fseek(fp_rho,blksize_rho,SEEK_SET);}
	
	Npart_snap_fin = Npart_snap;
	//Setting only gas particles
	if( type == 1 ){
	    Npart_snap_fin = Gheader.npart[0];}
	
	//Reading positions, velocities and id of each particle	
	for(i=0; i<Npart_snap_fin ;i++){
	    //Reading from binary files
	    fread(&pos[0],sizeof(float),3,fp_pos);
	    fread(&vel[0],sizeof(float),3,fp_vel);
	    fread(&id,sizeof(int),1,fp_ids);

	    //Positions
	    Part[i].pos[0] = pos[0];
	    Part[i].pos[1] = pos[1];
	    Part[i].pos[2] = pos[2];
 	    //Velocities
	    Part[i].vel[0] = vel[0];
	    Part[i].vel[1] = vel[1];
	    Part[i].vel[2] = vel[2];
 	    //IDs
	    Part[i].id = i;
	    
	    //GAS PROPERTIES
	    if( type == 1 ){
		//Reading from binary files
		fread(&mass,sizeof(float),1,fp_mss);
		fread(&energy,sizeof(float),1,fp_ene);
		fread(&rho,sizeof(float),1,fp_rho);
		
		//Mass
		Part[i].mass = mass;
		//Energy
		Part[i].energy = energy;
		//Density
		Part[i].rho = rho;
		//Pressure
		Part[i].pressure = pressure( Part[i].rho, Part[i].energy );
		//Temperature
		Part[i].temperature = temperature( Part[i].energy );}
	    
	    //Redshift
	    Part[i].z = Gheader.redshift;
	    //Time
	    Part[i].t = Gheader.time;
	    
	    global_acum++;}

	fclose(fp_pos);
	fclose(fp_vel);
	fclose(fp_ids);
	fclose(fp_mss);
	fclose(fp_ene);
	fclose(fp_rho);

	printf("==========================================================================\n");
	printf(" End of reading from gadget file %s, %d \n",filename, global_acum);
	printf("==========================================================================\n");}

    return 0;
}


/**************************************************************************************************
 NAME:	     ascii_data_gas
 FUNCTION:   This function write positions, velocities and ids and gas variables of each gas 
	     particle onto an ascci file
 INPUTS:     Struct with particles, number of particles, output filename, number of data sampling
 RETURN:     0
**************************************************************************************************/
int ascii_data_gas( struct part *parts,
		    char output[NMAX1],
		    int sampling )
{
    int i;
    FILE *out;
    
    //Creating output file
    out = fopen( output, "w" );

    //Writing data
    fprintf( out, 
    "#Id\tX\t\tY\t\tZ\t\tVX\t\tVY\t\tVZ\t\tMass\t\tU\t\tRHO\t\tP\t\tT\t\tz\t\ttime\n" );
    
    for( i=0;i<Gheader.npartTotal[0];i+=sampling ){
	fprintf( out, 
		 "%d\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n", 
		 parts[i].id,
		 parts[i].pos[X], parts[i].pos[Y], parts[i].pos[Z],
		 parts[i].vel[X], parts[i].vel[Y], parts[i].vel[Z],
		 parts[i].mass, parts[i].energy, parts[i].rho, parts[i].pressure, parts[i].temperature,
		 parts[i].z, parts[i].t);}
    
    fclose( out );
    
    return 0;
}


/**************************************************************************************************
 NAME:	     ascii_data_all
 FUNCTION:   This function write positions, velocities and ids of each particle onto an ascci file
 INPUTS:     Struct with particles, number of particles, output filename, number of data sampling,
	     type of particle to be read:   0-ALL  2-DM
 RETURN:     0
**************************************************************************************************/
int ascii_data_all( struct part *parts,
		    char output[NMAX1],
		    int sampling,
		    int type )
{
    int i;
    int N0, Nf;
    FILE *out;
    
    //Creating output file
    out = fopen( output, "w" );
    
    //All particles
    if( type == 0 ){
	N0 = 0;
	Nf = Npart_snap;}
    //DM particles
    if( type == 2 ){
	N0 = Gheader.npartTotal[0];
	Nf = Gheader.npartTotal[1]+Gheader.npartTotal[0];}

    //Writing data
    fprintf( out, "#Id\tX\t\tY\t\tZ\t\tVX\t\tVY\t\tVZ\t\tz\t\ttime\n" );
    
    for( i=N0;i<Nf;i+=sampling ){
	fprintf( out, "%d\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n", 
		 parts[i].id,
		 parts[i].pos[X], parts[i].pos[Y], parts[i].pos[Z],
		 parts[i].vel[X], parts[i].vel[Y], parts[i].vel[Z],
		 parts[i].z, parts[i].t);}
    
    fclose( out );
    
    return 0;
}


/**************************************************************************************************
 NAME:	     ascii_data_pos
 FUNCTION:   This function write positions of each particle onto an ascci file
 INPUTS:     Struct with particles, number of particles, output filename, number of data sampling,
	     type of particle to be read:   0-ALL      1-GAS      2-DM
 RETURN:     0
**************************************************************************************************/
int ascii_data_pos( struct part *parts,
		    char output[NMAX1],
		    int sampling,
		    int type )
{
    int i, it;
    int N0, Nf;
    FILE *out;
    
    //Creating output file
    out = fopen( output, "w" );
    
    //All particles
    if( type == 0 ){
	N0 = 0;
	Nf = Npart_snap;}
    //GAS particles
    if( type == 1 ){
	N0 = 0;
	Nf = Gheader.npartTotal[0];}
    //DM particles
    if( type == 2 ){
	N0 = Gheader.npartTotal[0];
	Nf = Gheader.npartTotal[1]+Gheader.npartTotal[0];}

    //Writing data
    fprintf( out, "%d\n", (Nf-N0)/sampling ); //Points in total
    fprintf( out, "%d\n", (Nf-N0)/sampling ); //Points in DM
    fprintf( out, "0\n" ); //Gas
    fprintf( out, "0.01\n" ); //Time
    fprintf( out, "0\n" ); //Nactive
    
    for( i=N0;i<Nf;i+=sampling ){
	fprintf( out, "%1.5e\t%1.5e\t%1.5e\n", 
		 parts[i].pos[X], parts[i].pos[Y], parts[i].pos[Z]);}
    
    fclose( out );
    
    return 0;
}


/**************************************************************************************************
 NAME:	     ascii_data_slide
 FUNCTION:   This function write positions of each particle within a given slide on an ascci file
 INPUTS:     Struct with particles, number of particles, output filename, number of data sampling,
	     type of particle to be read:   0-ALL      1-GAS      2-DM, axis, slide value, thick.
 RETURN:     0
**************************************************************************************************/
int ascii_data_slide( struct part *parts,
		      char output[NMAX1],
		      int sampling,
  		      int type,
  		      int axis,
		      float slide,
		      float dx )
{
    int i, it;
    int N0, Nf;
    FILE *out;
    
    //Creating output file
    out = fopen( output, "w" );
    
    //All particles
    if( type == 0 ){
	N0 = 0;
	Nf = Npart_snap;}
    //GAS particles
    if( type == 1 ){
	N0 = 0;
	Nf = Gheader.npartTotal[0];}
    //DM particles
    if( type == 2 ){
	N0 = Gheader.npartTotal[0];
	Nf = Gheader.npartTotal[1]+Gheader.npartTotal[0];}
    
    for( i=N0;i<Nf;i+=sampling )
	if( Part[i].pos[axis] >= slide && Part[i].pos[axis] < slide + dx ){
	    fprintf( out, "%1.5e\t%1.5e\t%1.5e\n", 
		    parts[i].pos[X], parts[i].pos[Y], parts[i].pos[Z]);}
    
    fclose( out );
    
    return 0;
}