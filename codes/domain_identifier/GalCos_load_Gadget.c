#include "GalCos_variables.h"

struct gadget_head
{
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
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 
							      256 Bytes */
};

struct gadget_head header1;

int GalCos_preload_Gadget(char *infile)
{
  
  int dummi, i, counter;
  FILE *fp_inp=NULL;
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&header1,sizeof(header1),1,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fclose(fp_inp);
  
  OMEGA_MATTER = header1.Omega0;
  OMEGALAMBDA  = header1.OmegaLambda;
  HUBBLEPARAM  = header1.HubbleParam;
  COSMIC_TIME  = header1.time;
  REDSHIFT     = header1.redshift;
  BoxSize      = header1.BoxSize;
  
  Npart_Total=0;
  for(i=0; i<6; i++)
    {
      if(task == 0) 
	printf(" * Header nall[%d] is: %d \n",i,header1.npartTotal[i]); fflush(stdout);
      
      Npart_Total = Npart_Total + header1.npartTotal[i];
    }
  
  if(task == 0) 
    printf("\n"); fflush(stdout);
  
  for(i=0; i<6; i++)
    { 
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	if(task == 0) 
	  printf(" * The mass of each particle is %d es %g\n",i,header1.mass[i]); fflush(stdout);
      
      if((header1.npart[i] != 0) && (header1.mass[i] == 0))
	if(task == 0) 
	  printf(" * There are individual mases for this particle set %d\n",i); fflush(stdout);
    }     
  
  counter=0;
  for(i=0; i<6; i++)
    {
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	counter++;
    }
  
  if(counter != 1)
    {
      printf("ERROR Multiple mass distribution or mistake reading file\n"); 
      exit(0);
    }
  
  PARTMASS = 0.0;
  
  for(i=0; i<6; i++)
    { 
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	PARTMASS = header1.mass[i];
    }  
  
  
  if(task == 0) 
    {
      printf("\n"); fflush(stdout);
      printf(" * Frame's Time... %g\n",header1.time); fflush(stdout);
      printf(" * Redshift... %g\n",header1.redshift); fflush(stdout);
      printf(" * Flagsfr... %d\n",header1.flag_sfr); fflush(stdout);
      printf(" * Flagfed... %d\n",header1.flag_feedback); fflush(stdout);
      printf(" * Flagcool... %d\n",header1.flag_cooling); fflush(stdout);
      printf(" * numfiles... %d\n",header1.num_files); fflush(stdout);
      printf(" * Boxsize... %g\n",header1.BoxSize); fflush(stdout);
      printf(" * Omega0... %g\n",header1.Omega0); fflush(stdout);
      printf(" * OmageLa... %g\n",header1.OmegaLambda); fflush(stdout);
      printf(" * Hubbleparam... %g\n",header1.HubbleParam); fflush(stdout);
      printf(" * Particle mass... %g\n",PARTMASS); fflush(stdout);
      printf(" * Total Number of particles... %d\n",Npart_Total); fflush(stdout);
    }
  
  return 0;
  
}


int GalCos_load_Gadget(int *unclust, int *sacum, char *infile)
{
  
  int dummi,i,j,counter,auxint,sorted_groups,maxsend,maxrecv,k;
  char buff[100];
  long int bloksize;
  FILE *aux_pf=NULL;
  FILE *fp_inp=NULL;
  struct part dp;
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  
  /* loading in task 0 */
  
  if(task == 0)
    {

      sprintf(buff,"%s%s",infile,".grp");      
      if( !(aux_pf = fopen(buff,"r")) ){
	    printf("load_gadget cannot open %s\n",buff);
	    exit(0);}
      fscanf(aux_pf,"%d",&auxint);
      
      if(auxint != Npart_Total)
	{
	  printf("Error reading files\n");
	  MPI_Finalize();
	  exit(0);
	}
      
      (*unclust) = 0;
      (*sacum) = 0;
      for(k=0; k<Npart_Total; k++) 
	{
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > (*sacum))
	    (*sacum) = sorted_groups;
	  if(sorted_groups == 0)
	    (*unclust) = (*unclust) + 1;
	}
      
      fclose(aux_pf);
      
      domain_info[task].istart = 0;
      domain_info[task].iend = Npart_Total;
      
      domain_info[task].Nparts_per_node = Npart_Total - (*unclust);
      printf("There are %d clustered particles\n",domain_info[task].Nparts_per_node); fflush(stdout);
      printf("There are %d FOF groups\n",(*sacum)); fflush(stdout);
      
      if(Particle != NULL) 
	free(Particle);
      
      if(domain_info[task].Nparts_per_node > 0)
	{
	  
	  Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
	  if(Particle == NULL){
	    printf("No memory available to load particles \n");
	    exit(0);
	  }
	  
	}
      
      //////////////////// Positions
      
      sprintf(buff,"%s%s",infile,".grp");
      aux_pf = fopen(buff,"r");
      fscanf(aux_pf,"%d",&auxint);
      
      bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
      fseek(fp_inp,bloksize,SEEK_SET);
      
      i=0;
      for(k=0; k<Npart_Total; k++) 
	{
	  
	  fread(&dp.pos[0],sizeof(float),3,fp_inp);
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > 0)
	    {
	      Particle[i].Cluster_ID = sorted_groups-1;
	      
	      for(j=0; j<3; j++)
		Particle[i].pos[j] = dp.pos[j];
	      //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
	      
	      Particle[i].mass = PARTMASS;
	      Particle[i].id = i;
	      Particle[i].Oid = k;
	      
	      i++;
	    }
	  
	}
      
      printf("i=%d %d\n",i,domain_info[task].Nparts_per_node); fflush(stdout);
      
      fclose(aux_pf);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      //////////////////// velocities
      
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      i=0;
      for(k=0; k<Npart_Total; k++) 
	{
	  fread(&dp.vel[0],sizeof(float),3,fp_inp);
	  
	  if((Particle[i].Oid == k) && (Particle[i].Cluster_ID >= 0))
	    {
	      for(j=0; j<3; j++) 
		Particle[i].vel[j] = dp.vel[j]; 
	      //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j]; 
	      
	      i++;
	    }
	  
	}
      
      //printf("i=%d %d\n",i,domain_info[task].Nparts_per_node); fflush(stdout);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
    }
    

  /* LODING DATA FOR TASKS != 0 */
  
  
  if(task != 0)
    {
      
      if(Particle != NULL) 
	free(Particle);
      
      if(domain_info[task].Nparts_per_node != 0)
	{
	  
	  Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
	  if(Particle == NULL){
	    printf("No memory available to load particles \n");
	    exit(0);
	  }
	  
	}
      
      //////////////////// Positions
      
      sprintf(buff,"%s%s",infile,".grp");
      aux_pf = fopen(buff,"r");
      fscanf(aux_pf,"%d",&auxint);
      
      // Positioning the pointers on the rigth position 
      
      for(k=0; k<domain_info[task].istart; k++)
	fscanf(aux_pf,"%d",&auxint);
      
      bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
      fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
      
      i=0;
      for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
      {
	  
	  fread(&dp.pos[0],sizeof(float),3,fp_inp);
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > 0)
	    Particle[i].Cluster_ID = sorted_groups-1;
	  else
	    Particle[i].Cluster_ID = EMPTY_FLAG;
	  
	  for(j=0; j<3; j++)
	    Particle[i].pos[j] = dp.pos[j];
	  //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
	  
	  Particle[i].mass = PARTMASS;
	  Particle[i].id = i;
	  Particle[i].Oid = k;
	  
	  i++;
      }
      
      fclose(aux_pf);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      //////////////////// velocities
      
      bloksize = sizeof(int)*5 + sizeof(struct gadget_head) + sizeof(float)*3*Npart_Total;
      fseek(fp_inp,bloksize + sizeof(float)*3*domain_info[task].istart,SEEK_SET);
      
      i=0;
      for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
	{
	  fread(&dp.vel[0],sizeof(float),3,fp_inp);
	  
	  for(j=0; j<3; j++)
	    Particle[i].vel[j] = dp.vel[j];
	  //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j];
	  
	  i++;
	}
      
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      (*unclust) = 0;
      (*sacum) = 0;
      
    }
  
  fclose(fp_inp);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;  
  
}


int reload_task0(char *infile)
{
  
  int dummi,i,j,counter,auxint,sorted_groups,maxsend,maxrecv,k;
  char buff[100];
  FILE *aux_pf=NULL;
  FILE *fp_inp=NULL;
  struct part dp;
  long int bloksize;
  
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  
  
  if(domain_info[task].Nparts_per_node != 0)
    {
      
      Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
      if(Particle == NULL){
	printf("No memory available to load particles \n");
	exit(0);
      }
      
    }
  
  //////////////////// Positions
  
  
  sprintf(buff,"%s%s",infile,".grp");
  aux_pf = fopen(buff,"r");
  fscanf(aux_pf,"%d",&auxint);
  
  // Positioning the pointers on the rigth position 
  
  bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
  fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
  
  i=0;
  for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
    {
      fread(&dp.pos[0],sizeof(float),3,fp_inp);
      fscanf(aux_pf,"%d",&sorted_groups);
      
      if(sorted_groups > 0)
	Particle[i].Cluster_ID = sorted_groups-1;
      else
	Particle[i].Cluster_ID = EMPTY_FLAG;
      
      for(j=0; j<3; j++) 
	Particle[i].pos[j] = dp.pos[j];
      //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
      
      Particle[i].mass = PARTMASS;
      Particle[i].id = i;
      Particle[i].Oid = k;
      
      i++;
    }
  
  fclose(aux_pf);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  //////////////////// velocities
  
  bloksize = sizeof(int)*5 + sizeof(struct gadget_head) + sizeof(float)*3*Npart_Total;
  fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  i=0;
  for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
    {
      
      fread(&dp.vel[0],sizeof(float),3,fp_inp);
      
      for(j=0; j<3; j++)
	Particle[i].vel[j] = dp.vel[j];
      //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j]; 
      
      i++;
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fclose(fp_inp);
  
  /*
    FILE *apf=fopen("ascii_snap.dat","w");
    for(i=0; i<Npart_Total; i++) 
    {
    fprintf(apf,"%g %g %g\n",Particle[i].pos[0],Particle[i].pos[1],Particle[i].pos[2]);
    }
    fclose(apf);
    exit(0);
  */
  
  return 0;  
  
}
