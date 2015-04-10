/* 
   This program takes a snapshot of a simulations, the list of FOF
   halos and a given mass threshold and computes the domain of each
   halo based on the virial mass and virial radius of the halo. With
   that, the code identifies the domain PARTICLES of every halo in
   that snapshot.
   
   Every particle in the box is asigned to a domain.

   this version of the code also runs the calculation of the domain of
   the halo with periodic boundary conditions.
   
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

#include "GalCos_variables.h"
#include "GalCos_variables.c"

int task,Number_Processes;

/* Data structure for the domain decomposition */

struct domain
{
  int Nparts_per_node;
  int istart;
  int iend;
} *domain_info;

void usage_main(void)
{
  
  if(task == 0)
    {
      printf("Domain_identifier = Identifies the domain of FOF halos.\n"); fflush(stdout);
      printf("Juan Carlos Munoz C.\n"); fflush(stdout);
      printf("Usage:./Domain_identifier <Gadget_snapshoth_file> <param_file>\n"); fflush(stdout);
    }
  
  MPI_Finalize();
  exit(0);
  
}

#include "GalCos_routines.c"
#include "GalCos_load_Gadget.c"

#include "GalCos_units.c"
#include "GalCos_cola.c"
#include "GalCos_clusterBound.c"

#include "GalCos_Halo_prop.c"
#include "GalCos_domain.c"
#include "kill_low_mass_halos.c"

struct data
{
  float mass;
  float pos[3];
  float Rvir;
  int Nmembers;
  int NDomain_particles; // Total number of particles in the domain = domain + envelope
  int IDcluster;
  int NDomain_particles_absolute; // Only the number of particles in the domain are here
};

int write_rescue(char *infile)
{
  
  int i,*Halo_particles=NULL,*Domain_particles=NULL,j;
  FILE *pf=NULL;
  char buf[200];
  struct data aux;
  
  sprintf(buf,"%s%s",infile,".rescue");
  pf=fopen(buf,"w");
  
  fwrite(&NCLUSTERS,sizeof(int),1,pf);
  
  for(i=0; i<NCLUSTERS; i++)
    {
      aux.mass = Halos[i].mass;
      aux.pos[0] = Halos[i].pos[0];
      aux.pos[1] = Halos[i].pos[1];
      aux.pos[2] = Halos[i].pos[2];
      aux.Nmembers = Halos[i].Nmembers;
      aux.NDomain_particles = Halos[i].NDomain_particles;
      aux.NDomain_particles_absolute = Halos[i].NDomain_particles_absolute;
      aux.IDcluster = Halos[i].IDcluster;
      aux.Rvir = Halos[i].Rvir;
      
      fwrite(&aux,sizeof(struct data),1,pf);
      
      Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
      
      for(j=0; j<Halos[i].Nmembers; j++)
	Halo_particles[j] = Halos[i].Halo_particles[j];
      
      fwrite(&Halo_particles[0],sizeof(int),Halos[i].Nmembers,pf);
      
      free(Halo_particles);
      
      Domain_particles = (int *) malloc((size_t) Halos[i].NDomain_particles*sizeof(int));
      
      for(j=0; j<Halos[i].NDomain_particles; j++)
	Domain_particles[j] = Halos[i].Domain_particles[j];
      
      fwrite(&Domain_particles[0],sizeof(int),Halos[i].NDomain_particles,pf);
      
      free(Domain_particles);
    }
  
  fclose(pf);
  
  return 0;
}


int main(int argc, char *argv[])
{
  
  int i,counter,iadvance,NGROUPS=0,UNCLUSTERED,Number_Processes,pos;
  int ihalo,Oistart,Oiend,ONparts_per_node,Nparts_in_node,NclusteredParts;
  int COLLECT_TAG=1,k,j,*Info_domains=NULL,sendbuff;
  int NDomain_parts_per_proc,Parts_in_Domain_per_node,l,m,Old_NDomain_particles;
  int Old_NDomain_particles_absolute,NDomain_parts_per_proc_absolute;
  float Mth;

  char *infile=NULL, *param_file=NULL;
  char buf[200];
  
  MPI_Status status;
  FILE *aux_pf=NULL;
  FILE *pf=NULL;
  
  struct halo *auxHalos=NULL;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  MPI_Comm_size(MPI_COMM_WORLD,&Number_Processes);
  
  if(task == 0)
    printf(" * Running on %d processors\n",Number_Processes); fflush(stdout);
  
  if(argc != 4)
    usage_main();
  
  infile = argv[1];
  param_file = argv[2];
  Mth = atof(argv[3]);
  
  Mth = pow(10.0,Mth)/(1.0e10); // passing the mass thrsehold to
				// internal units
  
  if(task == 0) 
    printf("Mth = %g\n",Mth); fflush(stdout);
  
  GalCos_units(param_file);
  GalCos_preload_Gadget(infile);
  
  
  /////////////////////////////////////////////////////////////////////////////
  /*                          DOMAIN DECOMPOSITION                           */
  /////////////////////////////////////////////////////////////////////////////
  
  
  domain_info = (struct domain *) malloc((size_t) Number_Processes*sizeof(struct domain));
  
  if(task == 0)
    {
      
      Nparts_in_node = (int) (Npart_Total/Number_Processes);
      
      printf("\n There are %d particles per node\n",Nparts_in_node); fflush(stdout);
      printf(" There are left particles %d\n",(Npart_Total%Number_Processes)); fflush(stdout);
      
      for(i=0; i<Number_Processes; i++)
	{
	  domain_info[i].istart = i*Nparts_in_node;
	  domain_info[i].iend = (i+1)*Nparts_in_node;
	  domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
	}
      
      i = Number_Processes-1;
      domain_info[i].iend = Npart_Total;
      domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
      
      i = 0;
      Oistart = domain_info[i].istart;
      Oiend = domain_info[i].iend;
      ONparts_per_node = domain_info[i].Nparts_per_node;
      
      for(i=0; i<Number_Processes; i++)
	printf(" task %d start at %d and ends at %d, nparts=%d\n",task,domain_info[i].istart,
	       domain_info[i].iend,domain_info[i].Nparts_per_node); fflush(stdout);
      
    }
  
  MPI_Bcast(&domain_info[0],Number_Processes*sizeof(struct domain),MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  //                 STARTING WITH FOF HALO IDENTIFICATION                           // 
  /////////////////////////////////////////////////////////////////////////////////////
  
  if(task == 0)
    {
      
      printf(" >> LOADING PARTICLE DATA FROM GADGET FILE...\n"); fflush(stdout);
      
      GalCos_load_Gadget(&UNCLUSTERED,&NCLUSTERS,infile);

      
      printf("\n %g Memory in particles = %g %d\n",sizeof(struct part)/(1024*1024.0),
	     Npart_Total*sizeof(struct part)/(1024*1024.0),Npart_Total); fflush(stdout);
      printf(" %g Memory in halos = %g\n",sizeof(struct halo)/(1024*1024.0),
	     60000*sizeof(struct halo)/(1024*1024.0)); fflush(stdout);
      printf(" THERE ARE %d UNCLUSTERED PARTICLES (%f percent)\n",UNCLUSTERED,
	     100.0*(1.0*UNCLUSTERED)/(Npart_Total*1.0)); fflush(stdout);
      printf(" There are %d fof halos\n",NCLUSTERS); fflush(stdout);
      
      
      /////////////////////////////////////////////////////////////////////////////////////
      //                 FILLING HALO STRUCTURES WITH PARTICLE ID's                      //
      /////////////////////////////////////////////////////////////////////////////////////
      
      /* INITILIZING GROUPS IN TASK 0 ONLY */
      
      Halos = (struct halo *) malloc((size_t) NCLUSTERS*sizeof(struct halo));
      if(Halos == NULL)
	{
	  printf("There are no memory left to allocate Halos\n");
	  MPI_Finalize();
	  exit(0);
	}
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  Halos[i].Nmembers = 0;
	  Halos[i].NDomain_particles = 0;
	  Halos[i].NDomain_particles_absolute = 0;
	  Halos[i].Halo_particles = NULL;
	  Halos[i].Domain_particles = NULL;
	  Halos[i].Domain_particles_index = NULL;
	}
      
      
      ////////////////////////////////////////////
      
      iadvance = 10000000;
      
      for(i=0; i<domain_info[task].Nparts_per_node; i++)
	{
	  
	  if((i%iadvance) == 0)
	    printf(" *Evaluated %d of %d\n",i,Npart_Total); fflush(stdout);
	  
	  if(Particle[i].Cluster_ID != EMPTY_FLAG)
	    {
	      ihalo = Particle[i].Cluster_ID;
	      Halos[ihalo].Nmembers++;
	      Halos[ihalo].Halo_particles = realloc(Halos[ihalo].Halo_particles, (size_t) Halos[ihalo].Nmembers*sizeof(int));
	      Halos[ihalo].Halo_particles[Halos[ihalo].Nmembers-1]= Particle[i].id;
	    }
	  
	}
      
      
      aux_pf = fopen("Mass_function_fof.dat","w");
      for(i=0; i<NCLUSTERS; i++)
	fprintf(aux_pf,"%g\n",Halos[i].Nmembers*PARTMASS);
      fclose(aux_pf);
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  Halos[i].ID_CenterHalo = Halos[i].Halo_particles[0];
	  Halos[i].IDcluster = i;
	  Halos[i].Domain_particles = NULL;
	  Halos[i].NDomain_particles = 0;
	}
      
      printf("There are %d FOF halos\n",NCLUSTERS); fflush(stdout);
      
      
      /////////////////////////////////////////////////////////////////////////////////////
      //                  CENTERING AND PERIODIC BOUNDARY CONDITIONS                     //
      /////////////////////////////////////////////////////////////////////////////////////
      
      iadvance = 10000;
      auxHalos = (struct halo *) malloc((size_t) NCLUSTERS*sizeof(struct halo));
      
      aux_pf = fopen("Halo_Catalog.dat", "w");

      k = 0;
      for(i=0; i<NCLUSTERS; i++)
	{
	  
	  if((i%iadvance) == 0)
	    printf(" *Computed Properties for %d\n",i); 
	  fflush(stdout);
	  
	  /* 
	     Doing corrections of periodic boundary, centering the
	     halo and shifthing particle coords and computing Mvir and
	     Rvir
	  */
	  periodic_boundary_corrections(i);
	  Halos[i].ID_CenterHalo = GalCos_clusterBound(Halos[i].Halo_particles,Halos[i].Nmembers,i);
	  Halo_Mvir_Rvir(i);
	  
	  fprintf(aux_pf,"%16.8f %16.8f %16.8f %16.8f %16.8f %12d\n", Halos[i].pos[0], Halos[i].pos[1], Halos[i].pos[2],
		  Halos[i].Mvir, Halos[i].Rvir, Halos[i].Nvir);

	  /* Killing halos with Mvir lower mass than Mth (Mvir < Mth) */
	  //kill_low_mass_halos(Mth,i);
	  
	  if(Halos[i].Nmembers > 0)
	    {
	      auxHalos[k].Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
	      
	      for(j=0; j<Halos[i].Nmembers; j++)
		auxHalos[k].Halo_particles[j] = Particle[Halos[i].Halo_particles[j]].Oid;
	      
	      auxHalos[k].pos[0] = Halos[i].pos[0];
	      auxHalos[k].pos[1] = Halos[i].pos[1];
	      auxHalos[k].pos[2] = Halos[i].pos[2];

	      auxHalos[k].vel[0] = Halos[i].vel[0];
	      auxHalos[k].vel[1] = Halos[i].vel[1];
	      auxHalos[k].vel[2] = Halos[i].vel[2];
	      
	      auxHalos[k].mass     = Halos[i].mass;
	      auxHalos[k].Rvir     = Halos[i].Rvir;
	      auxHalos[k].Mvir     = Halos[i].Mvir;
	      auxHalos[k].Nvir     = Halos[i].Nvir;
	      auxHalos[k].Nmembers = Halos[i].Nmembers;
	      
	      auxHalos[k].ID_CenterHalo = Halos[i].ID_CenterHalo;
	      auxHalos[k].NDomain_particles = Halos[i].NDomain_particles;
	      auxHalos[k].IDcluster = Halos[i].IDcluster;
	      auxHalos[k].NDomain_particles_absolute = Halos[i].NDomain_particles_absolute;
	      
	      k++;
	    }
	  
	  free(Halos[i].Halo_particles);
	  Halos[i].Halo_particles = NULL;
	  
	}
      
      free(Particle);
      fclose(aux_pf);
      
      /* updating info of survival halos */
      printf("\n The number of halos changed from %d to %d\n",NCLUSTERS,k); fflush(stdout);
      NCLUSTERS = k;
      Halos = (struct halo *) realloc(Halos, (size_t) NCLUSTERS*sizeof(struct halo));
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  Halos[i].mass   = auxHalos[i].mass;
	  Halos[i].pos[0] = auxHalos[i].pos[0];
	  Halos[i].pos[1] = auxHalos[i].pos[1];
	  Halos[i].pos[2] = auxHalos[i].pos[2];
	  
	  Halos[i].vel[0] = auxHalos[i].vel[0];
	  Halos[i].vel[1] = auxHalos[i].vel[1];
	  Halos[i].vel[2] = auxHalos[i].vel[2];
	  
	  Halos[i].Rvir = auxHalos[i].Rvir;
	  Halos[i].Mvir = auxHalos[i].Mvir;
	  Halos[i].Nvir = auxHalos[i].Nvir;
	  Halos[i].Nmembers = auxHalos[i].Nmembers;
	  
	  Halos[i].ID_CenterHalo = auxHalos[i].ID_CenterHalo;
	  Halos[i].NDomain_particles = auxHalos[i].NDomain_particles;
	  
	  Halos[i].IDcluster = i;
	  Halos[i].NDomain_particles_absolute = auxHalos[i].NDomain_particles_absolute;
	  
	  Halos[i].Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
	  
	  for(j=0; j<Halos[i].Nmembers; j++)
	    Halos[i].Halo_particles[j] = auxHalos[i].Halo_particles[j];
	  
	  free(auxHalos[i].Halo_particles);
	}
      
      
      /*
	aux_pf = fopen("Halo_Catalog.dat", "w");
	for(i=0; i<NCLUSTERS; i++)
	fprintf(aux_pf,"%16.8f %16.8f %16.8f %16.8f %16.8f %12d\n", Halos[i].pos[0], Halos[i].pos[1], Halos[i].pos[2],
	Halos[i].Mvir, Halos[i].Rvir, Halos[i].Nvir);
	fclose(aux_pf);
      */
      
      /////////////////////////////////////
      
      domain_info[task].istart = Oistart;
      domain_info[task].iend = Oiend;
      domain_info[task].Nparts_per_node = ONparts_per_node;
      
      reload_task0(infile);
      
      free(auxHalos);
      
      
      printf(" Done halo properties\n"); fflush(stdout);
      
    }//if task == 0
  
  //MPI_Finalize();
  //exit(0);
  
  MPI_Barrier(MPI_COMM_WORLD);
  //////////////////////////////////////////////////////////////
  /*                BROADCASTING HALO PROPERTIES              */
  //////////////////////////////////////////////////////////////
  
  if(task == 1)
    printf(" Loading particle data on all other tasks\n"); fflush(stdout);
  
  if(task > 0)
    GalCos_load_Gadget(&UNCLUSTERED,&NGROUPS,infile);
  
  MPI_Bcast(&NCLUSTERS,1,MPI_INT,0,MPI_COMM_WORLD);
  
  /*  INITILIZING HALOS IN OTHER TASKS > 0 */
  
  if(task > 0)
    {
      
      Halos = (struct halo *) malloc((size_t) NCLUSTERS*sizeof(struct halo));
      if(Halos == NULL)
	{
	  printf("There are no memory left to allocate Halos\n");
	  MPI_Finalize();
	  exit(0);
	}
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  Halos[i].Nmembers = 0;
	  Halos[i].NDomain_particles = 0;
	  Halos[i].NDomain_particles_absolute = 0;
	  Halos[i].Halo_particles = NULL;
	  Halos[i].Domain_particles = NULL;
	  Halos[i].Domain_particles_index = NULL;
	}
      
    }
  
  MPI_Bcast(&Halos[0],NCLUSTERS*sizeof(struct halo),MPI_BYTE,0,MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0)
    printf("Broadcasted halos %d was done...\n",NCLUSTERS); fflush(stdout);
  
  for(i=0; i<NCLUSTERS; i++)
    {
      
      if(task > 0)
	{
	  Halos[i].Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
	  
	  if(Halos[i].Halo_particles == NULL)
	    printf("Not able to allocate memory for particles in halo %d\n",i);
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&(Halos[i].Halo_particles[0]),Halos[i].Nmembers,MPI_INT,0,MPI_COMM_WORLD);
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0)
    printf("Broadcast of halos was done...\n"); fflush(stdout);
  
  
  /////////////////////////////////////////////////////////////////////////////////
  /*      AFTER REMOVING PARTICLES OF HALOS IN THE FOF HALO (FOR VIR
	  HALOS) I HAVE TO UPDATE THE HALO ID OF PARTICLES REMOVED
	  FROM THOSE HALOS  */
  /////////////////////////////////////////////////////////////////////////////////
  
  NclusteredParts = 0;
  for(i=0; i<NCLUSTERS; i++)
    {
      Halos[i].Domain_particles  = NULL;
      Halos[i].NDomain_particles = 0;
      NclusteredParts += Halos[i].Nmembers;
    }
  
  if(task == 0)
    printf("There are %d particles in halos (Nvir)\n",NclusteredParts); fflush(stdout);

  /* array with all particles in a halo */
  int *PartsInHalos = (int *) malloc((size_t) NclusteredParts*sizeof(int));
  
  i=0;
  for(j=0; j<NCLUSTERS; j++)
    {
      for(k=0; k<Halos[j].Nmembers; k++)
	{
	  PartsInHalos[i] = Halos[j].Halo_particles[k];
	  i++;
	}
      
    }
  
  gsl_sort_int(PartsInHalos,1,NclusteredParts);
  
  /* if the particle in the domain is not in the array, is not in a
     halo! */
  for(i=0; i<domain_info[task].Nparts_per_node; i++)
    {
      pos = mybsearch(PartsInHalos,NclusteredParts,Particle[i].Oid);
      if(pos == -1)
	Particle[i].Cluster_ID = EMPTY_FLAG;
    }
  
  
  ////////////////////////////////////////////////////////////////////////
  /*                  LOOKING FOR THE DOMAIN OF THE HALOS               */
  ////////////////////////////////////////////////////////////////////////
  
  counter = 0;
  for(i=0; i<domain_info[task].Nparts_per_node; i++)
    {
      
      if((counter%1000000) == 0)
	printf(" (%d) *Computed domains for %d\n",task,counter); fflush(stdout);
      
      if(Particle[i].Cluster_ID == EMPTY_FLAG)
	GalCos_domain(i);
      
      counter++;
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  

  /////////////////////////////////////////////////////////////////////////////
  /*               WRITING DOMAIN INFORMATION TO RESCUE FILE                 */
  /////////////////////////////////////////////////////////////////////////////
  
  for(l=0; l<Number_Processes; l++)
    {
      
      if(task == l)
	{
	  sprintf(buf,"%s%s",infile,".parts.rescue");
	  if(task == 0)
	    pf=fopen(buf,"w");
	  else
	    pf=fopen(buf,"a");
	  
	  for(i=0; i<domain_info[task].Nparts_per_node; i++)
	    fprintf(pf,"%d\n",Particle[i].Cluster_ID);
	  
	  fclose(pf);
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
      
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  for(i=0; i<NCLUSTERS; i++)
    Halos[i].NDomain_particles_absolute = Halos[i].NDomain_particles;
        
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  /////////////////////////////////////////////////////////////////////////////////
  /*      Sending Particle information from the other machines to the 0 task     */
  /////////////////////////////////////////////////////////////////////////////////
  
  for(i=1; i<Number_Processes; i++)
    {
      
      if(task == i)
	{
	  
	  /* Sending information of Halos */
	  
	  Parts_in_Domain_per_node = 0;
	  for(j=0; j<NCLUSTERS; j++)
	    Parts_in_Domain_per_node += Halos[j].NDomain_particles; 
	  
	  Info_domains = (int *) malloc((size_t) (Parts_in_Domain_per_node + 3*NCLUSTERS)*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      
	      Info_domains[k] = Halos[j].IDcluster;
	      k++;
	      Info_domains[k] = Halos[j].NDomain_particles;
	      k++;
	      Info_domains[k] = Halos[j].NDomain_particles_absolute;
	      k++;
	      
	      for(l=0; l<Halos[j].NDomain_particles; l++)
		{
		  Info_domains[k] = Halos[j].Domain_particles[l];
		  k++;
		}
	      
	    }
	  
	  sendbuff = Parts_in_Domain_per_node + 3*NCLUSTERS;
	  MPI_Ssend(&sendbuff,1,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  MPI_Ssend(&Info_domains[0],sendbuff,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  
	  free(Info_domains);
	  
	}
      
    }
  
  /* Now I have to collect the data of particles from the other tasks */
  
  if(task == 0)
    {
      
      for(i=1; i<Number_Processes; i++)
	{
	  
	  /* Receiving data of halos */
	  
	  MPI_Recv(&sendbuff,1,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  
	  printf("Receiving info of halos from %d %d\n",i,sendbuff); fflush(stdout);
	  
	  Info_domains = (int *) malloc((size_t) sendbuff*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  MPI_Recv(&Info_domains[0],sendbuff,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  printf("recibidos\n"); fflush(stdout);
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      ihalo = Info_domains[k];
	      k++;
	      
	      Old_NDomain_particles = Halos[ihalo].NDomain_particles;
	      NDomain_parts_per_proc = Info_domains[k];
	      Halos[ihalo].NDomain_particles += Info_domains[k];
	      k++;
	      
	      Old_NDomain_particles_absolute = Halos[ihalo].NDomain_particles_absolute;
	      NDomain_parts_per_proc_absolute = Info_domains[k];
	      Halos[ihalo].NDomain_particles_absolute += Info_domains[k];
	      k++;
	      
	      Halos[ihalo].Domain_particles = realloc(Halos[ihalo].Domain_particles,(size_t) Halos[ihalo].NDomain_particles*sizeof(int));
	      
	      m = Old_NDomain_particles;
	      for(l=0; l<NDomain_parts_per_proc; l++)
		{
		  Halos[ihalo].Domain_particles[m] = Info_domains[k];
		  k++;
		  m++;
		}
	    }
	  
	  free(Info_domains);
	  
	}
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0) 
    {
      
      aux_pf=fopen("Halo_Catalog.dat","w");
      fprintf(aux_pf,"#X\t\tY\t\tZ\t\tMvir\t\tRvir\t\tNmem\tNdom\tNdom_abs\n");
      for(i=0; i<NCLUSTERS; i++)
	fprintf(aux_pf,"%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%d\t%d\t%d\n",Halos[i].pos[0],Halos[i].pos[1],Halos[i].pos[2],
		Halos[i].Mvir,Halos[i].Rvir,Halos[i].Nmembers,Halos[i].NDomain_particles,
		Halos[i].NDomain_particles_absolute);
      fclose(aux_pf);
      
      printf("Writing info files\n");
      
      write_rescue(infile);
      
      printf("all is done!\n"); fflush(stdout);
      
    }
  

  MPI_Barrier(MPI_COMM_WORLD);
  free(Particle);
  free(Halos);
  
  MPI_Finalize();
  
  return 0;
  
}
