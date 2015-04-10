#include "GalCos_variables.h"
#include "BH_SPH_density.c"

int GalCos_domain(int ipart)
{
  
  int i,counter;
  float dist,mindist;
  
  counter = EMPTY_FLAG;
  mindist = 10000000.0*BoxSize;
  
  for(i=0; i<NCLUSTERS; i++)
    {
      
      
      //dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
      //Halos[i].pos[0],Halos[i].pos[1],Halos[i].pos[2]);
      
      /* For Periodic Boundary Corrections */
      dist = periodic_distance(Halos[i].pos,Particle[ipart].pos);
      
      dist = dist/Halos[i].Rvir;
      
      if(dist < mindist)
	{
	  mindist = dist;
	  counter = i;
	}
      
    }
  
  Particle[ipart].Cluster_ID = counter;
  
  Halos[counter].NDomain_particles++;
  Halos[counter].Domain_particles = realloc(Halos[counter].Domain_particles,(size_t) Halos[counter].NDomain_particles*sizeof(int));
  Halos[counter].Domain_particles_index = realloc(Halos[counter].Domain_particles_index,(size_t) Halos[counter].NDomain_particles*sizeof(int));
  
  Halos[counter].Domain_particles[Halos[counter].NDomain_particles-1] = Particle[ipart].Oid;
  Halos[counter].Domain_particles_index[Halos[counter].NDomain_particles-1] = Particle[ipart].id;
  
  return 0;
  
}

/* if is repeated return 1, if is not repeated returns 0*/
int is_not_repeated(int *array, int sizearray, int input)
{
  
  int i,FLAG=0;
  
  for(i=0; i<sizearray; i++)
    {
      if(array[i]  == input)
	{
	  FLAG=1;
	  break;
	}
    }
  
  return FLAG;
  
}


int GalCos_envelope(void)
{
  
  int i,j,k,counter,ipart,*SPH_NGB,kpart;
  float mindist,pos[3],ppos[3],x,y,z,Plenght;
  float Xpos,Ypos,Zpos,dx,dy,dz,*SPH_NGB_DISTANCES;
  
  TREENODEPTR root = NULL;
  
  SPH_NGB = (int *) malloc((size_t) NGB_MAX*sizeof(int));
  SPH_NGB_DISTANCES = (float *) malloc((size_t) NGB_MAX*sizeof(float));
  
  root = Build_tree_all(domain_info[task].Nparts_per_node);
  
  for(j=0; j<NCLUSTERS; j++)
    {
      
      Halos[j].NDomain_particles_absolute = Halos[j].NDomain_particles; // saving the information of domains
      
      for(i=0; i<Halos[j].NDomain_particles_absolute; i++)
	{
	  
	  ipart = Halos[j].Domain_particles_index[i];
	  
	  GalCos_SPH_NGB(ipart,SPH_NGB,SPH_NGB_DISTANCES,domain_info[task].Nparts_per_node,&root,NGB_MAX);
	  
	  for(k=0; k<NGB_MAX; k++)
	    {
	      kpart = SPH_NGB[k];
	      
	      if((Particle[kpart].Cluster_ID != j) && (is_not_repeated(Halos[j].Domain_particles,Halos[j].NDomain_particles,Particle[kpart].Oid) == 0))
		{
		  Halos[j].NDomain_particles++;
		  Halos[j].Domain_particles = realloc(Halos[j].Domain_particles,(size_t) Halos[j].NDomain_particles*sizeof(int));
		  Halos[j].Domain_particles[Halos[j].NDomain_particles-1] = Particle[kpart].Oid;
		}
	    }
	}
      
      if((j%1000)==0)
	printf("%d done halo %d\n",task,j);
      
    }
  
  return 0;
  
}
