#include<stdlib.h>
#include "GalCos_variables.h"

extern int DEAD_NODE;

struct queueNode *NGB_headPtr;
struct queueNode *NGB_tailPtr;

/* Looking for intersection between the search sphere and the node */

int GalCos_SPHintersect_node(float Lsearch,float *pos,float *posnode,float sizenode)
{
  
  int counter=0;
  float dist,Shalf,a,Npos[3],Plenght;
  
  dist = distance(pos[0],pos[1],pos[2],posnode[0],posnode[1],posnode[2]);
  Shalf = sizenode*0.5;
  a = 2*sqrt(2.0*Shalf*Shalf);
  
  if(dist <= (Lsearch + a))
    {
      counter=1;
      return counter;
    }

  /* Including corrections for periodic boundary conditions */
  
  Plenght = pow(0.5*BoxSize,2);
  
  Npos[0] = posnode[0];
  Npos[1] = posnode[1];
  Npos[2] = posnode[2];
  
  if(pow(posnode[0]-pos[0],2) > Plenght)
    {
      if(pos[0] > posnode[0])
	Npos[0] = posnode[0] + BoxSize;
      else
	Npos[0] = posnode[0] - BoxSize;
    }
  
  if(pow(posnode[1]-pos[1],2) > Plenght)
    {
      if(pos[1] > posnode[1])
	Npos[1] = posnode[1] + BoxSize;
      else
	Npos[1] = posnode[1] - BoxSize;
    }
  
  if(pow(posnode[2]-pos[2],2) > Plenght)
    {
      if(pos[2] > posnode[2])
	Npos[2] = posnode[2] + BoxSize;
      else
	Npos[2] = posnode[2] - BoxSize;
    }
  
  dist = distance(pos[0],pos[1],pos[2],Npos[0],Npos[1],Npos[2]);
    
  if(dist <= (Lsearch + a))
    {
      counter=1;
    }
  
  return counter;

}

/* This looks if the node is inside the sphere of search */

int GalCos_node_inside_node(float Lsearch,float *pos,float *posnode,float sizenode)
{
  
  int counter=0;
  float dist,Shalf,a,Npos[3],Plenght;
  
  dist = distance(pos[0],pos[1],pos[2],posnode[0],posnode[1],posnode[2]);
  Shalf = sizenode*0.5;
  a = 2*sqrt(2.0*Shalf*Shalf);
  
  if(Lsearch >= (dist+a))
    {
      counter=1;
      return counter;
    }
  
  if(sizenode >= Lsearch)
    {
      
      if( ((pos[0]+Lsearch) <= (posnode[0]+Shalf)) && ((pos[0]-Lsearch) >= (posnode[0]-Shalf)) )
	{
	  if( ((pos[1]+Lsearch) <= (posnode[1]+Shalf)) && ((pos[1]-Lsearch) >= (posnode[1]-Shalf)) )
	    {
	      if( ((pos[2]+Lsearch) <= (posnode[2]+Shalf)) && ((pos[2]-Lsearch) >= (posnode[2]-Shalf)) )
		{
		  counter=1;
		  return counter;
		}
	    }
	}
    }
  
  /* Including corrections for periodic boundary conditions */
  
  Plenght = pow(0.5*BoxSize,2);
  
  Npos[0] = posnode[0];
  Npos[1] = posnode[1];
  Npos[2] = posnode[2];
  
  if(pow(posnode[0]-pos[0],2) > Plenght)
    {
      if(pos[0] > posnode[0])
	Npos[0] = posnode[0] + BoxSize;
      else
	Npos[0] = posnode[0] - BoxSize;
    }
  
  if(pow(posnode[1]-pos[1],2) > Plenght)
    {
      if(pos[1] > posnode[1])
	Npos[1] = posnode[1] + BoxSize;
      else
	Npos[1] = posnode[1] - BoxSize;
    }
  
  if(pow(posnode[2]-pos[2],2) > Plenght)
    {
      if(pos[2] > posnode[2])
	Npos[2] = posnode[2] + BoxSize;
      else
	Npos[2] = posnode[2] - BoxSize;
    }
  
  dist = distance(pos[0],pos[1],pos[2],Npos[0],Npos[1],Npos[2]);
  
  if(Lsearch >= (dist+a))
    {
      counter=1;
      return counter;
    }
  
  if(sizenode >= Lsearch)
    {
      
      if( ((pos[0]+Lsearch) <= (Npos[0]+Shalf)) && ((pos[0]-Lsearch) >= (Npos[0]-Shalf)) )
	{
	  if( ((pos[1]+Lsearch) <= (Npos[1]+Shalf)) && ((pos[1]-Lsearch) >= (Npos[1]-Shalf)) )
	    {
	      if( ((pos[2]+Lsearch) <= (Npos[2]+Shalf)) && ((pos[2]-Lsearch) >= (Npos[2]-Shalf)) )
		{
		  counter=1;
		  return counter;
		}
	    }
	}
      
    }
  
  return counter;
  
}


void NGB_driver(int ipart,TREENODEPTR *p,float Lsearch)
{
  
  int flag_sup=0,i,jpart;
  float pos[3],Plenght,new_dist,Npos[3];
  
  if((*p) != NULL){
    
    pos[0] = Particle[ipart].pos[0];
    pos[1] = Particle[ipart].pos[1];
    pos[2] = Particle[ipart].pos[2];
    
    flag_sup = flag_sup + GalCos_SPHintersect_node(Lsearch,pos,(*p)->pos,(*p)->size);
    flag_sup = flag_sup + GalCos_node_inside_node(Lsearch,pos,(*p)->pos,(*p)->size);
    
    Plenght = pow(0.5*BoxSize,2);
    
    if(flag_sup > 0)
      {
	
	if((*p)->tag == 1){
	  
	  for(i=0; i<8; i++){
	    
	    //if ( ((*p)->NSon[i] >= 0) && ((*p)->NSon[i] < TREE_BUILD_NPARTICLES) )
	    if ( ((*p)->NSon[i] >= 0) && ((*p)->NSon[i] < domain_info[task].Nparts_per_node) )
	      {
		
		jpart = (*p)->NSon[i];
		if(jpart != ipart)
		  {
		    
		    /* Including corrections for periodic boundary conditions */
		    
		    Npos[0] = Particle[jpart].pos[0];
		    Npos[1] = Particle[jpart].pos[1];
		    Npos[2] = Particle[jpart].pos[2];
		    
		    if(pow(Particle[jpart].pos[0] - pos[0],2) > Plenght)
		      {
			if(pos[0] > Particle[jpart].pos[0])
			  Npos[0] = Particle[jpart].pos[0] + BoxSize;
			else
			  Npos[0] = Particle[jpart].pos[0] - BoxSize;
		      }
		    
		    if(pow(Particle[jpart].pos[1]-pos[1],2) > Plenght)
		      {
			if(pos[1] > Particle[jpart].pos[1])
			  Npos[1] = Particle[jpart].pos[1] + BoxSize;
			else
			  Npos[1] = Particle[jpart].pos[1] - BoxSize;
		      }
			  
		    if(pow(Particle[jpart].pos[2]-pos[2],2) > Plenght)
		      {
			if(pos[2] > Particle[jpart].pos[2])
			  Npos[2] = Particle[jpart].pos[2] + BoxSize;
			else
			  Npos[2] = Particle[jpart].pos[2] - BoxSize;
		      }
		    
		    new_dist = distance(pos[0],pos[1],pos[2],Npos[0],Npos[1],Npos[2]);
		    
		    if(new_dist <= Lsearch)
		      {
			Dencola(&NGB_headPtr,&NGB_tailPtr,jpart);
			encola(&NGB_cand,jpart);
			encola_float(&NGB_DISTcand,new_dist);
		      }
		    
		    
		  }
	      }
	    
	  }
	  
	  for(i=0; i<8; i++)
	    {
	      if((*p)->NSon[i] > TREE_BUILD_NPARTICLES)
		NGB_driver(ipart,&((*p)->sons[i]),Lsearch);
	      
	    }
	  
	}
	
      }
    
  }
  
}

int GalCos_SPH_NGB_Rfix(int ipart,int *Halo_particles,TREENODEPTR *p,float Lsearch)
{
  
  int i,iaux=0,Nmembers=0,ihalo;
  
  while(NGB_headPtr != NULL)
    iaux = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
  
  NGB_headPtr=NULL;
  NGB_tailPtr=NULL;
  
  NGB_driver(ipart,&(*p),Lsearch);
  Dsizecola(NGB_headPtr,&Nmembers);
  
  if(Nmembers > 0)
    {
      
      Halo_particles = (int *) malloc((size_t) Nmembers*sizeof(int));
      
      //for(i=0; i<Halos[ihalo].Nmembers; i++)
      for(i=0; i<Nmembers; i++)
	{
	  Halo_particles[i] = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
	}
      
    }
  
  while(NGB_headPtr != NULL)
    iaux = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
  
  NGB_headPtr=NULL;
  NGB_tailPtr=NULL;

  return Nmembers;  
}


void call_func(int ipart,TREENODEPTR *p,float Lsearch, int NGB_MAX)
{
  
  int Nmembers=0,iaux;
  
  NGB_driver(ipart,&(*p),Lsearch);
  Dsizecola(NGB_headPtr,&Nmembers);
  
  if(Nmembers < NGB_MAX)
    {
      
      Lsearch = Lsearch + Lsearch/2.0;
      
      while(NGB_headPtr != NULL)
	iaux = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
      
      NGB_headPtr=NULL;
      NGB_tailPtr=NULL;
      
      reinitialize_cola(&NGB_cand);
      reinitialize_cola_float(&NGB_DISTcand);
      
      call_func(ipart,&(*p),Lsearch,NGB_MAX);
      
    }
  
}



void GalCos_SPH_NGB(int ipart,int *SPH_NGB,float *SPH_NGB_DISTANCES,int NPARTICLES,TREENODEPTR *p,int NGB_MAX)
{
  
  int i,Nmembers=0,iaux;
  float N_dens,Lsearch;
  
  N_dens = 1.0*NPARTICLES/pow((*p)->size,3);
  Lsearch = 0.01*pow(1.0*NGB_MAX/N_dens,1.0/3.0);
  
  while(NGB_headPtr != NULL)
    iaux = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
  
  NGB_headPtr=NULL;
  NGB_tailPtr=NULL;
  
  call_func(ipart,&(*p),Lsearch,NGB_MAX);
  
  Dsizecola(NGB_headPtr,&Nmembers);
  gsl_fisort(Nmembers,NGB_DISTcand.inputs,NGB_cand.inputs);
  
  for(i=0; i<NGB_MAX; i++)
    {
      SPH_NGB[i] = Ddesencola(&NGB_headPtr,&NGB_tailPtr);
      SPH_NGB_DISTANCES[i] = NGB_DISTcand.inputs[i];
    }
  
  reinitialize_cola(&NGB_cand);
  reinitialize_cola_float(&NGB_DISTcand);
  
}

