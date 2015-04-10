#include "GalCos_variables.h"
#include "BHtree_pointers.c"

int GalCos_clusterBound(int *inputs,int NPARTICLES,int ihalo)
{
  
  int i,ipart,*indexpot=NULL,HaloIDCenter;
  float Xcm,Ycm,Zcm,VXcm,VYcm,VZcm,MASS;
  float *potencial=NULL;
  
  TREENODEPTR root=NULL;
  
  potencial = (float *) malloc((size_t) NPARTICLES*sizeof(float));
  if(potencial == NULL)
    {
      printf("Problem of memory allocation\n");
      exit(0);
    }
  
  indexpot = (int *) malloc((size_t) NPARTICLES*sizeof(int));
  if(indexpot == NULL)
    {
      printf("Problem of memory allocation\n");
      exit(0);
    }
  
  /////////////////////////////////////////////////////////////////////
  ///              COMPUTING GRAVITATIONAL POTENTIAL
  /////////////////////////////////////////////////////////////////////
  
  root = Build_tree(inputs,NPARTICLES);
  
  for(i=0; i<NPARTICLES; i++)
    {
      ipart = inputs[i];
      Particle[ipart].EP = 0.0;
      
      GalCos_tree_force(ipart,root);
      
      potencial[i] = Particle[ipart].EP;
      indexpot[i]  = ipart;
    }
  
  GalCos_Free_BHtree(&root);
  
  
  ///////////////////////////////////////////////////////////////////////////
  ///   LOCATING THE CENTER OF THE HALO AT THE POSITION OF THE MOST
  ///   BOUND PARTICLE THE VELOCITY OF THE HALO IS COMPUTED AS THE
  ///   VELOCITY OF THE CENTER OF MASS THE POSITION OF THE HALO IS
  ///   DEFINED AS THE POSITION OF THE MOST BOUND PARTICLE (LOWER
  ///   POTENTIAL ENERGY)
  ///////////////////////////////////////////////////////////////////////////
  
  
  gsl_fisort(NPARTICLES,potencial,indexpot);
  
  Xcm = Particle[indexpot[0]].pos[0];
  Ycm = Particle[indexpot[0]].pos[1];
  Zcm = Particle[indexpot[0]].pos[2];
  
  Halos[ihalo].pos[0] = Xcm;
  Halos[ihalo].pos[1] = Ycm;
  Halos[ihalo].pos[2] = Zcm;
  
  
  ////////////////////*********************
  
  
  VXcm = 0.0;
  VYcm = 0.0;
  VZcm = 0.0;
  MASS = 0.0;
  
  for(i=0; i<NPARTICLES; i++)
    {
      VXcm = VXcm + Particle[i].mass*Particle[i].vel[0];
      VYcm = VYcm + Particle[i].mass*Particle[i].vel[1];
      VZcm = VZcm + Particle[i].mass*Particle[i].vel[2];
      MASS = MASS + Particle[i].mass;
    }
  
  ///////////////////**********************
  
  VXcm = VXcm/MASS;
  VYcm = VYcm/MASS;
  VZcm = VZcm/MASS;
  
  Halos[ihalo].vel[0] = VXcm;
  Halos[ihalo].vel[1] = VYcm;
  Halos[ihalo].vel[2] = VZcm;
  
  HaloIDCenter = indexpot[0];
  
  free(potencial);
  free(indexpot);
  
  
  /* UPDATING THE POSITIONS OF THE PARTICLES RELATIVE TO THE CENTER OF
     THE HALO AND VELOCITIES RELATIVE TO VCM*/

  for(i=0; i<NPARTICLES; i++)
    {
      
      ipart = inputs[i];
      
      Particle[ipart].pos[0] = Particle[ipart].pos[0] - Xcm;
      Particle[ipart].pos[1] = Particle[ipart].pos[1] - Ycm;
      Particle[ipart].pos[2] = Particle[ipart].pos[2] - Zcm;
      
      Particle[ipart].vel[0] = Particle[ipart].vel[0] - VXcm;
      Particle[ipart].vel[1] = Particle[ipart].vel[1] - VYcm;
      Particle[ipart].vel[2] = Particle[ipart].vel[2] - VZcm;
      
    }
  
  return HaloIDCenter;
  
}
