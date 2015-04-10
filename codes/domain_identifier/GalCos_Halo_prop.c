#include "GalCos_variables.h"

void periodic_boundary_corrections(int IDhalo)
{
    
  int xcounter,ycounter,zcounter,i,ipart,*indexpot=NULL;
  int NPARTICLES,HaloIDCenter,k,dk,j,jpart;
  float Xpos,Ypos,Zpos,x,y,z,Plenght,dist;
  float *potencial=NULL;
  
  
  /*********************************************************/
  /* CHOOSING A CLOSE CANDIDATE TO THE CENTER OF THA HALO */
  
  NPARTICLES = 20;
  
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
  
  for(i=0; i<NPARTICLES; i++)
    {
      potencial[i] = 0;
      indexpot[i] = 0;
    }
  
  dk = (int)  (Halos[IDhalo].Nmembers/20);
  k=0;
  for(i=0; i<NPARTICLES; i++)
    {
      
      ipart = Halos[IDhalo].Halo_particles[k];
      Particle[ipart].EP = 0;
      potencial[i] = 0;
      
      for(j=(Halos[IDhalo].Nmembers-1); j>=0; j--)
	{
	  
	  jpart = Halos[IDhalo].Halo_particles[j];
	  
	  if(ipart != jpart)
	    {
	      dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
			      Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
	      
	      Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*(Particle[jpart].mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT);
	    }
        }
      
      potencial[i] = Particle[ipart].EP;
      indexpot[i] = ipart;
      Particle[ipart].EP = 0.0;
      k = k+dk;
      
    }
  
  gsl_fisort(NPARTICLES,potencial,indexpot);
  
  HaloIDCenter = indexpot[0];
  
  free(potencial);
  free(indexpot);
  
  Xpos = Particle[HaloIDCenter].pos[0];
  Ypos = Particle[HaloIDCenter].pos[1];
  Zpos = Particle[HaloIDCenter].pos[2];
  
  Plenght = pow(0.5*BoxSize,2);
  
  for(i=0; i<Halos[IDhalo].Nmembers; i++)
    {
      ipart = Halos[IDhalo].Halo_particles[i];
      
      x = Particle[ipart].pos[0];
      if(pow(Xpos-x,2) > Plenght)
	{
	  
	  if(Xpos > x)
	    Particle[ipart].pos[0] = Particle[ipart].pos[0] + BoxSize;
	  else
	    Particle[ipart].pos[0] = Particle[ipart].pos[0] - BoxSize;
	  
	}
      
      y = Particle[ipart].pos[1];
      if(pow(Ypos-y,2) > Plenght)
	{
	  
	  if(Ypos > y)
	    Particle[ipart].pos[1] = Particle[ipart].pos[1] + BoxSize;
	  else
	    Particle[ipart].pos[1] = Particle[ipart].pos[1] - BoxSize;
	  
	}
      
      z = Particle[ipart].pos[2];
      if(pow(Zpos-z,2) > Plenght)
	{
	  
	  if(Zpos > z)
	    Particle[ipart].pos[2] = Particle[ipart].pos[2] + BoxSize;
	  else
	    Particle[ipart].pos[2] = Particle[ipart].pos[2] - BoxSize;
	  
	}
      
    }
  
}


void Halo_Mvir_Rvir(int ihalo)
{
  
  float RhoCrit,TotalMass,density,M,*radius;
  int k,ipart,*counter_radius,FLAG,istart;
  
  radius = (float *) malloc((size_t) Halos[ihalo].Nmembers*sizeof(float));
  counter_radius = (int *) malloc((size_t) Halos[ihalo].Nmembers*sizeof(int));
  
  M = 0.0;
  for(k=0; k<Halos[ihalo].Nmembers; k++)
    {
      ipart = Halos[ihalo].Halo_particles[k];
      radius[k] = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],0.0,0.0,0.0);
      counter_radius[k] = ipart;
      M += PARTMASS;
    }
  
  gsl_fisort(Halos[ihalo].Nmembers,radius,counter_radius);
  
  /// The particles in the halo are now indexed in growing order of
  /// distances
  
  for(k=0; k<Halos[ihalo].Nmembers; k++)
    Halos[ihalo].Halo_particles[k] = counter_radius[k];
  
  RhoCrit = Virial_criterion(); // Delta_c times rho_c
  
  FLAG = 0;
  istart = 1;
  TotalMass = istart*PARTMASS;
  for(k=istart; k<Halos[ihalo].Nmembers; k++)
    {
      
      TotalMass += PARTMASS;
      density = (3.0*TotalMass)/(4.0*M_PI*radius[k]*radius[k]*radius[k]);
      
      if( density <= RhoCrit )
	{
	  Halos[ihalo].Rvir = radius[k];
	  Halos[ihalo].Mvir = TotalMass;
	  Halos[ihalo].Nvir = k;
	  FLAG = 1;
	  break;
	}
      
    }
  
  if(FLAG == 0)
    {
      Halos[ihalo].Mvir = M;
      Halos[ihalo].Rvir = radius[Halos[ihalo].Nmembers-1];
    }
  
  Halos[ihalo].mass = M;
  
  free(radius);
  free(counter_radius);
  
}
