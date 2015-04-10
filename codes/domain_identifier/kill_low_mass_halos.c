int kill_low_mass_halos(float Mth,int i)
{
  
  if(Halos[i].Mvir < Mth)
    {
      printf("killing halo %d with Mvir %g because Mth %g\n",i,Halos[i].Mvir,Mth);
      Halos[i].Nmembers = 0;
      Halos[i].Mvir = 0;
      Halos[i].Rvir = 0;
      Halos[i].Nvir = 0;
      Halos[i].NDomain_particles = 0;
      Halos[i].NDomain_particles_absolute = 0;
    }
  
  /* if Mvir >= Mth, Removing particles out of Mvir */
  
  if(Halos[i].Nmembers > 0)
    {
      Halos[i].Nmembers = Halos[i].Nvir;
      Halos[i].Halo_particles = realloc(Halos[i].Halo_particles, (size_t) Halos[i].Nmembers*sizeof(int));
    }
  
  return 0;
  
}
