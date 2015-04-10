#include<stdlib.h>
#include<stdio.h>
#include<math.h>

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


struct part
{
  float pos[3];
  float vel[3];
  //float mass;
  int id;
}*Part;


int Npart_Total;
struct gadget_head Gheader;


int read_head0(char *infile)
{
  
  int dummi,i;
  char label[4];
  FILE *fp_inp;
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(0);
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&label,sizeof(char),4,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&Gheader,sizeof(Gheader),1,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  /*
    REDSHIFT     = Gheader.redshift;
    OMEGA_MATTER = Gheader.Omega0;
    OMEGALAMBDA  = Gheader.OmegaLambda;
    HUBBLEPARAM  = Gheader.HubbleParam;
    COSMIC_TIME  = Gheader.time;
    BoxSize      = Gheader.BoxSize;
  */

  for(i=0; i<6; i++)
    printf(" * %d Particles of class %d\n",Gheader.npart[i],i); 
  
  printf("\n"); 
  
  for(i=0; i<6; i++)
    { 
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] != 0.0))
	printf(" * The mass of each particle is %d es %g\n",i,Gheader.mass[i]);
      
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] == 0.0))
	printf(" * There are individual mases for this particle set %d\n",i);
      
    }     
  
  printf("\n");
  
  printf(" * Frame's Time... %g\n",Gheader.time); 
  printf(" * Redshift... %g\n",Gheader.redshift);
  printf(" * Flagsfr... %d\n",Gheader.flag_sfr);
  printf(" * Flagfed... %d\n",Gheader.flag_feedback);
  
  printf("\n");
  
  Npart_Total = 0;
  for(i=0; i<6; i++)
    {
      printf(" * Header nall[%d] is: %d \n",i,Gheader.npartTotal[i]);
      Npart_Total += Gheader.npartTotal[i];
    }
  
  printf("\n");
  
  printf(" * Flagcool... %d\n",Gheader.flag_cooling);
  printf(" * numfiles... %d\n",Gheader.num_files);
  printf(" * Boxsize... %g\n",Gheader.BoxSize);
  printf(" * Omega0... %g\n",Gheader.Omega0);
  printf(" * OmageLa... %g\n",Gheader.OmegaLambda);
  printf(" * Hubbleparam... %g\n",Gheader.HubbleParam);
  printf(" * Npart_Total=%d\n",Npart_Total);
  
  Part = (struct part *) malloc((size_t) Npart_Total*sizeof(struct part));
  if(Part == NULL)
    {
      printf("No memory available for load dar particles (%s)\n",infile);
      exit(0);
    }
  
  printf("Using %f Mb in memory\n",Npart_Total*sizeof(struct part)/(1020*1024.0));

  return 0;
  
}


int main(int argc, char *argv[])
{
    
  int dummi,i,Npart_snap,global_acum,s,id,k;
  float pos[3],vel[3],aux[3];
  struct gadget_head header1;
  char label[4], buf[200], outfile[200];
  long int blksize_pos,blksize_vel,blksize_ids;
  int onlydm = atoi(argv[3]);
  
  float xmin,ymin,zmin,xmax,ymax,zmax;
  float vxmin,vymin,vzmin,vxmax,vymax,vzmax;
  

  FILE *fp_pos;
  FILE *fp_vel;
  FILE *fp_ids;
  FILE *fp_head;
  FILE *pf_out;
  
  char *snapbase = argv[1];
  int NSNAPS     = atoi(argv[2]);

  global_acum = 0;
  for(s=0; s<NSNAPS; s++)
    {
      //Filename of current file
      sprintf( buf,"%s.%d",snapbase,s );
      if( NSNAPS == 1 )
	  sprintf( buf,"%s",snapbase );

      if(s==0)
	read_head0(buf);

      if((fp_head=fopen(buf,"r"))==NULL)    
	{
	  printf("read_gadget cannot open %s",buf);
	  exit(0);
	}
      
      fp_pos=fopen(buf,"r");
      fp_vel=fopen(buf,"r");
      fp_ids=fopen(buf,"r");
      
      
      /* reading the header for this sub snap */
      fread(&dummi,sizeof(dummi),1,fp_head);
      fread(&label,sizeof(char),4,fp_head);
      fread(&dummi,sizeof(dummi),1,fp_head);
      fread(&dummi,sizeof(dummi),1,fp_head);

      fread(&dummi,sizeof(dummi),1,fp_head);
      fread(&header1,sizeof(header1),1,fp_head);
      fread(&dummi,sizeof(dummi),1,fp_head);
      
      fclose(fp_head);
      
      Npart_snap = 0;
      for(i=0; i<6; i++)
	{
	  printf(" * Header nall[%d] is: %d \n",i,header1.npart[i]);
	  Npart_snap += header1.npart[i];
	}
  
      
      //  Begining with the groups
      
      blksize_pos = 9*sizeof(int)  + 2*4*sizeof(char) + sizeof(header1);
      blksize_vel = 14*sizeof(int) + 3*4*sizeof(char) + sizeof(header1) + 3*Npart_snap*sizeof(float);
      blksize_ids = 19*sizeof(int) + 4*4*sizeof(char) + sizeof(header1) + 2*3*Npart_snap*sizeof(float);
      
      fseek(fp_pos,blksize_pos,SEEK_SET);
      fseek(fp_vel,blksize_vel,SEEK_SET);
      fseek(fp_ids,blksize_ids,SEEK_SET);
      
      for(i=0; i<Npart_snap ;i++)   
	{
	  fread(&pos[0],sizeof(float),3,fp_pos);
	  fread(&vel[0],sizeof(float),3,fp_vel);
	  fread(&id,sizeof(int),1,fp_ids);
 
	  Part[i].pos[0] = pos[0];
	  Part[i].pos[1] = pos[1];
	  Part[i].pos[2] = pos[2];
	  
	  Part[i].vel[0] = vel[0];
	  Part[i].vel[1] = vel[1];
	  Part[i].vel[2] = vel[2];
	  
	  Part[i].id = i;
	  
	  global_acum++;
	}
      
      fclose(fp_pos);
      fclose(fp_vel);
      fclose(fp_ids);
      
      printf("=====================================================\n");
      printf("   *-* End of reading from gadget file %s, %d \n",buf, global_acum);
      printf("=====================================================\n");
      
    }

  
  ///////////////////////////////////////////////////////////////////////
  /*                         writing in gadget 1 format !!             */
  ///////////////////////////////////////////////////////////////////////
  

  Gheader.npart[1] = Gheader.npartTotal[1];
  Gheader.num_files = 1;
  
  //Simulating an only dark-matter file
  if( onlydm == 1 ){
      for( i=0; i<6; i++ ){
// 	  printf( "******************MASS=%lf\n", Gheader.mass[i] );
	  Gheader.npart[i] = 0;
	  Gheader.npartTotal[i] = 0;}
    
      Gheader.npart[1] = Npart_Total;
      Gheader.npartTotal[1] = Npart_Total;}
  
  sprintf(outfile,"%s%s",snapbase,".FullSnap.gad1");
  
  printf("writing\n");

  if((pf_out=fopen(outfile,"w")) == NULL) {
    printf("cannot open %s\n",outfile);
    exit(0);
  }
  
  dummi = 256;
  
  fwrite(&dummi,sizeof(dummi),1,pf_out);
  fwrite(&Gheader,sizeof(Gheader),1,pf_out);
  fwrite(&dummi,sizeof(dummi),1,pf_out); // END GROUP  HEADER
  
  dummi = Npart_Total*3*sizeof(float);
  
  xmin=1.0e10;  ymin=1.0e10;  zmin=1.0e10; 
  xmax=-1.0e10; ymax=-1.0e10; zmax=-1.0e10;

  fwrite(&dummi,sizeof(dummi),1,pf_out); // INIT GROUP POSITIONS
  for(i=0; i<Npart_Total; i++){
    
    
    if(Part[i].pos[0] < xmin)
      xmin = Part[i].pos[0];

    if(Part[i].pos[1] < ymin)
      ymin = Part[i].pos[1];
    
    if(Part[i].pos[2] < zmin)
      zmin = Part[i].pos[2];
    
    if(Part[i].pos[0] > xmax)
      xmax = Part[i].pos[0];

    if(Part[i].pos[1] > ymax)
      ymax = Part[i].pos[1];
    
    if(Part[i].pos[2] > zmax)
      zmax = Part[i].pos[2];
    
    for(k=0; k<3; k++) 
      aux[k] = Part[i].pos[k];
    
    fwrite(aux,sizeof(float),3,pf_out);
  }
  fwrite(&dummi,sizeof(dummi),1,pf_out); // END GROUP POSITIONS
  
  printf(" xmin=%f ymin=%f zmin=%f\n",xmin,ymin,zmin);
  printf(" xmax=%f ymax=%f zmax=%f\n",xmax,ymax,zmax);
  
  vxmin=1.0e10;  vymin=1.0e10;  vzmin=1.0e10; 
  vxmax=-1.0e10; vymax=-1.0e10; vzmax=-1.0e10;
  
  fwrite(&dummi,sizeof(dummi),1,pf_out); // INIT GROUP VELOCITIES
  for(i=0; i<Npart_Total; i++){
    
    if(Part[i].vel[0] < vxmin)
      vxmin = Part[i].vel[0];

    if(Part[i].vel[1] < vymin)
      vymin = Part[i].vel[1];
    
    if(Part[i].vel[2] < vzmin)
      vzmin = Part[i].vel[2];

    if(Part[i].vel[0] > vxmax)
      vxmax = Part[i].vel[0];
    
    if(Part[i].vel[1] > vymax)
      vymax = Part[i].vel[1];
    
    if(Part[i].vel[2] > vzmax)
      vzmax = Part[i].vel[2];
    
    for(k=0; k<3; k++) 
      aux[k] = Part[i].vel[k];
    
    fwrite(aux,sizeof(float),3,pf_out);
  }
  fwrite(&dummi,sizeof(dummi),1,pf_out); //END GROUP VELOCITIES
  
  printf(" vxmin=%f vymin=%f vzmin=%f\n",vxmin,vymin,vzmin);
  printf(" vxmax=%f vymax=%f vzmax=%f\n",vxmax,vymax,vzmax);
  
  dummi = Npart_Total*sizeof(int);
  fwrite(&dummi,sizeof(dummi),1,pf_out); // INIT GROUP ID'S
  
  int idmin=2.0e9;
  int idmax=-2.0e9;
  
  for(i=0; i<Npart_Total; i++)
  {
      
      if(Part[i].id < idmin)
	  idmin = Part[i].id;
      
      if(Part[i].id > idmax)
	  idmax = Part[i].id;
      
      fwrite(&Part[i].id,sizeof(int),1,pf_out);
      
      if(i<0)
	  printf("Here is the problem!\n");
      
  }
  fwrite(&dummi,sizeof(dummi),1,pf_out); // END GROUP ID'S
  
  printf(" idmin=%d idmin=%d\n",idmin,idmin);
  
  fclose(pf_out);
  
  return 0;  
  
}
