#include "GalCos_variables.h"
//#include "GalCos_variables.c"

#include "cola_deitel.c"

int DEAD_NODE=-1;
int GLOBAL_COUNTER;
int MAX_ID_PARTS=0;
int *auxpointer=NULL,*copy_inputs=NULL,TREE_BUILD_NPARTICLES;

int get_node_center(TREENODEPTR father,float *pos,int i)
{
  
  /*
    Here is defined the numberin of nodes in the tree. starting with the
    lower left (i=0) growing in counter clockwise and from bottom to top
  */
  
  float Qsize;
  
  Qsize = 0.25*father->size;
  
  if(i == 0)
    {
      pos[0] = father->pos[0] - Qsize;
      pos[1] = father->pos[1] - Qsize;
      pos[2] = father->pos[2] - Qsize;
    }
  
  if(i == 1)
    {
      pos[0] = father->pos[0] + Qsize;
      pos[1] = father->pos[1] - Qsize;
      pos[2] = father->pos[2] - Qsize;
    }
  
  if(i == 2)
    {
      pos[0] = father->pos[0] + Qsize;
      pos[1] = father->pos[1] + Qsize;
      pos[2] = father->pos[2] - Qsize;
    }
  
  if(i == 3)
    {
      pos[0] = father->pos[0] - Qsize;
      pos[1] = father->pos[1] + Qsize;
      pos[2] = father->pos[2] - Qsize;
    }
  
  if(i == 4)
    {
      pos[0] = father->pos[0] - Qsize;
      pos[1] = father->pos[1] - Qsize;
      pos[2] = father->pos[2] + Qsize;
    }
  
  if(i == 5)
    {
      pos[0] = father->pos[0] + Qsize;
      pos[1] = father->pos[1] - Qsize;
      pos[2] = father->pos[2] + Qsize;
    }
  
  if(i == 6)
    {
      pos[0] = father->pos[0] + Qsize;
      pos[1] = father->pos[1] + Qsize;
      pos[2] = father->pos[2] + Qsize;
    }
  
  if(i == 7)
    {
      pos[0] = father->pos[0] - Qsize;
      pos[1] = father->pos[1] + Qsize;
      pos[2] = father->pos[2] + Qsize;
    }
  
  return 0;
}


int get_number_of_particles(float *pos,float size,int *IDpart,struct queueNode **FheadPtr,
			    struct queueNode **FtailPtr,int *Nimputs,
			    struct queueNode **SheadPtr,struct queueNode **StailPtr)
{
  
  float xmin,xmax,ymin,ymax,zmin,zmax,Hbox;
  int counter,i,ipart,FLAG,cn,NpartsInBox;
  
  Hbox = 0.5*size;
  
  xmin = pos[0] - Hbox;
  xmax = pos[0] + Hbox;
  
  ymin = pos[1] - Hbox;
  ymax = pos[1] + Hbox;
  
  zmin = pos[2] - Hbox;
  zmax = pos[2] + Hbox;
  
  counter = 0;
  cn = 0;
  
  NpartsInBox = *Nimputs;
  
  for(i=0; i<NpartsInBox; i++)
    {
      
      ipart = Ddesencola(FheadPtr,FtailPtr);
      FLAG = 0;
      
      if( (Particle[ipart].pos[0] >= xmin) && (Particle[ipart].pos[0] < xmax) )
	{
	  if( (Particle[ipart].pos[1] >= ymin) && (Particle[ipart].pos[1] < ymax) )
	    {
	      if( (Particle[ipart].pos[2] >= zmin) && (Particle[ipart].pos[2] < zmax) )
		{
		  
		  *IDpart = Particle[ipart].id;
		  FLAG=1;
		  
		  Dencola(SheadPtr,StailPtr,ipart);
		  counter++;
		  
		  if((ipart != Particle[ipart].id) || (ipart != *IDpart))
		    {
		      printf("mmm un problema!!\n");
		      exit(0);
		    }
		}
	    }
	}
      
      
      if(FLAG == 0)
	{
	  Dencola(FheadPtr,FtailPtr,ipart);
	  cn++;
	}
      
    }
  
  if((counter+cn) != NpartsInBox)
    {
      printf("Mistake building the tree...\n");
      exit(0);
    }

  *Nimputs=cn;
  
  return counter;
  
}


void alloca_node(TREENODEPTR *p,struct ptnode **father,float size,float *pos,int ith)
{
  
  int i,iparticle,Nparticles,iaux;
  float center[3],sizenode;
  
  if((*p) == NULL)
    {
      
      (*p) = (struct ptnode *) malloc((size_t) sizeof(struct ptnode));
      
      if((*p) != NULL)
	{
	  
	  (*p)->sons = NULL;
	  (*p)->headPtr=NULL;
	  (*p)->tailPtr=NULL;
	  
	  if(GLOBAL_COUNTER == MAX_ID_PARTS)
	    {
	      
	      (*p)->NumParts = TREE_BUILD_NPARTICLES;
	      
	      for(i=0; i<TREE_BUILD_NPARTICLES; i++)
		Dencola(&((*p)->headPtr),&((*p)->tailPtr),copy_inputs[i]);
	      
	      free(copy_inputs);
	      copy_inputs=NULL;
	      
	    }
	  else
	    {
	      
	      (*p)->NumParts = get_number_of_particles(pos,size,&iparticle,&((*father)->headPtr),
						       &((*father)->tailPtr),&((*father)->NumParts),
						       &((*p)->headPtr),&((*p)->tailPtr));
	    }
	  
	  (*p)->pos[0] = pos[0];
	  (*p)->pos[1] = pos[1];
	  (*p)->pos[2] = pos[2];
	  (*p)->size = size;
	  
	  Nparticles = (*p)->NumParts;
	  
	  /* Empty node */
	  if(Nparticles == 0)
	    {
	      (*father)->NSon[ith] = DEAD_NODE;
	      
	      while((*p)->headPtr != NULL)
		{
		  iaux = Ddesencola(&(*p)->headPtr,&(*p)->tailPtr);
		}
	      
	      free((*father)->sons[ith]);
	      (*father)->sons[ith]=NULL;
	      
	      return;
	    }
	  
	  /* leaf node... PARTICLE */
	  if(Nparticles == 1)
	    {
	      
	      (*father)->NSon[ith] = iparticle;
	      
	      while((*p)->headPtr != NULL)
		{
		  iaux = Ddesencola(&(*p)->headPtr,&(*p)->tailPtr);
		}
	      
	      free((*father)->sons[ith]);
	      (*father)->sons[ith]=NULL;
	      	      
	      return;
	    }
	  
	  /* Twig node */
	  if(Nparticles > 1)
	    {
	      (*p)->index = GLOBAL_COUNTER;
	      GLOBAL_COUNTER++;
	      
	      if((*p)->index > MAX_ID_PARTS)
		(*father)->NSon[ith] = (*p)->index;
	      
	      (*p)->tag = 1;
	      (*p)->mass = 0.0;
	      
	      (*p)->sons = (struct ptnode **) malloc((size_t) 8*sizeof(struct ptnode *));
	      if((*p)->sons == NULL)
		{
		  printf("No memory available to allocate tree\n");
		  exit(0);
		}
	      
	      for(i=0; i<8; i++)
		(*p)->sons[i]=NULL;
	      	      
	      for(i=0; i<8; i++)
		{
		  get_node_center(*p,center,i);
		  sizenode = 0.5*(*p)->size;
		  alloca_node(&((*p)->sons[i]),&(*p),sizenode,center,i);
		}
	      
	    }
	  
	}
      else
	{
	  printf("No memory available to allocate pointer\n");
	  exit(0);	    
	}
      
    }
  
}

void GalCos_tree_mass(TREENODEPTR *p)
{
  
  int i;
  
  if((*p) != NULL)
    {
      
      for(i=0; i<8; i++)
	{
	  if((*p)->NSon[i] > MAX_ID_PARTS)
	    GalCos_tree_mass(&((*p)->sons[i]));
	}
      
      //printf("mass %d\n",(*p)->index);
      
      if((*p)->tag == 1)
	{
	  
	  if((*p)->mass == 0)
	    {
	      for(i=0; i<8; i++)
		{
		  
		  if ( ((*p)->NSon[i] < MAX_ID_PARTS) && ((*p)->NSon[i] >= 0) ) // particle
		    (*p)->mass = (*p)->mass + Particle[(*p)->NSon[i]].mass;
		  else if((*p)->NSon[i] > MAX_ID_PARTS) // node
		    (*p)->mass = (*p)->mass + (*p)->sons[i]->mass;
		}
	    }
	  
	}
      
    }
  
}


void GalCos_tree_CM(TREENODEPTR *p)
{
  
  int i,ipart;
  
  if((*p) != NULL)
    {
      
      for(i=0; i<8; i++)
	{
	  if((*p)->NSon[i] > MAX_ID_PARTS)
	    GalCos_tree_CM(&((*p)->sons[i]));
	}
      
      if((*p)->tag == 1)
	{
	  
	  (*p)->pos[0] = 0.0;
	  (*p)->pos[1] = 0.0;
	  (*p)->pos[2] = 0.0;
	  
	  for(i=0; i<8; i++)
	    {
	      
	      if ( ((*p)->NSon[i] >= 0) && ((*p)->NSon[i] < MAX_ID_PARTS) ) // particle
		{
		  ipart=(*p)->NSon[i];
		  (*p)->pos[0] = (*p)->pos[0] + Particle[ipart].pos[0]*Particle[ipart].mass;
		  (*p)->pos[1] = (*p)->pos[1] + Particle[ipart].pos[1]*Particle[ipart].mass;
		  (*p)->pos[2] = (*p)->pos[2] + Particle[ipart].pos[2]*Particle[ipart].mass;
		}
	      else if((*p)->NSon[i] >= MAX_ID_PARTS) // node
		{
		  (*p)->pos[0] = (*p)->pos[0] + ((*p)->sons[i]->pos[0])*((*p)->sons[i]->mass);
		  (*p)->pos[1] = (*p)->pos[1] + ((*p)->sons[i]->pos[1])*((*p)->sons[i]->mass);
		  (*p)->pos[2] = (*p)->pos[2] + ((*p)->sons[i]->pos[2])*((*p)->sons[i]->mass);
		}
	    }
	  
	  (*p)->pos[0] = (*p)->pos[0]/(*p)->mass;
	  (*p)->pos[1] = (*p)->pos[1]/(*p)->mass;
	  (*p)->pos[2] = (*p)->pos[2]/(*p)->mass;
	  
	}
      
    }
  
}


void GalCos_Free_BHtree(TREENODEPTR *p)
{
  
  int i,iaux;
  
  if(*p != NULL)
    {
      
      
      if((*p)->tag == 1)
	{
	  
	  for(i=0; i<8; i++)
	    {
	      if((*p)->NSon[i] > MAX_ID_PARTS) 
		GalCos_Free_BHtree(&((*p)->sons[i]));
	    }
	  
	}
      
      while((*p)->headPtr != NULL)
	{
	  iaux = Ddesencola(&(*p)->headPtr,&(*p)->tailPtr);
	}
                        
      if((*p)->sons != NULL)
	free((*p)->sons);
      
      free((*p));
    }
  
  if(auxpointer != NULL)
    {
      free(auxpointer);
      auxpointer=NULL;
    }
  
}


/* 
   Computes the distance between two points but doing periodic
   boundary corrections with center in posa 
*/

float periodic_distance(float *posa,float *posb)
{
  
  float dist,Npos[3],Plenght;
  
  dist = distance(posa[0],posa[1],posa[2],posb[0],posb[1],posb[2]);
  
  if(dist <= 0.5*BoxSize)
    return dist;
  
  /* Including corrections for periodic boundary conditions */
  
  Plenght = pow(0.5*BoxSize,2);
  
  Npos[0] = posb[0];
  Npos[1] = posb[1];
  Npos[2] = posb[2];
  
  if(pow(posb[0]-posa[0],2) > Plenght)
    {
      if(posa[0] > posb[0])
	Npos[0] = posb[0] + BoxSize;
      else
	Npos[0] = posb[0] - BoxSize;
    }
  
  if(pow(posb[1]-posa[1],2) > Plenght)
    {
      if(posa[1] > posb[1])
	Npos[1] = posb[1] + BoxSize;
      else
	Npos[1] = posb[1] - BoxSize;
    }
  
  if(pow(posb[2]-posa[2],2) > Plenght)
    {
      if(posa[2] > posb[2])
	Npos[2] = posb[2] + BoxSize;
      else
	Npos[2] = posb[2] - BoxSize;
    }
  
  dist = distance(posa[0],posa[1],posa[2],Npos[0],Npos[1],Npos[2]);
  
  return dist;
  
}


void GalCos_tree_force(int ipart,TREENODEPTR p)
{
  
  float dist;
  int i,jpart;
  
  if(p != NULL)
    {
      
      if(p->tag == 1)
	{
	  
	  dist = periodic_distance(Particle[ipart].pos,p->pos);
	  
	  //dist=distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
	  //p->pos[0],p->pos[1],p->pos[2]);
	  
	  if((p->size/dist) <= BH_OPENING ) // Particle-Node
	    {
	      Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*((p->mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT));
	    }
	  else /* Opening the node -- walk across subnodes inside this node */
	    {
	      
	      /* 
		 First of all I will add the force from PP interaction
		 from the individual particles inside this node
	      */
	      
	      for(i=0; i<8; i++)
		{
		  if ( (p->NSon[i] >= 0) && (p->NSon[i] < MAX_ID_PARTS) )
		    {
		      jpart = p->NSon[i];
		      dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
				      Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
		      Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*((Particle[jpart].mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT));
		    }
		}
	      
	      for(i=0; i<8; i++)
		{
		  if(p->NSon[i] > MAX_ID_PARTS)
		    GalCos_tree_force(ipart,p->sons[i]);
		}
	      
	    }
	  
	}
       
    }
  
}

void GalCos_tree_print(TREENODEPTR p)
{
  
  int i;
  
  if(p != NULL)
    {
      
      if(p->tag == 1)
	{
	  printf("%f %f\n",p->NumParts*Particle[0].mass,p->mass);
	}
      
      
      if(p->NumParts > 1)
	{
	  
	  for(i=0; i<8; i++)
	    {
	      GalCos_tree_print(p->sons[i]);
	    }
	}
      
      /*
	float xmin,xmax,ymin,ymax,zmin,zmax,Hbox;
	FILE *pf;
	
	pf=fopen("plotar.gpl","a");
	
	Hbox=p->size/2.0;
	
	xmin=p->pos[0] - Hbox;
	xmax=p->pos[0] + Hbox;
	
	ymin=p->pos[1] - Hbox;
	ymax=p->pos[1] + Hbox;
	
	zmin=p->pos[2] - Hbox;
	zmax=p->pos[2] + Hbox;
	
	//fprintf(pf,"set ticslevel 0\n");
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmax,ymin,zmin);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmin,ymax,zmin);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmin,ymin,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmin,xmax,ymax,zmin);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmin,xmax,ymin,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymax,zmin,xmax,ymax,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmin,xmin,ymax,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmin,xmax,ymax,zmin);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmax,xmax,ymin,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmax,xmin,ymax,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmax,xmax,ymax,zmax);
	fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmax,xmax,ymax,zmax);
	//fprintf(pf,"pause -1\n");
	
	fclose(pf);
      */
    }
    
}


TREENODEPTR Build_tree(int *inputs, int NPARTICLES)
{
    
  float BOXSIZE,pos[3],xmin,xmax,ymin,ymax,zmin,zmax; 
  float xsize,ysize,zsize,BoxXCenter,BoxYCenter,BoxZCenter;
  int i,ipart;
  TREENODEPTR root=NULL;
    
  GLOBAL_COUNTER=0;
  TREE_BUILD_NPARTICLES=0;
  DEAD_NODE=-1;
  
  copy_inputs = (int *) malloc((size_t) NPARTICLES*sizeof(int));
  
  xmin = Particle[inputs[0]].pos[0]; ymin=Particle[inputs[0]].pos[1]; zmin=Particle[inputs[0]].pos[2];
  xmax = Particle[inputs[0]].pos[0]; ymax=Particle[inputs[0]].pos[1]; zmax=Particle[inputs[0]].pos[2];
  
  for(i=0; i<NPARTICLES; i++)
    {
      
      ipart = inputs[i];
      
      // min
      if(Particle[ipart].pos[0] < xmin)
	xmin = Particle[ipart].pos[0];
      
      if(Particle[ipart].pos[1] < ymin)
	ymin = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] < zmin)
	zmin = Particle[ipart].pos[2];
      
      // max
      if(Particle[ipart].pos[0] > xmax)
	xmax = Particle[ipart].pos[0];
      
      if(Particle[ipart].pos[1] > ymax)
	ymax = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] > zmax)
	zmax = Particle[ipart].pos[2];
      
      copy_inputs[i] = ipart;
      
    }
  
  xsize = (xmax - xmin);
  ysize = (ymax - ymin);
  zsize = (zmax - zmin);
  
  BoxXCenter = xmin + xsize*0.5;
  BoxYCenter = ymin + ysize*0.5;
  BoxZCenter = zmin + zsize*0.5;
  
  BOXSIZE = 0;
  
  if(xsize > BOXSIZE)
    BOXSIZE = xsize;
  if(ysize > BOXSIZE)
    BOXSIZE = ysize;
  if(zsize > BOXSIZE)
    BOXSIZE = zsize;
  
  if(BOXSIZE == 0)
    {
      printf("Problem computing the size root cell (N=%d) of the tree %g %g %g\n",
	     NPARTICLES,xsize,ysize,zsize); fflush(stdout);
      exit(0);
    }
  
  BOXSIZE = BOXSIZE + 0.01*BOXSIZE;
  
  pos[0] = BoxXCenter;
  pos[1] = BoxYCenter;
  pos[2] = BoxZCenter;
  
  
  /* INFORMATION */
  printf("\n");  fflush(stdout);
  printf("xmin=%g xmax=%g\n",xmin,xmax);  fflush(stdout);
  printf("ymin=%g ymax=%g\n",ymin,ymax); fflush(stdout);
  printf("zmin=%g zmax=%g\n",zmin,zmax); fflush(stdout);
  printf("Boxsize=%g  Nparts=%d\n",BOXSIZE,NPARTICLES); fflush(stdout);
  printf("x=%g y=%g z=%g\n",pos[0],pos[1],pos[2]); fflush(stdout);
  
  TREE_BUILD_NPARTICLES = NPARTICLES;
  DEAD_NODE = -1;
  
  MAX_ID_PARTS = 0;
  
  if(TREE_BUILD_NPARTICLES > domain_info[task].Nparts_per_node)
    MAX_ID_PARTS = TREE_BUILD_NPARTICLES;
  else
    MAX_ID_PARTS = domain_info[task].Nparts_per_node;
  
  GLOBAL_COUNTER = MAX_ID_PARTS;
  
  alloca_node(&root,NULL,BOXSIZE,pos,0);
  
  printf("Built %d tree nodes\n",GLOBAL_COUNTER-domain_info[task].Nparts_per_node);
  
  GalCos_tree_mass(&root);
  printf("done mass\n");
  GalCos_tree_CM(&root);
  printf("done cm\n");
  
  GLOBAL_COUNTER=0;
  
  return root;
  
}


TREENODEPTR Build_tree_all(int NPARTICLES)
{
    
    float BOXSIZE,pos[3],xmin,xmax,ymin,ymax,zmin,zmax; 
    float xsize,ysize,zsize,BoxXCenter,BoxYCenter,BoxZCenter;
    int i,ipart;
    TREENODEPTR root=NULL;
        
    copy_inputs = (int *) malloc((size_t) NPARTICLES*sizeof(int));
    
    xmin = Particle[0].pos[0]; ymin=Particle[0].pos[1]; zmin=Particle[0].pos[2];
    xmax = Particle[0].pos[0]; ymax=Particle[0].pos[1]; zmax=Particle[0].pos[2];
    
    for(i=0; i<NPARTICLES; i++)
    {
      
      ipart = i;
      
      // min
      if(Particle[ipart].pos[0] < xmin)
	xmin = Particle[ipart].pos[0];
      
      if(Particle[ipart].pos[1] < ymin)
	ymin = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] < zmin)
	zmin = Particle[ipart].pos[2];
      
      // max
      if(Particle[ipart].pos[0] > xmax)
	xmax = Particle[ipart].pos[0];
      
      if(Particle[ipart].pos[1] > ymax)
	ymax = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] > zmax)
	zmax = Particle[ipart].pos[2];
      
      copy_inputs[i] = ipart;
    }
    
    xsize = (xmax - xmin);
    ysize = (ymax - ymin);
    zsize = (zmax - zmin);
    
    BoxXCenter = xmin + xsize*0.5;
    BoxYCenter = ymin + ysize*0.5;
    BoxZCenter = zmin + zsize*0.5;
    
    BOXSIZE = 0;
    
    if(xsize > BOXSIZE)
	BOXSIZE = xsize;
    if(ysize > BOXSIZE)
      BOXSIZE = ysize;
    if(zsize > BOXSIZE)
      BOXSIZE = zsize;
    
    if(BOXSIZE == 0)
      {
	printf("Problem computing the size root cell of the tree %g %g %g\n",
	       xsize,ysize,zsize); fflush(stdout);
	exit(0);
      }
    
    BOXSIZE = BOXSIZE + 0.01*BOXSIZE;
    
    pos[0] = BoxXCenter;
    pos[1] = BoxYCenter;
    pos[2] = BoxZCenter;
    
    
    /* INFORMATION */
    printf("\n"); //fflush(stdout);
    printf("xmin=%g xmax=%g\n",xmin,xmax);  //fflush(stdout);
    printf("ymin=%g ymax=%g\n",ymin,ymax); //fflush(stdout);
    printf("zmin=%g zmax=%g\n",zmin,zmax); //fflush(stdout);
    printf("Boxsize=%g  Nparts=%d\n",BOXSIZE,NPARTICLES); //fflush(stdout);
    printf("x=%g y=%g z=%g\n",pos[0],pos[1],pos[2]); fflush(stdout);
    
    TREE_BUILD_NPARTICLES = NPARTICLES;
    DEAD_NODE = -1;
    
    MAX_ID_PARTS = 0;
    
    if(TREE_BUILD_NPARTICLES > domain_info[task].Nparts_per_node)
      MAX_ID_PARTS = TREE_BUILD_NPARTICLES;
    else
      MAX_ID_PARTS = domain_info[task].Nparts_per_node;
    
    GLOBAL_COUNTER = MAX_ID_PARTS;
    
    alloca_node(&root,NULL,BOXSIZE,pos,0);
    
    printf("Built %d tree nodes\n",GLOBAL_COUNTER-domain_info[task].Nparts_per_node);
    
    
    /* If this two options are enabled the tree can not be used to do
       NGB search because the positions of the nodes are going to be
       changed*/
    //GalCos_tree_mass(&root);
    //GalCos_tree_CM(&root);
    
    GLOBAL_COUNTER=0;
    
    return root;
    
}
