#ifndef INIT_GALCOS_H
#define INIT_GALCOS_H

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>

#define EMPTY_FLAG -1
#define SUCCES 0
#define FAILURE 1
#define BH_OPENING 0.7

/* Some physical and astrophysical constants (in cgs units) */

#define G_GRAVITY         6.672e-8
#define HUBBLE            3.2407789e-18	/* in h/sec */
#define HUBBLE_TIME       3.09e+17   // in sec/h

extern int Npart_Total,MINIMUM_MEMBERS,NGB_MAX,FOF_PROCEDURE;
extern int FLAG_SUBFIND,FLAG_WRITE_INDIVIDUAL_CLUSTERS;
extern int MINIMUM_NSUBSTRUCT,NCLUSTERS; 

extern float BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON;
extern float Llenght,COSMIC_TIME,b_Link;
extern float GRAV_SOFT,PARTMASS;

extern double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
extern double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
extern double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;


struct queueNode
{
  int data;
  struct queueNode *nextPtr;
};

struct ptnode
{
  int NumParts;             // Number of particles in this node
  short int tag;            // tag=1 if twig, tag=0 if leaf
  int index;
  int NSon[8];
  //struct tcola *partsInNode;
  struct queueNode *headPtr;
  struct queueNode *tailPtr;
  struct ptnode **sons;      // Pointers to the 8 possible sons of the node
  float size;               // Side size of the node
  float pos[3];
  float mass;  
};

typedef struct ptnode TREENODE;
typedef TREENODE *TREENODEPTR;

extern struct COLA
{
  int head;
  int tail;
  int cola_Nmembers;
  int *inputs;
}temp_cola;

extern struct COLA_float
{
  int head;
  int tail;
  int cola_Nmembers;
  float *inputs;
}NGB_DISTcand;

extern struct halo
{
  float mass;
  float pos[3];
  float vel[3];
  float Rvir,Mvir;
  int ID_CenterHalo;
  int Nmembers;
  int NDomain_particles;
  int IDcluster;
  int *Halo_particles;
  int *Domain_particles;
  int *Domain_particles_index;
  int NDomain_particles_absolute;
  int Nvir;
}*Halos;

extern struct part
{
  float pos[3];
  float vel[3];
  float mass;
  float EP;
  int id;
  int Oid;
  int Cluster_ID;
}*Particle;


void GalCos_InitHalo(int size_init);
float distance(float xi, float yi, float zi, float xj, float yj, float zj);
void GalCos_BuildCluster(struct COLA cola,int Ncluster,int HaloIDCenter);
void GalCos_Halo_prop(struct halo *Haloinp);

int gsl_fisort(int dimension,float *fvector,int *ivector);

#endif /* INIT_GALCOS_H */
