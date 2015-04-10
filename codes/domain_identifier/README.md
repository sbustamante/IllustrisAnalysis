README FOR DOMAIN_IDENTIFIER (June 09 of 2010)


Domain identifier is a program runing on MPI to compute the domain of
dark matter FOF halos identified in simulations. The domains are
defined as the set of particles that are nearest to a halo than to
other halos, when using distances weighted by the virial radius of the
halos. That is, a particle located at distance d(i) belongs to a halo
of virial radius Rvir(i) if d(i)/Rvir(i) is the minimum quantity when
computed across the overal set of halos.

The results are stored in a binary table, in wich the label of the
halo for every particle is stored.

The results of this program are needed for Densfield, who reads the
table and computes the Delaunay volumes of the particle distribution
and then build the density field of a populations of halos.

this version of the code make the constrain on the halo populations
considering only those halos with masess above Mth but using the
virial mass, not the fof mass.  This is especially important for the
study of the mass distribution of halos and environments in the study
of the evolution of the density profile.

JUAN CARLOS MUNOZ CUARTAS
ASTROPHYSIKALISCHES INSTITUT POTSDAM