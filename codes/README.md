Toolkit
==============================================
**Sebastian Bustamante**
*(Universidad de Antioquia)*

Description
-----------------------------------------------------------------------------------------
This is a set of miscellaneous functions for Gadget snapshot datafiles. The next functions
are available:

- **Cutter.out**: this functions selects particles within some given slide.
- **Ascii.out**: this functions stores the particle data with some sampling rate in an ascii format file.
- **Halo.out**: this functions calculates a FOF index of the particles using the washington FOF finder.
- **Evolution.py**: this functions makes a video of the evolution of a slide of the simulation.
- **PhaseDiagram.py**: this functions makes a video of the evolution of the phase diagram of gas particles.
    
    
Requirements
-----------------------------------------------------------------------------------------
For a proper running of these codes it is required local copies of the [Washington FOF halo finder](https://github.com/forero/HackFOF/tree/master/src)
and [Domain Identifier](https://github.com/sbustamante/IllustrisAnalysis/codes/domain_identifier).
This packages may require another requisites, consult their Readmes.