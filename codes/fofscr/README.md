

						FOF v1.1 

			A Group Finder for N-body Simulations

					October 26, 1994

Changes from v1.0:
	o 	Fixed bug in tree building, this bug only affected cases where
		a very small "bucket" size was chosen and the number of particles
		was not a power of two.

Included are:
	README
	Makefile
	kd.c
	kd.h
	main.c
	tipsydefs.h

To build:

	> make

To get further information contact:

	hpccsoft@astro.washington.edu

History:

This software has been modified to read particle data in ascii format (February 2013)
