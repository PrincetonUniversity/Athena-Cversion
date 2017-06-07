---
title: Particle History Dump
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Particle History Dump

Similar to the [gas history dump]({{site.baseurl}}/AthenaDocsUGhst), the particle history file contains a formatted table of a variety of volume averaged value associated with the particles. To add the particle history output, one simply sets the output format to "`phst`" as follows:

	<output1>
	out_fmt = phst		# particle history dump
	dt      = 0.2		# time step between output


At constant time increment `dt`, one line is created containing the bulk information about all particles in the simulation, including time, maximum particle density (`d_max`), stiffness of the particle-gas coupling (stiffmax, defined in equation 11 of [the ApJS method paper](http://adsabs.harvard.edu/abs/2010ApJS..190..297B)), energy loss rate, total particle mass, momentum and kinetic energy. In addition, for each type of particles, one line is written recording the averaged and dispersion of positions, velocities of all particles with the same type. In the following example, there are seven types of particles, so eight lines are produced every time. The new data will be appended to the old file, and an empty line is written to separate the output at different times.


	#Global scalars:
	#  [1]=time      [2]=d_max     [3]=stiffmax  [4]=Edot      [5]=mass      [6]=x1 Mom.   [7]=x2 Mom.   [8]=x3 Mom.   [9]=x1-KE     [10]=x2-KE    [11]=x3-KE  
	#Particle type dependent scalars:
	#  [1]=x1 avg    [2]=x2 avg    [3]=x3 avg    [4]=v1 avg    [5]=v2 avg    [6]=v3 avg    [7]=x1 var    [8]=x2 var    [9]=x3 var    [10]=v1 var   [11]=v2 var   [12]=v3 var 
	#
	   1.20000e+03   4.86787e+02   1.56761e+00   3.23519e-04   2.50663e-01  -1.96212e-04  -2.02956e-07   8.58908e-03   2.93034e-06   4.12140e-07   1.64247e-04
	   1.33760e-03  -3.15545e-04   0.00000e+00   1.17057e-03  -3.15957e-05   2.25529e-02   5.42656e-02   1.63370e-02   0.00000e+00   3.17432e-03   1.97289e-03   1.07885e-02
	   5.05719e-03  -1.45659e-04   0.00000e+00   1.23330e-03  -8.32332e-05   2.58420e-02   5.44111e-02   1.32023e-02   0.00000e+00   3.39601e-03   2.03703e-03   9.13660e-03
	   1.03176e-02  -3.79540e-04   0.00000e+00   1.29432e-03  -1.28075e-04   2.91581e-02   5.71123e-02   1.03522e-02   0.00000e+00   3.55413e-03   2.06492e-03   7.99096e-03
	   9.07115e-03   1.73965e-04   0.00000e+00   5.84681e-04  -1.68525e-04   3.17749e-02   5.20556e-02   7.99246e-03   0.00000e+00   3.79864e-03   2.19918e-03   6.80819e-03
	   2.30548e-02   5.84561e-04   0.00000e+00  -1.63216e-03  -1.21791e-04   3.69335e-02   5.86155e-02   5.09936e-03   0.00000e+00   4.16020e-03   1.95824e-03   6.52555e-03
	   4.51122e-02  -4.75647e-04   0.00000e+00  -2.08487e-03   1.23535e-04   4.70470e-02   4.07878e-02   1.67342e-03   0.00000e+00   3.30800e-03   9.62542e-04   5.41350e-03
	   1.66136e-02   3.76637e-04   0.00000e+00  -6.04525e-03   4.04017e-04   4.65500e-02   5.02019e-02   1.46032e-03   0.00000e+00   6.18724e-03   9.32662e-04   3.29413e-03


The user can also enroll new particle history variables in a similar way to [the gas case]({{site.baseurl}}/AthenaDocsUGUserExpress) by adding a function call to `dump_parhistory_enroll()` anywhere in the problem generator, which is defined as


	void dump_parhistory_enroll(Parfun_t pfun, char *label, int set);


Here `Parfun_t` is a function type defined in `/athena/src/athena.h`, with two arguments of types `GridS*` and `GrainS*`, which are used to calculate the desired quantity from each particle. The string "`label`" is intended for the name of this quantity. For the last variable, `set=0` means averaging this variable over all particles (hence in the first line of the output), while `set=1` means averaging this variable for each particle type (hence in the remaining lines of the output).
