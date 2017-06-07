---
title: Input Files with Particles
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Input File with Particles


For the particle code to be properly initialized (handled by the code), one needs to provide information in the block <particle> in the input file. The items of this block include


        <particle>
        
        partypes        = 1         # number of types of particles
        parnumgrid      = 8192      # number of particles per grid per type
        parnumcell      = 1         # number of particles per cell per type
        
        integrator      = 2         # particle integrator (1: explicit; 2: semi-implicit; 3: fully-implicit) 
        interp          = 2         # interpolation scheme (1: CIC; 2: TSC; 3: polynomial)
        tsmode          = 3         # stopping time calculation mode (1: General; 2: Epstein; 3: fixed);   
        
        tshuf           = 20        # time interval to shuffle the particles


**partypes**: the number of particle types, with the default being 1. It is passed to global variable `npartypes` and can NOT be modified in the [problem generator]({{site.baseurl}}/AthenaDocsParProb).

**parnumgrid** and **parnumcell**: in the particle initialization, the total number of particles (hence the size of the particle array in the grid) can either be specified from the number of particles per grid per type (`parnumgrid`), or from the number of particles per cell per type (`parnumcell`). If neither of them are provided, the default particle number is 1000. The number of particles can further be modified in the [problem generator]({{site.baseurl}}/AthenaDocsParProb).

**integrator**: three 2nd order accurate particle integrators can be chosen (explicit, semi-implicit and fully implicit). The semi-implicit (label=2) is the recommended (and default) integrator which gives better accuracy. However, the fully implicit integrator (label=3) should be used in the very stiff regime when the stopping time is less than the time step, in particular, when the particles are used as gas tracers. The particle integrator can also be assigned in the [problem generator]({{site.baseurl}}/AthenaDocsParProb), where different types of particles can be assigned with different particle integrators.

**interp**: three particle-gas interpolation schemes can be chosen, namely, CIC (cloud in a cell, label=1),  TSC (triangular shaped cloud, label=2) and QP (quadratic polynomial, label=3). TSC is recommended (and default) since it is of higher order accuracy and is positive definite.

**tsmode**: the method for calculating the stopping time. For simple applications, set `tsmode`=3 for fixed stopping time to be specified in the problem generator (tstop0). More realistically, set `tsmode`=2 for Epstein drag or tsmode=1 for the most general drag law, in which case some other conversion parameters need to be specified in the [problem generator]({{site.baseurl}}/AthenaDocsParProb).

**tshuf**: the time interval to sort the particles to improve the cache efficiency. By default it is set to 0 (do not shuffle).
