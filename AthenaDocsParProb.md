---
title: Problem Generators with Particles
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Problem Generator with Particles

To write a problem generator, one must take the following steps.

1. Set the total number of particles (`nparticle`), and reallocate the particle array if necessary. For the grid pointer pGrid, one example is given below


        if (pGrid->nparticle+2 > pGrid->arrsize)
            particle_realloc(pGrid, pGrid->nparticle+2);


    where `nparticle`+2 is the new array size. Note that the array size `arrsize` must be at least one more than `nparticle`.

2. Set the parameters needed to calculate the stopping time. When necessary, also set the particle size (`rad` defined in `Grain_Property`, see [Data Structure]({{site.baseurl}}/AthenaDocsParDataStruct)) for each particle type.

    For `tsmode`=3 (fixed stopping time), one simply sets the array `tstop0[]` for the stopping time of each particle type.

    For `tsmode`=2 (Epstein drag law), one needs to set the global array `grrhoa[]` for each particle type, so that the stopping time is calculated by `grrhoa/(d*cs)`, where `d` and `cs` are the gas density and sound speed.

    For `tsmode`=1 (general drag law), in addition to setting the array `grrhoa[]`, one further needs to set the global variable `alamcoeff`, so that the ratio between particle size and gas mean free path is given by `alamcoeff*rad*d`.

3. If feedback is included, set the mass for each type of particles (in `grproperty`, see [Data Structure]({{site.baseurl}}/AthenaDocsParDataStruct)). Note that particles in the code are "super-particles", each representing a swarm of real particles. The particle mass is the inertia of each super-particle, and has the same unit of gas density.

4. Specify the initial positions and velocities of all the particles, and assign a unique id for each particle.

5. The user has the freedom of adding user-defined forces to the particles in the following function


        void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                               const Real v1, const Real v2, const Real v3);


    The force vector `ft` is to be given as a function of particle position (e.g., in a static gravitational potential).  **Be warned:**
    if the force vector is a function of velocity in general this will require modifications to the integrator.

6. The user can also study the motion of particles in a fake gas velocity field. This can be achieved by letting the gas be in hydrostatic equilibrium (all the time), while creating a fake velocity field using the following function in the problem generator


        void gasvshift(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3);


    Here the fake gas velocities `u1`, `u2` and `u3` are to be set as a function of position `x1`, `x2`, `x3`.
