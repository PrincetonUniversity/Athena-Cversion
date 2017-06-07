---
title: Configuring & Compiling with Special Relativity
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Configuring

Configuring Athena to integrate the equations of relativistic (magneto-)hydrodynamics is accomplished by means of the configure script described in the User Guide. The simplest set of options necessary to configure the code in this fashion are:


        ./configure --enable-special-relativity --with-order=2p --with-integrator=vl --with-flux=hlle --with-gas=hydro


These choices will configure the code to integrate the equations of relativistic hydrodynamics using the Van Leer integrator, second order primitive variable reconstruction and HLLE fluxes. Note that the integrator must  be specified as Van Leer and reconstruction with limiting in the primitive variables must be selected.
These choices can be modified as follows:

### Optional physics (package=*choice*) controlled by configure

| Package | Choice | Comments|
|---------|--------|---------|
| Flux | exact | exact nonlinear Riemann solver (hydro only)|
| | hlle | HLLE Riemann solver|
| | hllc | HLLC Riemann solver|
| | hlld | HLLD Riemman solver|
| Order| 2p | second-order reconstruction with limiting in the primitive variables|
| | 3p | third-order reconstruction with limiting in the primitive variables|

In addition, the algorithmic stability can be greatly enhanced by activating the first-order flux-correction feature via `--enable-fofc`. This activates the various fallback procedures to correct unphysical conservative states described in 
[Beckwith & Stone, 2011](http://adsabs.harvard.edu/abs/2011ApJS..193....6B). Static mesh refinement and domain decomposition on parallel systems can be accessed in the same way as the Newtonian code by specifying `--enable-smr` and `--enable-mpi` respectively. Note that fargo, h-correction and hllalllwave features are not available for the relativistic algorithms.

Our recommended choices to integrate the equations of relativistic hydrodynamics are as follows:


        ./configure --enable-special-relativity --with-order=2p --with-integrator=vl --with-flux=hllc --enable-fofc --with-gas=hydro


Similarly, our recommended choices to integrate the equations of relativistic magnetohydrodynamics are as follows:


        ./configure --enable-special-relativity --with-order=2p --with-integrator=vl --with-flux=hlld --enable-fofc


Once code configuration is established, compilation is performed in exactly the same fashion as described in [Compiling]({{site.baseurl}}/AthenaDocsUGCompile)
