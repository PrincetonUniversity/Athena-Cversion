---
title: Introduction to Special Relativity in Athena
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Introduction


Much of the implementation of Newtonian MHD within Athena carries directly over to the relativistic case. Due to the increased complexity of the RMHD equations compared to their Newtonian counterparts, the algorithmic options are more limited for RMHD and possess additional complexity. A full description of the algorithms for RMHD implemented within Athena can be found in 'A Second Order Godunov Scheme for RMHD' ([Beckwith & Stone, 2011](http://adsabs.harvard.edu/abs/2011ApJS..193....6B)). Here we highlight the most important elements

### Integrator & Reconstruction Options

The characteristic decomposition of the RMHD equations necessary for the CTU integrator has not yet been implemented with Athena. Integration of the RMHD equations therefore relies on a variant of the MUSCL-Hancock type Van Leer integrator.  The CTU integrator is not available. For the same reason, reconstruction of variables from cell- to face-centers must be performed using the 'primitive' rather than 'characteristic' variant of the reconstruction algorithms. See Configuring & Compiling for further information.

### Primitive Variable Calculation

The relationship between conserved and primitive variables in RMHD is non-linear and cannot be written in a closed form. Athena implements a non-linear root finding algorithm for performing this inversion along with several fallback strategies (described below) to ensure that a physical primitive state is obtained from a given conserved state. Due to this non-linear relationship, the primitive variables should be output from the code for later analysis (see The Input File for further information).
The simulation must be initialized with a conserved state that corresponds to a physical primitive state. The easiest way to ensure this is the case is to initialize problems in the primitive variables and then convert these variables to a corresponding set of conserved variables using a supplied routine within Athena. Further details can be found in User Defined Problems.

### Riemann Solvers

The available solvers for RMHD are the HLLE/HLLC/HLLD Riemann solvers, implemented as described in Beckwith & Stone (2011). The same Riemann solvers with the addition of an exact Riemann solver are available for RHD. The Roe solver, FORCE flux and two-shock Riemann solvers are not available for either RHD or RMHD.

### First Order Flux Correction

The first order flux correction is available for both RHD and RMHD. The user is encouraged to use the code with this option enabled as it provides an extremely effective strategy for fixing unphysical cells (a common problem in both RHD and RMHD). When this option is activated in RMHD, additional fallbacks (the entropy fix and velocity fix) are also enabled. The code will print out diagnostics containing information regarding the number of times that these fixes are used each time-step (if this number is zero, then no information is printed).
