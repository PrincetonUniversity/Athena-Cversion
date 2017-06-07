---
title: Introduction to Particles
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Introduction to Particles

The particle code integrates a number of Lagrangian particles in the gas that are subject to the aerodynamic drag. The drag force to the particles is given by

<!---
![alt aerodynamic drag]({{site.baseurl}}/images/drag.png)
-->
\\[
F_{drag}=m{ {v_g-v_p}\over{t_{stop} } }
\\]

where \\(m\\) is the particle mass, \\(v_g\\) and \\(v_p\\) are the velocity of gas and particles, \\(t_{stop}\\) is the stopping time, characterizing the coupling between gas and particles. Particles can either respond to the gas motion passively, or also provide back reaction to the gas as momentum and energy feedback. Two particle integrators are designed to deal with particles with any \\(t_{stop}\\) from zero to infinity. A full description of the algorithms and test problems can be found in 'Particle-gas Dynamics with Athena: Method and Convergence' ([Bai & Stone, 2010](http://adsabs.harvard.edu/abs/2010ApJS..190..297B)). The particle integrators and the overall hybrid particle-gas scheme are both second order accurate in space and time.

Main applications of the particle-gas hybrid code include:

1. Aerodynamics of dust grains / solids in protoplanetary disks, with relevance on dust transport, heating, planetesimal/chondrule formation. Momentum feedback can be turned on or off depending on the application.

2. Tracking the motion of fluid parcels in any MHD simulations. In this case, simply use passive particles with \\(t_{stop}=0\\).

3. Proper modifications to the particle integrator can make it suitable to integrate particles that are subject to other forces, such as the Lorentz force (see [Lehe, Parrish & Quataert, 2009](http://adsabs.harvard.edu/abs/2009ApJ...707..404L)).

The current implementation of the particles works with MPI, but it **does NOT work with static mesh refinement**. Moreover, it **only works with the CTU integrator**.
