---
title: Geometry of Levels
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Geometry

The Mesh hierarchy in Athena is very flexible.  For example, it allows for an arbitrary number of levels, and
an arbitrary number of independent Domains at each level. However, there are some important restrictions on the geometry of levels, 
as described below.  

The code tests whether the Mesh hierarchy violates any of these restrictions at the start of execution of each run; 
it should quit with a warning message if a problem is found.

Restrictions on Domains touching the boundary of the root Domain
================================================================

Fine (level>0) Domains may touch the boundaries of the root (level=0) Domain (see the examples below).
In such cases, the physical boundary conditions
applied to the root Domain will also be applied to the finer level Domains.

<!---
![alt SMRDoms1]({{site.baseurl}}/images/SMRDoms1.png)
![alt SMRDoms2]({{site.baseurl}}/images/SMRDoms2.png)
![alt SMRDoms3]({{site.baseurl}}/images/SMRDoms3.png)
-->
<img width="32%" src="{{site.baseurl}}/images/SMRDoms1.png"/>
<img width="32%" src="{{site.baseurl}}/images/SMRDoms2.png"/>
<img width="32%" src="{{site.baseurl}}/images/SMRDoms3.png"/>

However, with periodic BCs, the finer level Domains must touch *both* sides of the root Domain, or none at all.  Thus, the first
two examples above are not allowed with periodic BCs, while the third is allowed.

Space is required between Domains at different levels
=====================================================

The edges of Domains at different levels cannot touch, unless they also touch the boundary of the root Domain.  In addition, there must be at
least *nghost* (generally 4) cells on the finest level between the boundaries of Domains at different levels.  This space is needed so that the ghost 
zones on the finer level can be set by a prolongation of *active cells* (not ghost zones) on the coarser level.  The first example below is
allowed (except with periodic BCs), and the following two are not allowed.

<!---
![alt SMRDoms4]({{site.baseurl}}/images/SMRDoms4.png)
![alt SMRDoms5]({{site.baseurl}}/images/SMRDoms5.png)
![alt SMRDoms6]({{site.baseurl}}/images/SMRDoms6.png)
-->
<img width="32%" src="{{site.baseurl}}/images/SMRDoms4.png"/>
<img width="32%" src="{{site.baseurl}}/images/SMRDoms5.png"/>
<img width="32%" src="{{site.baseurl}}/images/SMRDoms6.png"/>

Domains at the same level cannot touch
======================================

In order to provide space for the ghost zones of each Domain to be set by prolongation, Domains at the same level cannot touch, or
be closer than *nghost* (usually 4) cells on that Domain.  If Domains touch, then in principle the ghost zones can be set by mapping them to the
corresponding active cells on the neighboring Domain, however this logic is complicated with face-centered magnetic fields in MHD
and has not been implemented (yet).  None of the examples below are allowed.

<!---
![alt SMRDoms7]({{site.baseurl}}/images/SMRDoms7.png)
![alt SMRDoms8]({{site.baseurl}}/images/SMRDoms8.png)
-->
<img width="32%" src="{{site.baseurl}}/images/SMRDoms7.png"/>
<img width="32%" src="{{site.baseurl}}/images/SMRDoms8.png"/>

Maximum refinement ratio between successive levels is two
=========================================================

The ratio of grid spacing between successive levels is always restricted to be two.  Any arbitrary number
of levels is allowed to reach any desired resolution, but each level in the Mesh hierarchy must involve a
decrease in the size of the grids cells of two and only two.
