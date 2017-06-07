---
title: Introduction to SMR
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Introduction

The following are the most important features of the SMR algorithm implemented in Athena.

Static, not Adaptive, Mesh Refinement
=====================================

The mesh refinement does not adapt to the solution, instead the number and location of mesh levels is defined by the
user at the start through parameters in the input file.  This is ideal for problems in which regions requiring
higher resolution are known *a priori* and fixed (for example, the midplane of an accretion disk, the center of
a cluster of galaxies, or the beam of a supersonic jet), because it allows for much better load balancing and much better
efficiency on parallel computers than AMR. 

Block Adaptive
==============

The refined region must contain blocks of at least 4 cells in each dimension (although blocks this small should never be
used in practice).  There is no upper limit on the size of refined regions (other than the overall size of the root Domain).
However, the maximum refinement ratio between levels is limited to two.

Uniform Time step
=================

All levels in the mesh hierarchy are integrated with the same time step.  Most AMR codes integrate the fine resolution
regions with a proportionally smaller time step, however this has two problems:

 * In MHD, the smallest time step in the calculation might be set by regions with the coarsest resolution, if the density in
   these regions is small (producing a large Alfven speed).  For example, in stratified accretion disks, the highest Alfven
   speed is in the upper regions (corona) of the disk, where the resolution is the lowest.
 * Load balancing on large number of cores is problematic with variable time steps.  The fine grid solution requires the coarse grid
   to be integrated first to provide boundary conditions.  If every processor does not contain some fine grid cells, then those processors
   will sit idle waiting for data.  It may be simpler and more efficient just to use those processors to integrate every cell at the same time step.
 * Enforcing the divergence-free constraint on the magnetic field is much more complicated with variable time steps.

Arbitrary Parallelization
=========================

The mesh hierarchy can be spread across multiple processors in any arbitrary fashion.  For example, suppose there are *N* levels to be
updated by *N* processors.  Each level could be assigned to one processor, or each level could be decomposed into *N* Grids which are
each updated by one processor, or some mixture of the two is allowed.  The only restriction is that the total number of Grids must be an
integer multiple of *N*.
