---
title: Mesh
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Mesh

In order to enable calculations with *Static Mesh Refinement (SMR)*, and parallel calculations using multiple processors with MPI,
the computational volume in Athena
is organized into a hierarchy of *Grids* and *Domains*.  The complete hierarchy, including all the Domains and
Grids in the calculation, is called the **Mesh**.

**Domains** are contiguous regions of the computational volume at a particular level of refinement.  There may be
more than one Domain at the same level, covering different parts of the overall volume.  Every calculation must have at
least one *root Domain*.

**Grids** are the smallest regions of a Domain updated on a single processor.  They are different from Domains only
in MPI parallel calculations.  Grids contain an arbitrary number of **cells**, discrete volumes on which the
conserved variables are averaged, stored, and updated.

An important restriction is that **Athena only works with a regular, uniform, logically rectangular array of cells** on each Grid.

In summary, the *Mesh* is divided into *Domains*, which are divided into *Grids*, which contain the *cells* on which
the solution is discretized.  

Dividing the Mesh into Domains
==============================

The figure below shows an example of how the computational volume can be divided into *Domains*:

![alt CompMesh]({{site.baseurl}}/images/CompMesh.png)

The entire computational volume is always covered by the lowest level *root Domain*.  In this example, two other regions in the volume are covered by
two separate level=1 Domains, which have twice the resolution of the root level.  Finally, one of the two level=1 Domains is
further refined into a single level=2 Domain.

Note the following.
 * There is always one, and only one, root Domain.  Since Athena is written in C, counting starts at zero, so the root Domain has level=0.
 * If SMR is *not* used, there is only one Domain, the root.
 * With SMR, there can be an arbitrary number of Domains at any level greater than zero.
 * With SMR, each Domain at every level can be further refined with an arbitrary number of Domains at the next highest level.
See the [User Guide to SMR]({{site.baseurl}}/AthenaDocsSMR) for more details on the Mesh hierarchy with SMR.

Dividing Domains into Grids
===========================

The figure below demonstrates the relationship between *Domains* and *Grids* in a *Mesh*.

![alt Grids]({{site.baseurl}}/images/Grids.png)

The region outlined in bold black denotes the computational volume, spanned by the single root Domain.  In an MPI parallel calculation,
this Domain has been further sub-divided into *Grids*.
Each Grid contains a patch of cells at the same resolution updated by a single processor.  In this example, the root Domain
has been divided into 3x3=9 Grids.

The bold red region denotes a level=1 Domain.  It has also been sub-divided into 3x3=9 Grids.  Note the size and
dimensions of the Grids at the two levels is arbitrary and need not be the same.  Domains can be divided into Grids along
each dimension independently, allowing for slab, pencil, or block decompositions.  See 
[Athena with Parallel Processors]({{site.baseurl}}/AthenaDocsUGParallel) in the [User Guide]({{site.baseurl}}/AthenaDocsUG) for details
of specifying different decompositions with MPI.

Note the following.
 * For serial calculations without SMR, there will be only one Grid, and it will span the entire root level Domain.
 * For parallel calculations without SMR, the number of Grids must equal the number of processors.  
 * Each Domain can be decomposed into Grids independently of every other Domain.  For example
 block decomposition can be used in some Domains, and slab in others.
 * Grids at different levels can overlap in an arbitrary (and potentially complicated) manner.
See the [User Guide to SMR]({{site.baseurl}}/AthenaDocsSMR) for more details on using mesh refinement with MPI.


 
