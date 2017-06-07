---
title: Domains
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Domains

The Root Domain
===============

Every calculation must contain one *root Domain*.  This Domain has indices $(nl,nd)=(0,0)$ in the 
[Mesh structure]({{site.baseurl}}/AthenaDocsPGMeshS), e.g. it can be referenced as `Mesh.Domain[0][0]`.

The figure below shows the relationship between the physical coordinates (x1,x2) and the coordinates of the root Domain.
Although this figure (and all the figures in this section) has been drawn for a 2D calculation, it is straightforward
to generalize (restrict) it to 3D (1D). 

![alt RootDomainX]({{site.baseurl}}/images/RootDomainX.png)

The region *inside* the bold black line is the physical volume of the calculation.  It spans the region 
`RootMinX[i]` ≤ `x[i]` ≤ `RootMaxX[i]`, for `i=0,2`.
To be more precise, in each direction `i`
 * `RootMinX[i]` is the *left-interface* of the first actively updated cell
 * `RootMaxX[i]` is the *right-interface* of the last actively updated cell.
Again, to emphasize, *these values are the locations of cell interfaces, not the locations of cell-centers*.

It is important that the location of the origin of x1-x2 coordinate system is arbitrary with respect to the region spanned by
the root Domain.  This allows more freedom when including physical effects that depend on the coordinates, for example a
gravitational potential due to a point mass located at the origin.

Other Domains
=============

For SMR calculations, there will be more than just the root Domain.  The figure below shows how the location of these other
Domains are specified in the x1-x2 coordinate system.

![alt DomainsX]({{site.baseurl}}/images/DomainsX.png)

The bold black lines outline the region spanned by the root Domain.  The bold red line outlines the region spanned by
a Domain at the next level.  It spans the region 
`MinX[i]` ≤ `x[i]` ≤ `MaxX[i]`, for `i=0,2`.  As in the case of the root Domain, in each direction `i`
 * `MinX[i]` is the *left-interface* of the first actively updated cell on each Domain.
 * `MaxX[i]` is the *right-interface* of the last actively updated cell on each Domain.
