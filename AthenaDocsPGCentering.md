---
title: Centering of Variables
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Centering

A final, important characteristic of the Mesh used in Athena is the centering of the conserved variables in Grid cells.
The figure below shows this centering.  The volume-averaged conserved variables $U$ = (density, momentum, total energy)
are stored at cell-centers.  The components of the area-averaged magnetic field are stored at the corresponding cell faces.
Additional copies of the volume-averaged cell-centered magnetic field are stored in $U$, however these are not
the fundamental description of the field (they are computed from the face-centered fields).

![alt Centering]({{site.baseurl}}/images/Centering.png)

The cell-centered variables are stored in the [Conserved Variable structure]({{site.baseurl}}/AthenaDocsPGConsPrim), which is
a 3D array in the [Grid structure]({{site.baseurl}}/AthenaDocsPGGridS).  Components of this array are addressed using integer
indices with the $i$ (x1) index incremented fastest: $U_{k,j,i}$.

The face centered magnetic fields are stored as 3D arrays directly in the [Grid structure]({{site.baseurl}}/AthenaDocsPGGridS).
The array element on the *left-interface* is given the index $i$, that is $B_{x,i-1/2,j,k} = B1i[k][j][i]$.
