---
title: The <Domain> Block(s)
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Domain Block(s)

Parameters in these blocks control the properties of the *Domains*, and with MPI how each Domain is decomposed into *Grids*.  Without SMR, there
should be only one Domain block, and it must have level=0.  With SMR, there should be one block for
each Domain at each refinement level, and there must be one level=0 Domain (see the section on SMR in the [User Guide]({{site.baseurl}}/AthenaDocsUG) for details of 
Domains and Grids within the Mesh). 

With MPI, each Domain can be decomposed in any coordinate
direction, and into an arbitrary number of Grids.  This allows slab, pencil,
and block decompositions.  Different Domains can have different decompositions,
and the same processors can be used to update Grids on different Domains.
Athena can also automatically figure out the optimal domain
decomposition that minimizes the amount of data communicated; this option can be
chosen with `AutoWithNProc` parameter in each Domain block.
Again, see the section on SMR in the [User Guide]({{site.baseurl}}/AthenaDocsUG) for details of parallelization of 
Domains with SMR.

Example:

	<domain1>
	level           = 0         # refinement level this Domain (root=0)
	Nx1             = 8         # Number of zones in X-direction
	x1min           = -0.5      # minimum value of X
	x1max           = 0.5       # maximum value of X
	bc_ix1          = 4         # boundary condition flag for inner-I (X1)
	bc_ox1          = 4         # boundary condition flag for outer-I (X1)
	NGrid_x1        = 1         # with MPI, number of Grids in X1 coordinate
	AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
	
	Nx2             = 64        # Number of zones in Y-direction
	x2min           = -4.0      # minimum value of Y
	x2max           = 4.0       # maximum value of Y
	bc_ix2          = 4         # boundary condition flag for inner-J (X2)
	bc_ox2          = 4         # boundary condition flag for outer-J (X2)
	NGrid_x2        = 1         # with MPI, number of Grids in X2 coordinate
	AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
	
	Nx3             = 96        # Number of zones in X3-direction
	x3min           = -6.0      # minimum value of X3
	x3max           = 6.0       # maximum value of X3
	bc_ix3          = 4         # boundary condition flag for inner-K (X3)
	bc_ox3          = 4         # boundary condition flag for outer-K (X3)
	NGrid_x3        = 1         # with MPI, number of Grids in X3 coordinate
	AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
	
	<domain2>
	level           = 1         # refinement level this Domain (root=0)
	Nx1             = 16        # Number of zones in X1-direction
	Nx2             = 128       # Number of zones in X2-direction
	Nx3             = 96        # Number of zones in X3-direction
	iDisp           = 0         # i-displacement measured in cells of this level
	jDisp           = 0         # j-displacement measured in cells of this level
	kDisp           = 48        # k-displacement measured in cells of this level
	AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
	NGrid_x1        = 1         # with MPI, number of Grids in X1 coordinate
	NGrid_x2        = 1         # with MPI, number of Grids in X2 coordinate
	NGrid_x3        = 1         # with MPI, number of Grids in X3 coordinate
	`

**level:** Refinement level of this Domain.  There must be one and only one root 
(level=0) Domain block.  It does not need to be given first.

**Nx1, Nx2, Nx3:** number of grid cells in the x1-, x2-, and x3-directions in each
Domain.

**x1min, x2min, x3min:** x1-, x2-, x3-coordinate of left-edge of first cell for the
root (level=0) Domain.  Only specified for root level.

**x1max, x2max, x3max:** x1-, x2-, x3-coordinate of right-edge of last cell for the 
root (level=0) Domain.
The root level domain in the x1-direction spans `x1max-x1min`, the grid
spacing is `dx1=(x1max-x1min)/Nx1`, and the center of the
first cell is located at `x1=x1min + dx1/2`. Also, `x1max > x1min`
is required.  Similarly for the
x2- and x3-directions.  Only specified for root level.

**bc_ix1, bc_ox1:** integer flags for boundary conditions applied
at *inner* (left) and *outer* (right) edges of grid.  Currently three
values are implemented: 1 = reflecting, 2 = outflow (projection),
and 4 = periodic.  See [Boundary Conditions]({{site.baseurl}}/AthenaDocsUGBCs) for more information.
Similarly for the x2- and x3-directions.

**iDisp, jDisp, kDisp:**  Offset of origin of Domain from origin of root level.
Given in units of cells at this level.  Only specified for levels>0.


**NGrid_x1, NGrid_x2, NGrid_x3:** Number of MPI blocks in the x1-, x2-, x3-directions in this Domain.


**AutoWithNProc:** Set to number of processers desired for this Domain for automatic domain decomposition that optimizes efficiency.  In this case, the `NGrid_x1`, etc.
values are ignored.
