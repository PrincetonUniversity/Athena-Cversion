---
title: MPI with SMR
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/MPI with SMR

The SMR algorithm in Athena allows for very flexible decompositions using MPI.  Each Domain at every level can
be decomposed into Grids (each updated by a single processor) independently.  The code automatically figures out the
overlap between Grids on different processors, and handles the communication necessary for the restriction and prolongation 
operations.  

*Since different decompositions require differing amounts of data to be communicated, users are strongly recommended
to experiment in order to find the most efficient layout for their specific problem.*

The decomposition of Domains into Grids with MPI is specified using parameters in the `<domain>` blocks in the
input file (see [Domain Blocks]({{site.baseurl}}/AthenaDocsUGDomainBlck) in the User Guide).  For example, to specify two levels
in a 3D calculation, with 8 Grids per level in a 2x2x2 configuration in the root level, and a 4x2x1 configuration in
the level=1 Domain, use

        <domain1>
        level           = 0         # refinement level this Domain (root=0)
        Nx1             = 128       # Number of zones in X1-direction
        x1min           = 0.0       # minimum value of X1
        x1max           = 3.0       # maximum value of X1
        bc_ix1          = 4         # boundary condition flag for inner-I (X1)
        bc_ox1          = 4         # boundary condition flag for outer-I (X1)
        NGrid_x1        = 2         # with MPI, number of Grids in X1 coordinate
        AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
        
        Nx2             = 64        # Number of zones in X2-direction
        x2min           = 0.0       # minimum value of X2
        x2max           = 1.5       # maximum value of X2
        bc_ix2          = 4         # boundary condition flag for inner-J (X2)
        bc_ox2          = 4         # boundary condition flag for outer-J (X2)
        NGrid_x2        = 2         # with MPI, number of Grids in X2 coordinate
        
        Nx3             = 64        # Number of zones in X3-direction
        x3min           = 0.0       # minimum value of X3
        x3max           = 1.5       # maximum value of X3
        bc_ix3          = 4         # boundary condition flag for inner-K (X3)
        bc_ox3          = 4         # boundary condition flag for outer-K (X3)
        NGrid_x3        = 2         # with MPI, number of Grids in X3 coordinate
        
        <domain2>
        level           = 1         # refinement level this Domain (root=0)
        Nx1             = 128       # Number of zones in X1-direction
        Nx2             = 64        # Number of zones in X2-direction
        Nx3             = 64        # Number of zones in X3-direction
        iDisp           = 32        # i-displacement measured in cells of this level
        jDisp           = 16        # j-displacement measured in cells of this level
        kDisp           = 16        # k-displacement measured in cells of this level
        AutoWithNProc   = 0         # set to Nproc for auto domain decomposition
        NGrid_x1        = 4         # with MPI, number of Grids in X1 coordinate
        NGrid_x2        = 2         # with MPI, number of Grids in X2 coordinate
        NGrid_x3        = 1         # with MPI, number of Grids in X3 coordinate

This will create Grids of size 64x32x32 on the root level, and 32x32x64 on level=1.  

This configuration could be run on either 8 or 16 processors.  In the former case, each processor would run one Grid from each level.
In the latter case, each Grid would be run on a separate processor.
