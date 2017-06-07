---
title: Mesh Structure
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Mesh Structure

The overall hierarchy of Domains and Grids is stored in the `MeshS` structure defined in `/athena/src/athena.h`:

        typedef struct Mesh_s{
          Real RootMinX[3];   /* min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
          Real RootMaxX[3];   /* max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
          Real dx[3];     /* cell size on root Domain [0,1,2]=[x1,x2,x3] */
          Real time, dt;  /* current time and timestep for entire Mesh */
          int Nx[3];      /* # of zones in each dir on root Domain [0,1,2]=[x1,x2,x3] */
          int nstep;                 /* number of integration steps taken */
          int BCFlag_ix1, BCFlag_ox1;  /* BC flag on root domain for inner/outer x1 */
          int BCFlag_ix2, BCFlag_ox2;  /* BC flag on root domain for inner/outer x2 */
          int BCFlag_ix3, BCFlag_ox3;  /* BC flag on root domain for inner/outer x3 */
          int NLevels;               /* overall number of refinement levels in mesh */
          int *DomainsPerLevel;      /* number of Domains per level (DPL) */
          DomainS **Domain;          /* array of Domains, indexed over levels and DPL */
          char *outfilename;           /* basename for output files containing -id#  */
        }MeshS;

Note that Domains are stored as a 2D array whose dimensions are the number of levels `nl`, and the number
of Domains per level `nd`.  This array is allocated (in `init_mesh.c`) in a way that allows `nd` to
be different for each `nl`.
