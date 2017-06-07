---
title: Memory and Arrays
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Memory and Arrays

Athena uses dynamic memory allocation, meaning that all arrays are
created at run time using `malloc` or `calloc`.  This has the 
great advantage that problems
of different dimensions (1D, 2D, or 3D), and problems of different sizes,
can be run without the need for recompiling.

The indexing convention used in Athena requires that all multidimensional
arrays must have the x1 index as the fastest incrementing index.
In this way, data in cell `i-1` should be contiguous with data in cell
`i`, whereas data in cell `j-1` will be `Nx1+2*nghost` data elements
away from the data in cell `j`, and the data in cell `k-1` will be
`(Nx1+2*nghost)*(Nx2+2*nghost)` data elements away from the data in
cell `k`.  In order to achieve this in C, arrays must be referenced as
`A[k][j][i]`.  The file `ath_array.c` contains functions for
creating (allocating) and destroying (deallocating) 2D and 3D arrays.  For example,
to allocate a 3D array of `ConsS` structures with dimensions Nx3 x Nx2 x Nx1, use

	if((MyArray = (ConsS***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(ConsS)) == NULL) 
	    ath_error("Failed to allocate MyArray\n");

If the calloc fails, `ath_error` will print a warning to stdout and terminate execution.

On each Grid, the first active cell is
labeled `is` and the last is labeled `ie` (similarly for `js` and `je`,
and `ks` and `ke`).  Thus, to stride across all active zones requires
triply nested `for` loops:

	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
	        ....
	      }
	    }
	  }

Note the inner loop should always be over the `i`-index.

Arrays of the dependent variables are always allocated as 3D (i.e., as having 3 indices), however the `k` (in 1D and 2D) and/or
`j` (in 1D) indices may have only a single element.  The code automatically sets `ks=ke=0` in 1D 
and 2D, `js=je=0` in 1D, thus arrays in 1D calculations can either be referenced in a single loop

	  for (i=is; i<=ie; i++) {
	    MyArray[0][0][i] = ...
	  }

or in a single loop with `j=js` or `je`, and `k=ks` or `ke`

	  for (i=is; i<=ie; i++) {
	    MyArray[ks][js][i] = ...
	  }

or in a triply-nested loop

	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
	        MyArray[k][j][i] = ...
	      }
	    }
	  }

In the last case, the `j` and `k` indices will only have the single value (zero). 
Thus, such loops will work (will not address arrays out-of-bounds) for a problem of any
dimensions.
