---
title: Passive Scalars
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Passive Scalars

An arbitrary number of passive scalars which will be advected with the flow can be added to any simulation, using

	% configure --with-nscalars=#

where `#` is an integer specifying the number if scalars desired.

These scalars can be used to follow mixing in the flow.  In addition, the `Userwork_in_loop()` function in the
problem generator can be used to add source and sink terms to the advection equations for these scalars.  This allows
chemical or ionization/recombination models to be added (albeit only at first-order, using operator splitting).

The scalars are stored as an array `s[n]` in the [Conserved Variable Structure]({{site.baseurl}}/AthenaDocsPGConsPrim) which in turn is
stored in the [Grid Structure]({{site.baseurl}}/AthenaDocsPGGridS) in the file `/src/athena.h`

Outputs such as binary and vtk dumps will automatically include all passive scalars.  Images or outputs of specific members
of the scalars array can be added using [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress).

*Not all Riemann solvers have been extended to work with passive scalars.*  The Roe, HLLC, and HLLD solvers work with passive scalars.
