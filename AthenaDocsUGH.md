---
title: H-correction
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/H-correction

In many kinds of Godunov schemes (including Athena), strong planar shocks in multidimensions that propagate along a grid direction
can be subject to a numerical instability that causes large amplitude perturbations in the shock front at the grid scale.
This is often called the *carbuncle instability*.

Fixing the carbuncle instability requires adding extra dissipation *in the direction perpendicular to the shock only*.  We have
found the best way to do this is using the H-correction, see Appendix C in the [ApJS Method Paper](http://adsabs.harvard.edu/abs/2008ApJS..178..137S)
for details.

To use the H-correction, use the Roe flux, and enable the correction during configure using:

	% configure --with-flux=roe --enable-h-correction


*The H-correction is not needed for most problems.*  It is only required when there are strong planar shocks propagating parallel to
grid lines for a significant fraction of a dynamical time.
