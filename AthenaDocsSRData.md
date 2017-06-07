---
title: Data Output with Special Relativity
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Data Output

All of the output options described in [Overview]({{site.baseurl}}/AthenaDocsUGOutputOverview) are available for relativistic (magneto-)hydrodynamics. We recommend that the user outputs the primitive variables along with any other quantities necessary, due to the non-linear relationship between conserved and primitive variables in relativistic (magneto-)hydrodynamics. All interesting physical quantities can be constructed directly from these variables. Note that the velocity contained in the primitive data dumps is the velocity three-vector. Note also  that the magnetic fields output in the primitive data dumps are the cell-centered fields calculated from the face-centered fields evolved by the induction equation. The magnetic fields used to calculate the Lorentz force can be derived from these quantities in combination with the velocity three-vector.
