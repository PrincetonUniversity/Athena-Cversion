---
title: Non-ideal MHD
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Non-ideal MHD


A variety of non-ideal MHD effects can be added, including Ohmic resistivity, the Hall effect, and ambipolar diffusion.  In all cases,
configure the code with

	% configure --enable-resistivity

Each effect is added at first-order via operator splitting.

The update is explicit in time, so that a very restrictive CFL constraint on the time step will be used.

Ohmic resistivity
=================

Hall effect
===========

Ambipolar diffusion
===================
