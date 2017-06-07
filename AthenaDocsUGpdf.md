---
title: pdf Files
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/PDF Files

Contain a pdf (probability distribution function) of selected variables as a formatted table (essentially a .tab file).  A variety of
statistics are computed, including the mean and variance of the distribution, the average and standard deviation, and the skewness and
kurtosis of the distribution.

A new sequentially numbered
file is created whenever the integration time exceeds an integer multiple of <output>/dt. At the end of execution, the lesser of tlim/dt or <time>/nlim 
sequentially numbered files will be created. Also writes a 
history-type file with name `prb_stat`.*id*, where *id* is the string identifier for this output, with each row in the file
containing statistics of each .pdf output.

This is a specialized format that is useful
only for some problems, where the statistics of fluctuations in some variable is important (for example, 3D turbulence).

An example of an `<output>` block which creates a pdf of the density is given below.

	<output2>
	out_fmt = pdf                # Binary data dump
	out     = d                  # density
	id      = d                  # string identifier
	dt      = 1.0                # time increment between outputs

Created by the function output_pdf.c; this function should be studied carefully to understand more about this format and what
information it contains.
