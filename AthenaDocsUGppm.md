---
title: Image Files
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Image Files

Two dimensional images of any variable can be created by Athena in two formats: ppm (portable pixel map) and pgm (portable gray map).  The former
should be used for color images (a corresponding palette must be specified), and the latter for grayscale images.

Since images are more compact than storing the floating point data directly (e.g. as binary or vtk files), the image formats can be used
to create high time resolution animations of any variable.  However, accuracy is lost by not storing the floating point data, therefore image files
cannot be used for further quantitative analysis.

The advantage of using the ppm image format over the many others available (e.g. gif, jpg, tif, etc.) is that no external libraries
are required to create the images.  Moreover, since they are not compressed, there is no loss of information, on the other hand this often makes
ppm images much larger than other formats.

It is easy to post-process the image files created by Athena into other formats, including mpeg (or other format) animations, using
image conversion software such as [ImageMagick](http://www.imagemagick.org).

Both ppm and pgm image files contain a snapshot of the data at a particular time.  New files are created whenever the integration time exceeds an integer multiple of <output>/dt. At the end of execution, the lesser of tlim/dt or <time>/nlim sequentially numbered files will be created.

ppm image files
===============

An example of an <output> block in an input file that generates ppm images of the pressure is given below:

	<output3>
	out_fmt = ppm
	dt      = 0.004
	out     = P
	id      = P
	dmin    = 0.01
	dmax    = 0.70
	palette = rainbow

Note the extra parameters required for ppm image files.  

**dmin** and **dmax** specify the global minimum and maximum applied to all
the images (image colors are scaled to these values).  If these values are not specified, each image is auto-scaled to the min/max at that
time, however this is not generally useful for animations.  Athena will print the min/max over all time for each output variable in a
diagnostic message at the end of execution, so low resolution runs can be used to provide a first guess for the proper values of `dmin/dmax`.

**palette** specifies which color palette to use, examples are given in [Color Palettes in Athena]({{site.baseurl}}/AthenaDocsPalettes).

Slicing is very useful for creating two-dimensional images from 3D runs, see [Specifying Output Slices]({{site.baseurl}}/AthenaDocsUGSlices), and the example
using the pgm format below.

New (user-defined) variables can be output as images, see [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress), and the example using the pgm format below.

pgm image files
===============

Grayscale images can be created using the pgm format, rather than specifying a grayscale palette in the ppm format.  An example of an
<output> block which creates a pgm image of a 2D slice in a 3D calculation using a user-defined variable is given below.

	<output3>
	out_fmt = pgm                # pgm image file
	out     = dVy                # user-defined variable: fluctuations in V3
	id      = dVy                # file id
	usr_expr_flag = 1            # user-defined variable
	dt      = 62.831853          # time step between output of delta V3
	dmin    = -0.0006            # min value for imaging delta V3
	dmax    =  0.0006            # max value for imaging delta V3
	x2      = 0.0                # slice in X-Z plane at Y=0

The variable dVy would have to be calculated in a special function added to the problem generator file.
