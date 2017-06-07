---
title: Editing the Input File
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[Tutorial]({{site.baseurl}}/AthenaDocsTut)/Input File

The previous two examples (Sod shocktube in 1D hydrodynamics, and Orszag-Tang vortex in 2D MHD)
used the default input files in the `/athena/tst` directories.  Suppose, however, one would like
to change the grid resolution, or output files, or any other parameter of the runs.  This is accomplished
by editing the input files.

Changing the Grid Resolution
============================

For example, the grid resolution in the 2D Orszag-Tang vortex test can be changed by editing the input
file using the following steps.

1. First, make sure there is an executable for the Orszag-Tang problem in the `/athena/bin` directory.  If necessary,
   repeat step 1 from the [2D MHD]({{site.baseurl}}/AthenaDocsTutOT1) tutorial

        % make clean
        % configure --with-problem=orszag-tang --with-order=3
        % make all


2. Next, copy the input file into the `/bin` directory.

        % cd bin
        % cp ../tst/2D-mhd/athinput.orszag-tang  athinput.new


3. Now, edit the new input file to change the grid resolution.  This is set by the `Nx1` and `Nx2` parameters in the `<domain1>`
   block.  Therefore, change the lines:

        ...
        Nx1             = 192       # Number of zones in X1-direction
        ...
        Nx2             = 192       # Number of zones in X2-direction
        ...

    to any other values you like.  Make the resolution larger will make the code run slower, so try *decreasing* the resolution to 

        ...
        Nx1             = 128       # Number of zones in X1-direction
        ...
        Nx2             = 128       # Number of zones in X2-direction
        ...


4. Run the code using this new input file

        % athena -i athinput.new

    This should generate the same data files as before, but on a grid with resolution 128^2^.  The code should have run about four times 
    faster compared to the default resolution of 192^2^

Adding a New Output
===================

Suppose you would like to output the data in the 2D Orszag-Tang vortex test in vtk, instead of binary, format, and you would
like to add output to make a movie of the magnetic energy.  You need to edit existing, and add new, `<output>` blocks in the input file.
Use the following steps.

1. If you haven't already, configure and compile the code, and copy the input file to the `/bin` directory.

        % make clean
        % configure --with-problem=orszag-tang --with-order=3
        % make all
        % cd bin
        % cp ../tst/2D-mhd/athinput.orszag-tang  athinput.new


2. Next, edit the `<output2>` block in the input file to change the output format from:

        out_fmt = bin               # Binary data dump

    to

        out_fmt = vtk               # vtk data dump


3. Add a new `<output>` block for ppm images of the magnetic energy.  A useful time interval is one that creates several hundred images
   (so the time evolution is smooth).  Since we don't know what the min/max of the magnetic energy will be, we can use autoscaling of the images.
   Thus, the new output block might be

        <output5>
        out_fmt = ppm
        dt      = 0.004
        out     = ME
        id      = ME
        palette = heat


4. Now increase the number of output blocks to be read, by changing the `maxout` parameter in the `<job>` block to

        maxout       = 5            # Read output blocks numbered from 1 -> maxout


5. Now run the code.  The last few lines of output (when run on a 192^2^ grid) should be something like

        cycle=655 time=9.999868e-01 next dt=1.320914e-05 last dt=1.449403e-03
        cycle=656 time=1.000000e+00 next dt=0.000000e+00 last dt=1.320914e-05
        
        terminating on time limit
          tlim= 1.000000e+00   nlim= 100000
          time= 1.000000e+00  cycle= 656
        
        zone-cycles/cpu-second = 2.251027e+05
        
        elapsed wall time = 1.074434e+02 sec.
        
        zone-cycles/wall-second = 2.250746e+05
        Global min/max for P: 0.0100164 0.755822
        Global min/max for d: 0.0442782 0.609765
        Global min/max for ME: 8.51744e-08 0.411589
        
        Simulation terminated on Thu Apr 29 11:56:58 2010

    Note the global min/max of ME are reported.  If you animate the `OrszagTang.*.ME.ppm` images, they will seem to flicker
    (especially near the beginning, where the magnetic energy is evolving rapidly) since each
    image uses a different scaling (min/max).  So edit the `<output5>` block to add appropriate values for the min/max of ME
    reported above:

        <output5>
        out_fmt = ppm
        dt      = 0.004
        out     = ME
        id      = ME
        dmin    = 0.0
        dmax    = 0.4
        palette = heat

    Run the code again, and now the movie of the ME should be smooth.  Try editing the `<output5>` block to use different
    palettes, and see what the resulting images are like.

6. Note the code no longer outputs `.bin` files, but these have been replaced by `.vtk` files.  Athena vtk files can be
   read directly by [VisIt](http://www.llnl.gov/visit).  Try plotting the density, pressure, velocity, etc using VisIt to read the
   .vtk files.  You can also make movies of any quantity you like from the files using VisIt.

Try modifying other parameters of the runs by editing the input files.
