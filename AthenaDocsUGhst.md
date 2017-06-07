---
title: Hitory Files
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/History Files

Contain a formatted
table of a variety of volume averaged values.  A new line in the
table is created whenever the integration time exceeds an integer multiple of `<output>/dt`. At the end of execution, the history file will
contain the lesser of tlim/dt or `<time>/nlim` lines, forming a time-history of these
quantities.

The following example shows an `<output>` block in an
input file that generates history files:

        <output2>
        out_fmt = hst            # history dump
        dt      = 0.01           # time increment between outputs


History files contain a text header to label the columns.  An example is given below.

        #   [1]=time      [2]=dt         [3]=mass       [4]=total E    [5]=x1 Mom.    [6]=x2 Mom.    [7]=x3 Mom.    [8]=x1-KE      [9]=x2-KE      [10]=x3-KE   
        #
           0.000000e+00   5.282214e-03   5.625000e-01   1.375000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
           1.207067e-02   2.946638e-03   5.625000e-01   1.375000e+00   1.086360e-02   0.000000e+00   0.000000e+00   3.608664e-03   0.000000e+00   0.000000e+00
           2.081229e-02   2.841588e-03   5.625000e-01   1.375000e+00   1.873106e-02   0.000000e+00   0.000000e+00   6.758749e-03   0.000000e+00   0.000000e+00
           3.217688e-02   2.818091e-03   5.625000e-01   1.375000e+00   2.895919e-02   0.000000e+00   0.000000e+00   1.087849e-02   0.000000e+00   0.000000e+00
           4.063156e-02   2.817390e-03   5.625000e-01   1.375000e+00   3.656840e-02   0.000000e+00   0.000000e+00   1.394250e-02   0.000000e+00   0.000000e+00
           5.191824e-02   2.826396e-03   5.625000e-01   1.375000e+00   4.672641e-02   0.000000e+00   0.000000e+00   1.802746e-02   0.000000e+00   0.000000e+00
           6.040296e-02   2.831914e-03   5.625000e-01   1.375000e+00   5.436267e-02   0.000000e+00   0.000000e+00   2.110439e-02   0.000000e+00   0.000000e+00
           7.174008e-02   2.838001e-03   5.625000e-01   1.375000e+00   6.456607e-02   0.000000e+00   0.000000e+00   2.522883e-02   0.000000e+00   0.000000e+00
           8.025705e-02   2.840966e-03   5.625000e-01   1.375000e+00   7.223135e-02   0.000000e+00   0.000000e+00   2.832463e-02   0.000000e+00   0.000000e+00
           9.162639e-02   2.843850e-03   5.625000e-01   1.375000e+00   8.246375e-02   0.000000e+00   0.000000e+00   3.244672e-02   0.000000e+00   0.000000e+00
           1.001601e-01   2.845360e-03   5.625000e-01   1.375000e+00   9.014413e-02   0.000000e+00   0.000000e+00   3.554578e-02   0.000000e+00   0.000000e+00


Note that successive lines in the history file are not spaced exactly by `<output>/dt` in time, but rather an output occurs
as soon as the current integration time exceeds the next output time (integer multiple of `<output>/dt`).  For example,
if the integration time step is large compared to
`<output>/dt`, an output will occur once every cycle.

History files are created by the function `dump_history.c`.
The data is appended to the file each time the `dump_history()`
function is called.  This means if a history is present from a previous
calculation, **the new data will be appended to the old file.**

New (user-defined) variables can be added as additional columns to history files by following the steps in
[User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress).  The maximum number of new variables that can be added
is controlled by the macro `MAX_USR_H_COUNT` defined in `defs.h.in` (currently 30).
