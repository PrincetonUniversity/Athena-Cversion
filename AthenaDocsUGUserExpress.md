---
title: User-defined Output Variables
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/User Variables

New Variables added to .tab, .bin, .vtk, or image (.ppm or .pgm) outputs
========================================================================

Often it is useful to output a variable other than the pre-defined values accessed by the
`<output>/out` parameter
(see [Output Blocks]({{site.baseurl}}/AthenaDocsUGOutputBlck)).  For example,
one might like to create ppm images of the current density, J=Curl(B).  This
can be easily accomplished using the `<output>/usr_expr_flag` parameter.  The
following steps are required.

1. Write a function that computes the desired variable, and add it to the file that contains the 
   problem generator.  The function must be of
   type `Real`, and the argument list must contain the `GridS`
   structure and the indices of the grid cell.
   (For details of the data contained in the `GridS`
   structure, see the [Programmer's Guide]({{site.baseurl}}/AthenaDocsPG).)
   The following example, taken from `/athena/src/prob/field_loop.c`, computes the
   z-component of the current density at cell *i,j,k*.

        static Real current(const GridS *pG, const int i, const int j, const int k)
        {  return ((pG->B2i[k][j][i]-pG->B2i[k][j][i-1])/pG->dx1 -
                   (pG->B1i[k][j][i]-pG->B1i[k][j-1][i])/pG->dx2); 
        }

2. Use the function `get_usr_expr()` in the
   problem generator file to set a pointer to this new function if the
   string in `<output>/out` has the appropriate value.  As an example, the following code sets the output variable to be
   the current density computed using the function given above if the
   `<output>/out` string is `J3`

        ConsFun_t get_usr_expr(const char *expr)
        {
          if(strcmp(expr,"J3")==0) return current;
          return NULL;
        }

3. To create a movie of the current density, use an output block in the input
   file in which `<output>/out=J3`, and `<output>/usr_expr_flag=1`,
   with the other valid parameters set as appropriate (to control, for
   example, the time interval between outputs, min/max scaling, etc.), for example

        <output3>
        out_fmt = ppm       # ppm image
        out     = J3
        id      = J3
        usr_expr_flag = 1
        palette = rainbow
        dt      = 0.004     # time step between images
        dmin    = -0.04     # min value for imaging J3
        dmax    =  0.08     # max value for imaging J3

Alternatively, one could output the current density as floating-point numbers in vtk format using `out_fmt=vtk`, or as a formatted table using
`out_fmt=tab`.

New Variables added to .hst dumps
=================================

It is also possible to add new volume-averaged quantities that are output as additional columns in history (.hst) dumps.
To add a new, user-defined history variable, use the following steps.

1. Write a function that computes the desired variable, and add it to the file that contains the 
   problem generator, identical to step 1 above for all other output types.

2. Enroll this new function by adding a call to `dump_history_enroll()` anywhere in the problem generator.  This function has two arguments,
   the first is the name of the function created in step 1, the second is a string used as a header for the new column in the history file.  For example,
   the following four lines of code in the problem generator enroll four new history variables.

        dump_history_enroll(hst_Bx, "<Bx>");
        dump_history_enroll(hst_By, "<By>");
        dump_history_enroll(hst_Bz, "<Bz>");
        dump_history_enroll(hst_BxBy, "<-Bx By>");

These variables are computed by the user-defined functions `hst_Bx()`, `hst_By()`, `hst_Bz()`, and `hst_BxBy()` respectively,
and create new columns with headers `<Bx>`, `<By>`, `<Bz>`, and `<-Bx By>` respectively.
