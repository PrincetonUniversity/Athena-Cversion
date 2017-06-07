---
title: Output Structure
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Output Structure

A final important data structure defined in `/athena/src/athena.h` is used
for outputs.  Each `<output>` block in an input file creates
a new element in an array whose members are the following structure:

	typedef struct Output_s{
	  int n;          /* the N from the <outputN> block of this output */
	  Real dt;        /* time interval between outputs  */
	  Real t;         /* next time to output */
	  int num;        /* dump number (0=first) */
	  char *out;      /* variable (or user fun) to be output */
	  char *id;       /* filename is of the form <basename>[.idump][.id].<ext> */
	
	/* level and domain number of output (default = [-1,-1] = output all levels) */
	
	  int nlevel, ndomain;
	
	/* variables which describe data min/max */
	  Real dmin,dmax;   /* user defined min/max for scaling data */
	  Real gmin,gmax;   /* computed global min/max (over all output data) */
	  int sdmin,sdmax;  /* 0 = auto scale, otherwise use dmin/dmax */
	
	/* variables which describe coordinates of output data volume */
	  int ndim;       /* 3=cube 2=slice 1=vector 0=scalar */
	  int reduce_x1;  /* flag to denote reduction in x1 (0=no reduction) */
	  int reduce_x2;  /* flag to denote reduction in x2 (0=no reduction) */
	  int reduce_x3;  /* flag to denote reduction in x3 (0=no reduction) */
	  Real x1l, x1u;  /* lower/upper x1 range for data slice  */
	  Real x2l, x2u;  /* lower/upper x2 range for data slice  */
	  Real x3l, x3u;  /* lower/upper x3 range for data slice  */
	
	/* variables which describe output format */
	  char *out_fmt;  /* output format = {bin, tab, hdf, hst, pgm, ppm, ...} */
	  char *dat_fmt;  /* format string for tabular type output, e.g. "%10.5e" */
	  char *palette;  /* name of palette for RGB conversions */
	  float *rgb;     /* array of RGB[256*3] values derived from palette */
	  float *der;     /* helper array of derivatives for interpolation into RGB */
	
	/* pointers to output functions; data expressions */
	  VOutFun_t out_fun; /* output function pointer */
	  VResFun_t res_fun; /* restart function pointer */
	  ConsFun_t expr;     /* pointer to expression that computes quant for output */
	
	}OutputS;

As can be seen, all of the information associated with the output
(variables to be written, format, frequency of output, etc.) are stored
in this structure.  The functions in the file `output.c` should be
consulted for more details.
