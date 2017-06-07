---
title: User-defined Output Formats
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/User Formats

It is also fairly easy to add entirely new data output formats, beyond those
described in the [User Guide]({{site.baseurl}}/AthenaDocsUG).  This might include complex data formats that rely on external libraries, like HDF5, etc.,
or simple output like writing data to stdout.

To create a new user-defined output format, follow these steps:

First, write a function that creates the desired output, and add it to the file that contains the 
problem generator.  The function must be of type `void`, and the argument list must contain the `MeshS` and `OutputS`
structures.(For details of the data contained in these
structures, see the [Programmer's Guide]({{site.baseurl}}/AthenaDocsPG).)  As an example, consider a function that outputs
the minimum and maximum of the density on the root (level=0) Domain to stdout.  The following would work:

	static void my_output(MeshS *pM, OutputS *pOut)
	{
	  GridS *pGrid=pM->Domain[0][0].Grid;
	  int nx1,nx2;
	  Real dmin, dmax;
	  Real **data2d=NULL; /* 2D array of data to be dumped */
	
	/* Allocate memory for and compute 2D array of data */
	  data2d = OutData2(pGrid,pOut,&nx1,&nx2);
	  if (data2d == NULL) return; /* data not in range of Grid */
	
	/* Store the global min / max, for output at end of run */
	  minmax2(data2d,nx2,nx1,&dmin,&dmax);
	  pOut->gmin = MIN(dmin,pOut->gmin);
	  pOut->gmax = MAX(dmax,pOut->gmax);
	
	  printf("dmin/dmax = %e %e\n",dmin,dmax);
	
	  free_2d_array(data2d);
	  return;
	}

Next, use the function `get_usr_out_fun()` in the
problem generator file to set a pointer to this new output function if the
string in `<output>/name` has the appropriate value.  As an example, suppose we name 
the new output format we created above *jim*.  Then the `get_usr_out_fun` should 
contain the following:

	VOutFun_t get_usr_out_fun(const char *name)
	{
	  if(strcmp(name,"jim")==0) return my_output;
	  return NULL;
	}

Finally, add an output block to the input file which sets the parameter `name` to the desired string (`jim` in this example)

	<output4>
	name    = jim                # user-defined data dump
	out     = d
	id      = jim
	dt      = 628.31853          # time increment between outputs

When the code is run, the min/max of the density will be printed to stdout at a time interval of `dt`.  Note that the min/max of any other
variable could be specified by changing the parameter `<output4>/out`.  New user-defined variables could also be output by following the
directions in [Adding New Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress).  Indeed, all of the functionality of the output blocks is available
for this new output type.
