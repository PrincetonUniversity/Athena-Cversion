#include "copyright.h"
/*==============================================================================
 * FILE: dump_dx.c
 *
 * PURPOSE: Function to write dumps that can be read by OpenDX.  DX dumps
 *   consist of two files, 'filename.dx' and 'filename.bin'.  This function
 *   writes the first (basically a header file with information about the grid),
 *   and calls dump_binary() to write the second (the actual data in binary
 *   format).  The filenames must have the same consecutive numbering system.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   dump_dx -
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* dump_dx:    */

void dump_dx(Grid *pGrid, Output *pOut)
{
  int n, offset, dnum = pOut->num; 
  FILE *pfile;
  int nzones, nx1, nx2, nx3;
  char *fname;
  char *array_name[] = 
    {"Mass Density",
     "1-Momentum",
     "2-Momentum",
     "3-Momentum"
#ifndef ISOTHERMAL
     ,"Energy Density"
#endif
#ifdef MHD
     ,"1-Field","2-Field","3-Field"
#endif
    };

#ifdef WRITE_GHOST_CELLS
  nx1 = pGrid->Nx1 > 1 ? pGrid->Nx1 + 2*nghost : pGrid->Nx1;
  nx2 = pGrid->Nx2 > 1 ? pGrid->Nx2 + 2*nghost : pGrid->Nx2;
  nx3 = pGrid->Nx3 > 1 ? pGrid->Nx3 + 2*nghost : pGrid->Nx3;
#else
  nx1 = pGrid->Nx1;
  nx2 = pGrid->Nx2;
  nx2 = pGrid->Nx3;
#endif
  nzones = nx1*nx2*nx3;

/* call dump_binary() to output data in corresponding .bin file */
  dump_binary(pGrid, pOut);

/* reconstruct binary dump filename, and use to construct .dx filename */
  if((fname = fname_construct(pGrid->outfilename,num_digit,dnum,NULL,"bin")) 
     == NULL){
    ath_error("[dump_binary]: No filename -- no binary dump\n");
    return;
  }
  if ((pfile = ath_fopen(fname,0,0,NULL,"dx","w")) == NULL) {
    ath_error("[write_dx_header]: Unable to open dx header file\n");
    return; 
  }

/* start writing descriptive information to .dx file */

  fprintf(pfile,"# The default data type is ");
  if(ath_big_endian()){
    fprintf(pfile,"most-significant-byte first and binary\n");
    fprintf(pfile,"data mode msb binary\n\n");
  }
  else{
    fprintf(pfile,"least-significant-byte first and binary\n");
    fprintf(pfile,"data mode lsb binary\n\n");
  }

  fprintf(pfile,"# object 1 is the grid cell corner positions\n");
  fprintf(pfile,"object 1 class gridpositions counts %d %d %d\n",
          nx3+1,nx2+1,nx1+1);
  fprintf(pfile,"origin 0 0 0\n");
  fprintf(pfile,"delta  0 0 1\n");
  fprintf(pfile,"delta  0 1 0\n");
  fprintf(pfile,"delta  1 0 0\n\n");

  fprintf(pfile,"# object 2 are the regular connections, quads\n");
  fprintf(pfile,"object 2 class gridconnections counts %d %d %d\n\n",
          nx3+1,nx2+1,nx1+1);

#ifdef MHD
  fprintf(pfile,"# object 3 is the cell interface positions\n");
  fprintf(pfile,"object 3 class gridpositions counts %d %d %d\n",nx3,nx2,nx1);
  fprintf(pfile,"origin 0 0 0.5\n");
  fprintf(pfile,"delta  0 0 1\n");
  fprintf(pfile,"delta  0 1 0\n");
  fprintf(pfile,"delta  1 0 0\n\n");

  fprintf(pfile,"# object 4 is the cell interface positions\n");
  fprintf(pfile,"object 4 class gridpositions counts %d %d %d\n",nx3,nx2,nx1);
  fprintf(pfile,"origin 0 0.5 0\n");
  fprintf(pfile,"delta  0 0 1\n");
  fprintf(pfile,"delta  0 1 0\n");
  fprintf(pfile,"delta  1 0 0\n\n");

  fprintf(pfile,"# object 5 is the cell interface positions\n");
  fprintf(pfile,"object 5 class gridpositions counts %d %d %d\n",nx3,nx2,nx1);
  fprintf(pfile,"origin 0.5 0 0\n");
  fprintf(pfile,"delta  0 0 1\n");
  fprintf(pfile,"delta  0 1 0\n");
  fprintf(pfile,"delta  1 0 0\n\n");
#endif /* MHD */

  fprintf(pfile,"# objects 6 on are the data, which ");
  fprintf(pfile,"are in a one-to-one correspondence\n");
  fprintf(pfile,"# with the connections (\"dep\" on connections).\n\n");

/* Initialize the byte offset to skip over the first 2 fwrites 
 * (3 ints and 2 floats) */
  offset = 3*sizeof(int) + 2*sizeof(float);

  for(n=0; n<NVAR; n++){
    fprintf(pfile,"# %s\n",array_name[n]);
    fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	    n+5,nzones);
    fprintf(pfile,"data file \"%s\",%d\n",fname,offset);
    fprintf(pfile,"attribute \"dep\" string \"connections\"\n\n");
    offset += nzones*sizeof(float);
  }

  for(n=0; n<72; n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Construct a field for each of the objects above.\n");

  for(n=0; n<NVAR; n++){
    fprintf(pfile,"object \"%s\" class field\n",array_name[n]);
    fprintf(pfile,"component \"positions\" value 1\n");
    fprintf(pfile,"component \"connections\" value 2\n");
    fprintf(pfile,"component \"data\" value %d\n\n",n+5);
  }

  for(n=0; n<72; n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Construct a string object");
  fprintf(pfile," which stores the names of the fields.\n");
  fprintf(pfile,"# This way we can select an individual field");
  fprintf(pfile," before importing the data.\n\n");

  fprintf(pfile,"object \"Variables\" class string\n");
  for(n=0;n<NVAR;n++){
    fprintf(pfile," \"%s\"\n",array_name[n]);
  }

  fprintf(pfile,"\nend\n");

  fclose(pfile);
  return;
}
