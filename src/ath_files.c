#include "copyright.h"
/*==============================================================================
 * FILE: ath_files.c
 *
 * PURPOSE: Functions for creating descriptive output filenames, and for
 *   opening files with such filenames.  Filename has form:
 *        <basename>[.idump][.id].<ext>
 *   where basename     basename of file (usually problem name, e.g. "Sod")
 *         dlen         number of digits to use for numeric extension (1..10)
 *         idump        optional dump number (0,1,2,.....)
 *                      if(dlen > 0 and idump >= 0) we use the dump number
 *                      <idump> uses C-format descriptor "%0<dlen>d"
 *         id           optional additional identifier
 *         ext          file extension, e.g. ".tab", ".bin", ".dx", ".vtk"
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   fname_construct()
 *   ath_fopen()
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* fname_construct: creates filename with form <basename>[.idump][.id].<ext>
 *   Used by most of the dump_*() and output_*() functions.
 */

char *fname_construct(const char *basename, const int dlen, const int idump, 
                      const char *id, const char *ext)
{
  char fmt[80];
  size_t size;
  char *fname;

/* 2 = "." separating basename and extension + NUL terminator */
  size = 2 + strlen(basename) + strlen(ext);
  if(dlen > 0) size += 1 + dlen; /* additional 1 for the "." */
  if(id != NULL) size += 1 + strlen(id); /* additional 1 for the "." */

  if((fname = (char*)malloc(size*sizeof(char))) == NULL){
    ath_perr(-1,"[fname_construct]: malloc returned a NULL pointer\n");
    return NULL;
  }

/* Build the filename */
  if(id){
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s.%%s",dlen);
      sprintf(fname,fmt,basename,idump,id,ext);
    }
    else sprintf(fname,"%s.%s.%s",basename,id,ext);
  }
  else{
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s",dlen);
      sprintf(fname,fmt,basename,idump,ext);
    }
    else sprintf(fname,"%s.%s",basename,ext);
  }

  return fname;
}

/*----------------------------------------------------------------------------*/
/* ath_fopen: open a file for writing with given basename, optional idump 
 *   and identifier, followed by a standard file extension.  It returns the
 *   usual FILE pointer which is used in subsequent I/O routines until the
 *   file is closed with fclose().
 *   Used, for example, for "hst", "tab", "vtk", "pdf", and "rst" dumps/outputs
 */

FILE *ath_fopen(const char *basename, const int dlen, const int idump, 
                const char *id, const char *ext, const char *mode)
{
  char fmt[80];
  char filename[256];
  FILE *fp;

/* Build the filename */
  if(id){
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s.%%s",dlen);
      sprintf(filename,fmt,basename,idump,id,ext);
    }
    else sprintf(filename,"%s.%s.%s",basename,id,ext);
  }
  else{
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s",dlen);
      sprintf(filename,fmt,basename,idump,ext);
    }
    else sprintf(filename,"%s.%s",basename,ext);
  }

  fp = fopen(filename,mode);
  return fp;
}

/*----------------------------------------------------------------------------*/
/* ath_fwrite:  */

size_t ath_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  n = fwrite(ptr,size,nmemb,stream);

  if (n != nmemb) ath_perr(-1,"ath_fwrite: could write enough data\n");
  return n;
}

