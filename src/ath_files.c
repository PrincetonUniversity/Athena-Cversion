#include "copyright.h"
/*==============================================================================
 * FILE: ath_files.c
 *
 * PURPOSE: Functions for creating descriptive output filenames, and for
 *   opening files with such filenames.  Filename has form:
 *        [path]<basename>[.idump][.id].<ext>
 *   where path         optional path
 *         basename     basename of file (usually problem name, e.g. "Sod")
 *         dlen         number of digits to use for numeric extension (1..10)
 *         idump        optional dump number (0,1,2,.....)
 *                      if(dlen > 0 and idump >= 0) we use the dump number
 *                      <idump> uses C-format descriptor "%0<dlen>d"
 *         id           optional additional identifier
 *         ext          file extension, e.g. ".tab", ".bin", ".dx", ".vtk"
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ath_fname()
 *   ath_fopen()
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* ath_fname: creates filename with form [path]<basename>[.idump][.id].<ext>
 *   Used by most of the dump_*() and output_*() functions.
 */

char *ath_fname(const char *path, const char *basename,
                const int dlen, const int idump, 
                const char *id, const char *ext)
{
  char fmt[80];
  size_t size, pathlen = 0;
  char *cp, *fname;

/* 2 = "." following the basename + NUL terminator */
  size = 2 + strlen(basename) + strlen(ext);
  if(path != NULL){
    pathlen = strlen(path);
    size += 1 + pathlen; /* add 1 for possible '/' */
  }
  if(dlen > 0)
    size += (dlen > 10 ? dlen : 10) + 1; /* add 1 for the "." */
  if(id != NULL) size += 1 + strlen(id); /* add 1 for the "." */

  if((fname = (char*)malloc(size*sizeof(char))) == NULL){
    ath_perr(-1,"[ath_fname]: malloc returned a NULL pointer\n");
    return NULL;
  }

/* Build the filename */
  cp = fname;
  /* Start with the optional path */
  if(path != NULL && pathlen > 0){
    strcpy(cp, path);
    cp = &(cp[pathlen]); /* point cp at the '\0' terminator */

    /* Append a '/' if necessary */
    if(cp[-1] != '/'){
      cp[0] = '/';
      cp[1] = '\0';
      cp = &(cp[1]);
    }
  }

/* Append the name of the file */
  if(id){
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s.%%s",dlen);
      sprintf(cp,fmt,basename,idump,id,ext);
    }
    else sprintf(cp,"%s.%s.%s",basename,id,ext);
  }
  else{
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s",dlen);
      sprintf(cp,fmt,basename,idump,ext);
    }
    else sprintf(cp,"%s.%s",basename,ext);
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

FILE *ath_fopen(const char *path, const char *basename,
                const int dlen, const int idump, 
                const char *id, const char *ext, const char *mode)
{
  char fmt[80];
  size_t size, pathlen = 0;
  char *cp, filename[1024];
  FILE *fp;

  /* 2 = "." following the basename + NUL terminator */
  size = 2 + strlen(basename) + strlen(ext);
  if(path != NULL){
    pathlen = strlen(path);
    size += 1 + pathlen; /* add 1 for possible '/' */
  }
  if(dlen > 0)
    size += (dlen > 10 ? dlen : 10) + 1; /* add 1 for the "." */
  if(id != NULL) size += 1 + strlen(id); /* add 1 for the "." */

  /* If the file name is too large, allocate space for it... */
  if(size > 1024){
    char *fname;

    ath_perr(0,"[sim_fopen]: Warning filename length %zu > 1024\n",size);

    fname = ath_fname(path, basename, dlen, idump, id, ext);
    if(fname == NULL){
      ath_perr(-1,"[sim_fopen]: Error allocating a filename\n");
      return NULL;
    }

    fp = fopen(fname, mode);
    free(fname);

    return fp;
  }

/* Build the filename */
  cp = filename;
  /* Start with the optional path */
  if(path != NULL && pathlen > 0){
    strcpy(cp, path);
    cp = &(cp[pathlen]); /* point cp at the '\0' terminator */

    /* Append a '/' if necessary */
    if(cp[-1] != '/'){
      cp[0] = '/';
      cp[1] = '\0';
      cp = &(cp[1]);
    }
  }

/* Append the name of the file */
  if(id){
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s.%%s",dlen);
      sprintf(cp,fmt,basename,idump,id,ext);
    }
    else sprintf(cp,"%s.%s.%s",basename,id,ext);
  }
  else{
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,"%%s.%%0%dd.%%s",dlen);
      sprintf(cp,fmt,basename,idump,ext);
    }
    else sprintf(cp,"%s.%s",basename,ext);
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

