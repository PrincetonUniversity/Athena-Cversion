#include "copyright.h"
/*==============================================================================
 * FILE: ath_array.c
 *
 * PURPOSE: Functions that construct and destruct 1D, 2D and 3D arrays of any
 *   type. Elements in these arrays can be accessed in the standard fashion of
 *   array[in][im] (where in = 0 -> nr-1 and im = 0 -> nc-1) as if it
 *   were an array of fixed size at compile time with the statement:
 *      "Type" array[nr][nc];
 *  Equally so for 3D arrys:  array[nt][nr][nc]       -- TAG -- 8/2/2001
 *
 * EXAMPLE usage of 3D construct/destruct functions for arrays of type Real:
 *   array = (Real ***)calloc_3d_array(nt,nr,nc,sizeof(Real));
 *   free_3d_array(array);
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   calloc_1d_array() - creates 1D array
 *   calloc_2d_array() - creates 2D array
 *   calloc_3d_array() - creates 3D array
 *   calloc_4d_array() - creates 4D array
 *   calloc_5d_array() - creates 5D array
 *   calloc_6d_array() - creates 6D array
 *   calloc_7d_array() - creates 7D array
 *   free_1d_array()   - destroys 1D array
 *   free_2d_array()   - destroys 2D array
 *   free_3d_array()   - destroys 3D array
 *   free_4d_array()   - destroys 4D array
 *   free_5d_array()   - destroys 5D array
 *   free_6d_array()   - destroys 6D array
 *   free_7d_array()   - destroys 7D array
 *============================================================================*/

#include <stdlib.h>
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* calloc_1d_array: construct 1D array = array[nc]  */

void* calloc_1d_array(size_t nc, size_t size)
{
  void *array;
  
  if ((array = (void *)calloc(nc,size)) == NULL) {
    ath_error("[calloc_1d] failed to allocate memory (%d of size %d)\n",
              (int)nc,(int)size);
    return NULL;
  }
  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_2d_array: construct 2D array = array[nr][nc]  */

void** calloc_2d_array(size_t nr, size_t nc, size_t size)
{
  void **array;
  size_t i;

  if((array = (void **)calloc(nr,sizeof(void*))) == NULL){
    ath_error("[calloc_2d] failed to allocate mem for %d pointers\n",(int)nr);
    return NULL;
  }

  if((array[0] = (void *)calloc(nr*nc,size)) == NULL){
    ath_error("[calloc_2d] failed to allocate memory (%d X %d of size %d)\n",
              (int)nr,(int)nc,(int)size);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nr; i++){
    array[i] = (void *)((unsigned char *)array[0] + i*nc*size);
  }

  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_3d_array: construct 3D array = array[nt][nr][nc]  */

void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size)
{
  void ***array;
  size_t i,j;

  if((array = (void ***)calloc(nt,sizeof(void**))) == NULL){
    ath_error("[calloc_3d] failed to allocate memory for %d 1st-pointers\n",
              (int)nt);
    return NULL;
  }

  if((array[0] = (void **)calloc(nt*nr,sizeof(void*))) == NULL){
    ath_error("[calloc_3d] failed to allocate memory for %d 2nd-pointers\n",
              (int)(nt*nr));
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nt; i++){
    array[i] = (void **)((unsigned char *)array[0] + i*nr*sizeof(void*));
  }

  if((array[0][0] = (void *)calloc(nt*nr*nc,size)) == NULL){
    ath_error("[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n",
              (int)nt,(int)nr,(int)nc,(int)size);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(j=1; j<nr; j++){
    array[0][j] = (void **)((unsigned char *)array[0][j-1] + nc*size);
  }

  for(i=1; i<nt; i++){
    array[i][0] = (void **)((unsigned char *)array[i-1][0] + nr*nc*size);
    for(j=1; j<nr; j++){
      array[i][j] = (void **)((unsigned char *)array[i][j-1] + nc*size);
    }
  }

  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_4d_array: construct 4D array = array[ni][nj][nk][nl]  */

void**** calloc_4d_array(size_t ni, size_t nj, size_t nk, size_t nl, size_t size)
{
  void ****array;
  size_t i,j,k;

  if((array = (void ****)calloc(ni,sizeof(void***))) == NULL){
    fprintf(stderr,"[calloc_4d] failed to allocate memory for %d 1st-pointers\n",
              (int)ni);
    return NULL;
  }

  if((array[0] = (void ***)calloc(ni*nj,sizeof(void**))) == NULL){
    fprintf(stderr,"[calloc_4d] failed to allocate memory for %d 2nd-pointers\n",
              (int)ni*nj);
    free((void *)array);
    return NULL;
  }

  if((array[0][0] = (void **)calloc(ni*nj*nk,sizeof(void*))) == NULL){
    fprintf(stderr,"[calloc_4d] failed to allocate memory for %d 3rd-pointers\n",
              (int)ni*nj*nk);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0] = (void *)calloc(ni*nj*nk*nl,size)) == NULL){
    fprintf(stderr,"[calloc_4d] failed to alloc. memory (%d X %d X %d X %d of size %d)\n",
	    (int)ni,(int)nj,(int)nk,(int)nl,(int)size);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<ni; i++){
    array[i] = (void ***)((unsigned char *)array[i-1] + nj*sizeof(void**));
  }

  for(j=1; j<nj; j++){
    array[0][j] = (void **)((unsigned char *)array[0][j-1] + nk*sizeof(void*));
  }

  for(i=1; i<ni; i++){
    array[i][0] = (void **)((unsigned char *)array[i-1][0] + nj*nk*sizeof(void*));
    
    for(j=1; j<nj; j++){
      array[i][j] = (void **)((unsigned char *)array[i][j-1] + nk*sizeof(void*));
    }
  }

  for(k=1; k<nk; k++){
    array[0][0][k] = (void *)((unsigned char *)array[0][0][k-1] + nl*size);
  }

  for(j=1; j<nj; j++){
    array[0][j][0] = (void *)((unsigned char *)array[0][j-1][0] + nk*nl*size);
    
    for(k=1; k<nk; k++){
      array[0][j][k] = (void *)((unsigned char *)array[0][j][k-1] + nl*size);
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0] = (void *)((unsigned char *)array[i-1][0][0] + nj*nk*nl*size);
    
    for(k=1; k<nk; k++)
      array[i][0][k] = (void *)((unsigned char *)array[i][0][k-1] + nl*size);
    
    for(j=1; j<nj; j++){
      array[i][j][0] = (void *)((unsigned char *)array[i][j-1][0] + nk*nl*size);

      for(k=1; k<nk; k++)
	array[i][j][k] = (void *)((unsigned char *)array[i][j][k-1] + nl*size);      
    }
  }


  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_5d_array: construct 5D array = array[ni][nj][nk][nl][nm]  */

void***** calloc_5d_array(size_t ni, size_t nj, size_t nk, size_t nl, size_t nm,
			  size_t size)
{
  void *****array;
  size_t i,j,k,l;

  if((array = (void *****)calloc(ni,sizeof(void****))) == NULL){
    fprintf(stderr,"[calloc_5d] failed to allocate memory for %d 1st-pointers\n",
              (int)ni);
    return NULL;
  }

  if((array[0] = (void ****)calloc(ni*nj,sizeof(void***))) == NULL){
    fprintf(stderr,"[calloc_5d] failed to allocate memory for %d 2nd-pointers\n",
              (int)ni*nj);
    free((void *)array);
    return NULL;
  }

  if((array[0][0] = (void ***)calloc(ni*nj*nk,sizeof(void**))) == NULL){
    fprintf(stderr,"[calloc_5d] failed to allocate memory for %d 3rd-pointers\n",
              (int)ni*nj*nk);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0] = (void **)calloc(ni*nj*nk*nl,sizeof(void*))) == NULL){
    fprintf(stderr,"[calloc_5d] failed to allocate memory for %d 4th-pointers\n",
              (int)(ni*nj*nk*nl));
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0] = (void *)calloc(ni*nj*nk*nl*nm,size)) == NULL){
    fprintf(stderr,"[calloc_5d] failed to alloc. memory (%d X %d X %d X %d X %d of size %d)\n",
	    (int)ni,(int)nj,(int)nk,(int)nl,(int)nm,(int)size);
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<ni; i++){
    array[i] = (void ****)((unsigned char *)array[i-1] + nj*sizeof(void***));
  }

  for(j=1; j<nj; j++){
    array[0][j] = (void ***)((unsigned char *)array[0][j-1] + nk*sizeof(void**));
  }

  for(i=1; i<ni; i++){
    array[i][0] = (void ***)((unsigned char *)array[i-1][0] + nj*nk*sizeof(void**));
    for(j=1; j<nj; j++){
      array[i][j] = (void ***)((unsigned char *)array[i][j-1] + nk*sizeof(void**));
    }
  }

  for(k=1; k<nk; k++){
    array[0][0][k] = (void **)((unsigned char *)array[0][0][k-1] + nl*sizeof(void*));
  }

  for(j=1; j<nj; j++){
    array[0][j][0] = (void **)((unsigned char *)array[0][j-1][0] + nk*nl*sizeof(void*));
    for(k=1; k<nk; k++){
      array[0][j][k] = (void **)((unsigned char *)array[0][j][k-1] + nl*sizeof(void*));
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0] = (void **)((unsigned char *)array[i-1][0][0] + nj*nk*nl*sizeof(void*));
    
    for(k=1; k<nk; k++)
      array[i][0][k] = (void **)((unsigned char *)array[i][0][k-1] + nl*sizeof(void*));
    
    for(j=1; j<nj; j++){
      array[i][j][0] = (void **)((unsigned char *)array[i][j-1][0] + nk*nl*sizeof(void*));

      for(k=1; k<nk; k++)
	array[i][j][k] = (void **)((unsigned char *)array[i][j][k-1] + nl*sizeof(void*));      
    }
  }

  for(l=1; l<nl; l++)
    array[0][0][0][l] = (void *)((unsigned char *)array[0][0][0][l-1] + nm*size);
  

  for(k=1; k<nk; k++){
    array[0][0][k][0] = (void *)((unsigned char *)array[0][0][k-1][0] + nl*nm*size);

    for(l=1; l<nl; l++)
      array[0][0][k][l] = (void *)((unsigned char *)array[0][0][k][l-1] + nm*size);    
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0] = (void *)((unsigned char *)array[0][j-1][0][0] + nk*nl*nm*size);

    for(l=1; l<nl; l++)
      array[0][j][0][l] = (void *)((unsigned char *)array[0][j][0][l-1] + nm*size);
    
    for(k=1; k<nk; k++){
      array[0][j][k][0] = (void *)((unsigned char *)array[0][j][k-1][0] + nl*nm*size);

      for(l=1; l<nl; l++)
	array[0][j][k][l] = (void *)((unsigned char *)array[0][j][k][l-1] + nm*size);      
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0] = (void *)((unsigned char *)array[i-1][0][0][0] + nj*nk*nl*nm*size);

    for(l=1; l<nl; l++)
      array[i][0][0][l] = (void *)((unsigned char *)array[i][0][0][l-1] + nm*size);

    for(k=1; k<nk; k++){
      array[i][0][k][0] = (void *)((unsigned char *)array[i][0][k-1][0] + nl*nm*size);
      
      for(l=1; l<nl; l++)
	array[i][0][k][l] = (void *)((unsigned char *)array[i][0][k][l-1] + nm*size);   
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0] = (void *)((unsigned char *)array[i][j-1][0][0] + nk*nl*nm*size);

      for(l=1; l<nl; l++)
	array[i][j][0][l] = (void *)((unsigned char *)array[i][j][0][l-1] + nm*size);

      for(k=1; k<nk; k++){
	array[i][j][k][0] = (void *)((unsigned char *)array[i][j][k-1][0] + nl*nm*size);
	
	for(l=1; l<nl; l++)
	  array[i][j][k][l] = (void *)((unsigned char *)array[i][j][k][l-1] + nm*size);
      }
    }
  }


  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_6d_array: construct 6D array = array[ni][nj][nk][nl][nm][nn]  */

void****** calloc_6d_array(size_t ni, size_t nj, size_t nk, size_t nl, size_t nm,
			   size_t nn, size_t size)
{
  void ******array;
  size_t i,j,k,l,m;

  if((array = (void ******)calloc(ni,sizeof(void*****))) == NULL){
    fprintf(stderr,"[calloc_6d] failed to allocate memory for %d 1st-pointers\n",
              (int)ni);
    return NULL;
  }

  if((array[0] = (void *****)calloc(ni*nj,sizeof(void****))) == NULL){
    fprintf(stderr,"[calloc_6d] failed to allocate memory for %d 2nd-pointers\n",
              (int)ni*nj);
    free((void *)array);
    return NULL;
  }

  if((array[0][0] = (void ****)calloc(ni*nj*nk,sizeof(void***))) == NULL){
    fprintf(stderr,"[calloc_6d] failed to allocate memory for %d 3rd-pointers\n",
              (int)ni*nj*nk);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0] = (void ***)calloc(ni*nj*nk*nl,sizeof(void**))) == NULL){
    fprintf(stderr,"[calloc_6d] failed to allocate memory for %d 4th-pointers\n",
              (int)(ni*nj*nk*nl));
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0] = (void **)calloc(ni*nj*nk*nl*nm,sizeof(void*))) == NULL){
    fprintf(stderr,"[calloc_6d] failed to allocate memory for %d 5th-pointers\n",
              (int)(ni*nj*nk*nl*nm));
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0][0] = (void *)calloc(ni*nj*nk*nl*nm*nn,size)) == NULL){
    fprintf(stderr,"[calloc_6d] failed to alloc. memory (%d X %d X %d X %d X %d X %d of size %d)\n",
	    (int)ni,(int)nj,(int)nk,(int)nl,(int)nm,(int)nn,(int)size);
    free((void *)array[0][0][0][0]);
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<ni; i++){
    array[i] = (void *****)((unsigned char *)array[i-1] + nj*sizeof(void****));
  }

  for(j=1; j<nj; j++){
    array[0][j] = (void ****)((unsigned char *)array[0][j-1] + nk*sizeof(void***));
  }

  for(i=1; i<ni; i++){
    array[i][0] = (void ****)((unsigned char *)array[i-1][0] + nj*nk*sizeof(void***));
    for(j=1; j<nj; j++){
      array[i][j] = (void ****)((unsigned char *)array[i][j-1] + nk*sizeof(void***));
    }
  }

  for(k=1; k<nk; k++){
    array[0][0][k] = (void ***)((unsigned char *)array[0][0][k-1] + nl*sizeof(void**));
  }

  for(j=1; j<nj; j++){
    array[0][j][0] = (void ***)((unsigned char *)array[0][j-1][0] + nk*nl*sizeof(void**));
    for(k=1; k<nk; k++){
      array[0][j][k] = (void ***)((unsigned char *)array[0][j][k-1] + nl*sizeof(void**));
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0] = (void ***)((unsigned char *)array[i-1][0][0] + nj*nk*nl*sizeof(void**));
    
    for(k=1; k<nk; k++)
      array[i][0][k] = (void ***)((unsigned char *)array[i][0][k-1] + nl*sizeof(void**));
    
    for(j=1; j<nj; j++){
      array[i][j][0] = (void ***)((unsigned char *)array[i][j-1][0] + nk*nl*sizeof(void**));

      for(k=1; k<nk; k++)
	array[i][j][k] = (void ***)((unsigned char *)array[i][j][k-1] + nl*sizeof(void**));      
    }
  }


  for(l=1; l<nl; l++)
    array[0][0][0][l] = (void **)((unsigned char *)array[0][0][0][l-1] + nm*sizeof(void*));
  

  for(k=1; k<nk; k++){
    array[0][0][k][0] = (void **)((unsigned char *)array[0][0][k-1][0] + nl*nm*sizeof(void*));

    for(l=1; l<nl; l++)
      array[0][0][k][l] = (void **)((unsigned char *)array[0][0][k][l-1] + nm*sizeof(void*));    
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0] = (void **)((unsigned char *)array[0][j-1][0][0] + nk*nl*nm*sizeof(void*));

    for(l=1; l<nl; l++)
      array[0][j][0][l] = (void **)((unsigned char *)array[0][j][0][l-1] + nm*sizeof(void*));
    
    for(k=1; k<nk; k++){
      array[0][j][k][0] = (void **)((unsigned char *)array[0][j][k-1][0] + nl*nm*sizeof(void*));

      for(l=1; l<nl; l++)
	array[0][j][k][l] = (void **)((unsigned char *)array[0][j][k][l-1] + nm*sizeof(void*));      
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0] = (void **)((unsigned char *)array[i-1][0][0][0] + nj*nk*nl*nm*sizeof(void*));

    for(l=1; l<nl; l++)
      array[i][0][0][l] = (void **)((unsigned char *)array[i][0][0][l-1] + nm*sizeof(void*));

    for(k=1; k<nk; k++){
      array[i][0][k][0] = (void **)((unsigned char *)array[i][0][k-1][0] + nl*nm*sizeof(void*));
      
      for(l=1; l<nl; l++)
	array[i][0][k][l] = (void **)((unsigned char *)array[i][0][k][l-1] + nm*sizeof(void*));   
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0] = (void **)((unsigned char *)array[i][j-1][0][0] + nk*nl*nm*sizeof(void*));

      for(l=1; l<nl; l++)
	array[i][j][0][l] = (void **)((unsigned char *)array[i][j][0][l-1] + nm*sizeof(void*));

      for(k=1; k<nk; k++){
	array[i][j][k][0] = (void **)((unsigned char *)array[i][j][k-1][0] + nl*nm*sizeof(void*));
	
	for(l=1; l<nl; l++)
	  array[i][j][k][l] = (void **)((unsigned char *)array[i][j][k][l-1] + nm*sizeof(void*));
      }
    }
  }

  for(m=1; m<nm; m++)
    array[0][0][0][0][m] = (void *)((unsigned char *)array[0][0][0][0][m-1] + nn*size);


  for(l=1; l<nl; l++) {
    array[0][0][0][l][0] = (void *)((unsigned char *)array[0][0][0][l-1][0] + nm*nn*size);

    for(m=1; m<nm; m++)
      array[0][0][0][l][m] = (void *)((unsigned char *)array[0][0][0][l][m-1] + nn*size);
  }

  for(k=1; k<nk; k++){
    array[0][0][k][0][0] = (void *)((unsigned char *)array[0][0][k-1][0][0] + nl*nm*nn*size);

    for(m=1; m<nm; m++)
      array[0][0][k][0][m] = (void *)((unsigned char *)array[0][0][k][0][m-1] + nn*size);

    for(l=1; l<nl; l++) {
      array[0][0][k][l][0] = (void *)((unsigned char *)array[0][0][k][l-1][0] + nm*nn*size);   

      for(m=1; m<nm; m++)
	array[0][0][k][l][m] = (void *)((unsigned char *)array[0][0][k][l][m-1] + nn*size);  	
    } 
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0][0] = (void *)((unsigned char *)array[0][j-1][0][0][0] + nk*nl*nm*nn*size);

    for(m=1; m<nm; m++)
      array[0][j][0][0][m] = (void *)((unsigned char *)array[0][j][0][0][m-1] + nn*size);

    for(l=1; l<nl; l++) {
      array[0][j][0][l][0] = (void *)((unsigned char *)array[0][j][0][l-1][0] + nm*nn*size);
    
      for(m=1; m<nm; m++)
	array[0][j][0][l][m] = (void *)((unsigned char *)array[0][j][0][l][m-1] + nn*size);
    }

    for(k=1; k<nk; k++){
      array[0][j][k][0][0] = (void *)((unsigned char *)array[0][j][k-1][0][0] + nl*nm*nn*size);

      for(m=1; m<nm; m++)
	array[0][j][k][0][m] = (void *)((unsigned char *)array[0][j][k][0][m-1] + nn*size);

      for(l=1; l<nl; l++) {
	array[0][j][k][l][0] = (void *)((unsigned char *)array[0][j][k][l-1][0] + nm*nn*size);    

	for(m=1; m<nm; m++)
	  array[0][j][k][l][m] = (void *)((unsigned char *)array[0][j][k][l][m-1] + nn*size);  
      } 
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0][0] = (void *)((unsigned char *)array[i-1][0][0][0][0] + nj*nk*nl*nm*nn*size);

    for(m=1; m<nm; m++)
      array[i][0][0][0][m] = (void *)((unsigned char *)array[i][0][0][0][m-1] + nn*size);

    for(l=1; l<nl; l++) {
      array[i][0][0][l][0] = (void *)((unsigned char *)array[i][0][0][l-1][0] + nm*nn*size);

      for(m=1; m<nm; m++)
	array[i][0][0][l][m] = (void *)((unsigned char *)array[i][0][0][l][m-1] + nn*size);
    }

    for(k=1; k<nk; k++){
      array[i][0][k][0][0] = (void *)((unsigned char *)array[i][0][k-1][0][0] + nl*nm*nn*size);

      for(m=1; m<nm; m++)
	array[i][0][k][0][m] = (void *)((unsigned char *)array[i][0][k][0][m-1] + nn*size);

      for(l=1; l<nl; l++) {
	array[i][0][k][l][0] = (void *)((unsigned char *)array[i][0][k][l-1][0] + nm*nn*size);    

	for(m=1; m<nm; m++)
	  array[i][0][k][l][m] = (void *)((unsigned char *)array[i][0][k][l][m-1] + nn*size);  
      } 
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0][0] = (void *)((unsigned char *)array[i][j-1][0][0][0] + nk*nl*nm*nn*size);

      for(m=1; m<nm; m++)
	array[i][j][0][0][m] = (void *)((unsigned char *)array[i][j][0][0][m-1] + nn*size);
      
      for(l=1; l<nl; l++) {
	array[i][j][0][l][0] = (void *)((unsigned char *)array[i][j][0][l-1][0] + nm*nn*size);
    
	for(m=1; m<nm; m++)
	  array[i][j][0][l][m] = (void *)((unsigned char *)array[i][j][0][l][m-1] + nn*size);
      }

      for(k=1; k<nk; k++){
	array[i][j][k][0][0] = (void *)((unsigned char *)array[i][j][k-1][0][0] + nl*nm*nn*size);

	for(m=1; m<nm; m++)
	  array[i][j][k][0][m] = (void *)((unsigned char *)array[i][j][k][0][m-1] + nn*size);
	
	for(l=1; l<nl; l++) {
	  array[i][j][k][l][0] = (void *)((unsigned char *)array[i][j][k][l-1][0] + nm*nn*size);    
	  
	  for(m=1; m<nm; m++)
	    array[i][j][k][l][m] = (void *)((unsigned char *)array[i][j][k][l][m-1] + nn*size);  
	} 
      }
    }
  }


  return array;
}

/*----------------------------------------------------------------------------*/
/* calloc_7d_array: construct 7D array = array[ni][nj][nk][nl][nm][nn][no]  */

void******* calloc_7d_array(size_t ni, size_t nj, size_t nk, size_t nl, size_t nm, 
			    size_t nn, size_t no, size_t size)
{
  void *******array;
  size_t i,j,k,l,m,n;

  if((array = (void *******)calloc(ni,sizeof(void******))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 1st-pointers\n",
              (int)ni);
    return NULL;
  }

  if((array[0] = (void ******)calloc(ni*nj,sizeof(void*****))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 2nd-pointers\n",
              (int)ni*nj);
    free((void *)array);
    return NULL;
  }

  if((array[0][0] = (void *****)calloc(ni*nj*nk,sizeof(void****))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 3rd-pointers\n",
              (int)ni*nj*nk);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0] = (void ****)calloc(ni*nj*nk*nl,sizeof(void***))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 4th-pointers\n",
              (int)(ni*nj*nk*nl));
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0] = (void ***)calloc(ni*nj*nk*nl*nm,sizeof(void**))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 5th-pointers\n",
              (int)(ni*nj*nk*nl*nm));
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0][0] = (void **)calloc(ni*nj*nk*nl*nm*nn,sizeof(void*))) == NULL){
    fprintf(stderr,"[calloc_7d] failed to allocate memory for %d 6th-pointers\n",
              (int)(ni*nj*nk*nl*nm*nn));
    free((void *)array[0][0][0][0]);
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  if((array[0][0][0][0][0][0] = (void *)calloc(ni*nj*nk*nl*nm*nn*no,size)) == NULL){
    fprintf(stderr,"[calloc_7d] failed to alloc. memory (%d X %d X %d X %d X %d X %d X %d of size %d)\n",
	    (int)ni,(int)nj,(int)nk,(int)nl,(int)nm,(int)nn,(int)no,(int)size);
    free((void *)array[0][0][0][0][0]);
    free((void *)array[0][0][0][0]);
    free((void *)array[0][0][0]);
    free((void *)array[0][0]);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<ni; i++){
    array[i] = (void ******)((unsigned char *)array[i-1] + nj*sizeof(void*****));
  }

  for(j=1; j<nj; j++){
    array[0][j] = (void *****)((unsigned char *)array[0][j-1] + nk*sizeof(void****));
  }

  for(i=1; i<ni; i++){
    array[i][0] = (void *****)((unsigned char *)array[i-1][0] + nj*nk*sizeof(void****));
    for(j=1; j<nj; j++){
      array[i][j] = (void *****)((unsigned char *)array[i][j-1] + nk*sizeof(void****));
    }
  }

  for(k=1; k<nk; k++){
    array[0][0][k] = (void ****)((unsigned char *)array[0][0][k-1] + nl*sizeof(void***));
  }

  for(j=1; j<nj; j++){
    array[0][j][0] = (void ****)((unsigned char *)array[0][j-1][0] + nk*nl*sizeof(void***));
    for(k=1; k<nk; k++){
      array[0][j][k] = (void ****)((unsigned char *)array[0][j][k-1] + nl*sizeof(void***));
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0] = (void ****)((unsigned char *)array[i-1][0][0] + nj*nk*nl*sizeof(void***));
    
    for(k=1; k<nk; k++)
      array[i][0][k] = (void ****)((unsigned char *)array[i][0][k-1] + nl*sizeof(void***));
    
    for(j=1; j<nj; j++){
      array[i][j][0] = (void ****)((unsigned char *)array[i][j-1][0] + nk*nl*sizeof(void***));

      for(k=1; k<nk; k++)
	array[i][j][k] = (void ****)((unsigned char *)array[i][j][k-1] + nl*sizeof(void***));      
    }
  }


  for(l=1; l<nl; l++)
    array[0][0][0][l] = (void ***)((unsigned char *)array[0][0][0][l-1] + nm*sizeof(void**));
  

  for(k=1; k<nk; k++){
    array[0][0][k][0] = (void ***)((unsigned char *)array[0][0][k-1][0] + nl*nm*sizeof(void**));

    for(l=1; l<nl; l++)
      array[0][0][k][l] = (void ***)((unsigned char *)array[0][0][k][l-1] + nm*sizeof(void**));    
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0] = (void ***)((unsigned char *)array[0][j-1][0][0] + nk*nl*nm*sizeof(void**));

    for(l=1; l<nl; l++)
      array[0][j][0][l] = (void ***)((unsigned char *)array[0][j][0][l-1] + nm*sizeof(void**));
    
    for(k=1; k<nk; k++){
      array[0][j][k][0] = (void ***)((unsigned char *)array[0][j][k-1][0] + nl*nm*sizeof(void**));

      for(l=1; l<nl; l++)
	array[0][j][k][l] = (void ***)((unsigned char *)array[0][j][k][l-1] + nm*sizeof(void**));      
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0] = (void ***)((unsigned char *)array[i-1][0][0][0] + nj*nk*nl*nm*sizeof(void**));

    for(l=1; l<nl; l++)
      array[i][0][0][l] = (void ***)((unsigned char *)array[i][0][0][l-1] + nm*sizeof(void**));

    for(k=1; k<nk; k++){
      array[i][0][k][0] = (void ***)((unsigned char *)array[i][0][k-1][0] + nl*nm*sizeof(void**));
      
      for(l=1; l<nl; l++)
	array[i][0][k][l] = (void ***)((unsigned char *)array[i][0][k][l-1] + nm*sizeof(void**));   
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0] = (void ***)((unsigned char *)array[i][j-1][0][0] + nk*nl*nm*sizeof(void**));

      for(l=1; l<nl; l++)
	array[i][j][0][l] = (void ***)((unsigned char *)array[i][j][0][l-1] + nm*sizeof(void**));

      for(k=1; k<nk; k++){
	array[i][j][k][0] = (void ***)((unsigned char *)array[i][j][k-1][0] + nl*nm*sizeof(void**));
	
	for(l=1; l<nl; l++)
	  array[i][j][k][l] = (void ***)((unsigned char *)array[i][j][k][l-1] + nm*sizeof(void**));
      }
    }
  }

  for(m=1; m<nm; m++)
    array[0][0][0][0][m] = (void **)((unsigned char *)array[0][0][0][0][m-1] + nn*sizeof(void*));


  for(l=1; l<nl; l++) {
    array[0][0][0][l][0] = (void **)((unsigned char *)array[0][0][0][l-1][0] + nm*nn*sizeof(void*));

    for(m=1; m<nm; m++)
      array[0][0][0][l][m] = (void **)((unsigned char *)array[0][0][0][l][m-1] + nn*sizeof(void*));
  }

  for(k=1; k<nk; k++){
    array[0][0][k][0][0] = (void **)((unsigned char *)array[0][0][k-1][0][0] + nl*nm*nn*sizeof(void*));

    for(m=1; m<nm; m++)
      array[0][0][k][0][m] = (void **)((unsigned char *)array[0][0][k][0][m-1] + nn*sizeof(void*));

    for(l=1; l<nl; l++) {
      array[0][0][k][l][0] = (void **)((unsigned char *)array[0][0][k][l-1][0] + nm*nn*sizeof(void*));   

      for(m=1; m<nm; m++)
	array[0][0][k][l][m] = (void **)((unsigned char *)array[0][0][k][l][m-1] + nn*sizeof(void*));  	
    } 
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0][0] = (void **)((unsigned char *)array[0][j-1][0][0][0] + nk*nl*nm*nn*sizeof(void*));

    for(m=1; m<nm; m++)
      array[0][j][0][0][m] = (void **)((unsigned char *)array[0][j][0][0][m-1] + nn*sizeof(void*));

    for(l=1; l<nl; l++) {
      array[0][j][0][l][0] = (void **)((unsigned char *)array[0][j][0][l-1][0] + nm*nn*sizeof(void*));
    
      for(m=1; m<nm; m++)
	array[0][j][0][l][m] = (void **)((unsigned char *)array[0][j][0][l][m-1] + nn*sizeof(void*));
    }

    for(k=1; k<nk; k++){
      array[0][j][k][0][0] = (void **)((unsigned char *)array[0][j][k-1][0][0] + nl*nm*nn*sizeof(void*));

      for(m=1; m<nm; m++)
	array[0][j][k][0][m] = (void **)((unsigned char *)array[0][j][k][0][m-1] + nn*sizeof(void*));

      for(l=1; l<nl; l++) {
	array[0][j][k][l][0] = (void **)((unsigned char *)array[0][j][k][l-1][0] + nm*nn*sizeof(void*));    

	for(m=1; m<nm; m++)
	  array[0][j][k][l][m] = (void **)((unsigned char *)array[0][j][k][l][m-1] + nn*sizeof(void*));  
      } 
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0][0] = (void **)((unsigned char *)array[i-1][0][0][0][0] + nj*nk*nl*nm*nn*sizeof(void*));

    for(m=1; m<nm; m++)
      array[i][0][0][0][m] = (void **)((unsigned char *)array[i][0][0][0][m-1] + nn*sizeof(void*));

    for(l=1; l<nl; l++) {
      array[i][0][0][l][0] = (void **)((unsigned char *)array[i][0][0][l-1][0] + nm*nn*sizeof(void*));

      for(m=1; m<nm; m++)
	array[i][0][0][l][m] = (void **)((unsigned char *)array[i][0][0][l][m-1] + nn*sizeof(void*));
    }

    for(k=1; k<nk; k++){
      array[i][0][k][0][0] = (void **)((unsigned char *)array[i][0][k-1][0][0] + nl*nm*nn*sizeof(void*));

      for(m=1; m<nm; m++)
	array[i][0][k][0][m] = (void **)((unsigned char *)array[i][0][k][0][m-1] + nn*sizeof(void*));

      for(l=1; l<nl; l++) {
	array[i][0][k][l][0] = (void **)((unsigned char *)array[i][0][k][l-1][0] + nm*nn*sizeof(void*));    

	for(m=1; m<nm; m++)
	  array[i][0][k][l][m] = (void **)((unsigned char *)array[i][0][k][l][m-1] + nn*sizeof(void*));  
      } 
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0][0] = (void **)((unsigned char *)array[i][j-1][0][0][0] + nk*nl*nm*nn*sizeof(void*));

      for(m=1; m<nm; m++)
	array[i][j][0][0][m] = (void **)((unsigned char *)array[i][j][0][0][m-1] + nn*sizeof(void*));
      
      for(l=1; l<nl; l++) {
	array[i][j][0][l][0] = (void **)((unsigned char *)array[i][j][0][l-1][0] + nm*nn*sizeof(void*));
    
	for(m=1; m<nm; m++)
	  array[i][j][0][l][m] = (void **)((unsigned char *)array[i][j][0][l][m-1] + nn*sizeof(void*));
      }

      for(k=1; k<nk; k++){
	array[i][j][k][0][0] = (void **)((unsigned char *)array[i][j][k-1][0][0] + nl*nm*nn*sizeof(void*));

	for(m=1; m<nm; m++)
	  array[i][j][k][0][m] = (void **)((unsigned char *)array[i][j][k][0][m-1] + nn*sizeof(void*));
	
	for(l=1; l<nl; l++) {
	  array[i][j][k][l][0] = (void **)((unsigned char *)array[i][j][k][l-1][0] + nm*nn*sizeof(void*));    
	  
	  for(m=1; m<nm; m++)
	    array[i][j][k][l][m] = (void **)((unsigned char *)array[i][j][k][l][m-1] + nn*sizeof(void*));  
	} 
      }
    }
  }

  for(n=1; n<nn; n++)
    array[0][0][0][0][0][n] = (void *)((unsigned char *)array[0][0][0][0][0][n-1] + no*size);

  for(m=1; m<nm; m++) {
    array[0][0][0][0][m][0] = (void *)((unsigned char *)array[0][0][0][0][m-1][0] + nn*no*size);

    for(n=1; n<nn; n++)
      array[0][0][0][0][m][n] = (void *)((unsigned char *)array[0][0][0][0][m][n-1] + no*size);
  }

  for(l=1; l<nl; l++) {
    array[0][0][0][l][0][0] = (void *)((unsigned char *)array[0][0][0][l-1][0][0] + nm*nn*no*size);

    for(n=1; n<nn; n++)
      array[0][0][0][l][0][n] = (void *)((unsigned char *)array[0][0][0][l][0][n-1] + no*size);

    for(m=1; m<nm; m++) {
      array[0][0][0][l][m][0] = (void *)((unsigned char *)array[0][0][0][l][m-1][0] + nn*no*size);

      for(n=1; n<nn; n++)
	array[0][0][0][l][m][n] = (void *)((unsigned char *)array[0][0][0][l][m][n-1] + no*size);
    }
  }

  for(k=1; k<nk; k++){
    array[0][0][k][0][0][0] = (void *)((unsigned char *)array[0][0][k-1][0][0][0] + nl*nm*nn*no*size);

    for(n=1; n<nn; n++)
      array[0][0][k][0][0][n] = (void *)((unsigned char *)array[0][0][k][0][0][n-1] + no*size);

    for(m=1; m<nm; m++) {
      array[0][0][k][0][m][0] = (void *)((unsigned char *)array[0][0][k][0][m-1][0] + nn*no*size);

      for(n=1; n<nn; n++)
	array[0][0][k][0][m][n] = (void *)((unsigned char *)array[0][0][k][0][m][n-1] + no*size);
    }

    for(l=1; l<nl; l++) {
      array[0][0][k][l][0][0] = (void *)((unsigned char *)array[0][0][k][l-1][0][0] + nm*nn*no*size);   

      for(n=1; n<nn; n++)
	array[0][0][k][l][0][n] = (void *)((unsigned char *)array[0][0][k][l][0][n-1] + no*size);

      for(m=1; m<nm; m++) {
	array[0][0][k][l][m][0] = (void *)((unsigned char *)array[0][0][k][l][m-1][0] + nn*no*size);  

	for(n=1; n<nn; n++)
	  array[0][0][k][l][m][n] = (void *)((unsigned char *)array[0][0][k][l][m][n-1] + no*size);
      }	
    } 
  }

  for(j=1; j<nj; j++){
    array[0][j][0][0][0][0] = (void *)((unsigned char *)array[0][j-1][0][0][0][0] + nk*nl*nm*nn*no*size);

    for(n=1; n<nn; n++)
      array[0][j][0][0][0][n] = (void *)((unsigned char *)array[0][j][0][0][0][n-1] + no*size);

    for(m=1; m<nm; m++) {
      array[0][j][0][0][m][0] = (void *)((unsigned char *)array[0][j][0][0][m-1][0] + nn*no*size);

      for(n=1; n<nn; n++)
	array[0][j][0][0][m][n] = (void *)((unsigned char *)array[0][j][0][0][m][n-1] + no*size);
    }

    for(l=1; l<nl; l++) {
      array[0][j][0][l][0][0] = (void *)((unsigned char *)array[0][j][0][l-1][0][0] + nm*nn*no*size);

      for(n=1; n<nn; n++)
	array[0][j][0][l][0][n] = (void *)((unsigned char *)array[0][j][0][l][0][n-1] + no*size);    

      for(m=1; m<nm; m++) {
	array[0][j][0][l][m][0] = (void *)((unsigned char *)array[0][j][0][l][m-1][0] + nn*no*size);

	for(n=1; n<nn; n++)
	  array[0][j][0][l][m][n] = (void *)((unsigned char *)array[0][j][0][l][m][n-1] + no*size); 
      }
    }

    for(k=1; k<nk; k++) {
      array[0][j][k][0][0][0] = (void *)((unsigned char *)array[0][j][k-1][0][0][0] + nl*nm*nn*no*size);

      for(n=1; n<nn; n++)
	array[0][j][k][0][0][n] = (void *)((unsigned char *)array[0][j][k][0][0][n-1] + no*size);

      for(m=1; m<nm; m++) {
	array[0][j][k][0][m][0] = (void *)((unsigned char *)array[0][j][k][0][m-1][0] + nn*no*size);

	for(n=1; n<nn; n++)
	  array[0][j][k][0][m][n] = (void *)((unsigned char *)array[0][j][k][0][m][n-1] + no*size);
      }

      for(l=1; l<nl; l++) {
	array[0][j][k][l][0][0] = (void *)((unsigned char *)array[0][j][k][l-1][0][0] + nm*nn*no*size);    

	for(n=1; n<nn; n++)
	  array[0][j][k][l][0][n] = (void *)((unsigned char *)array[0][j][k][l][0][n-1] + no*size);

	for(m=1; m<nm; m++) {
	  array[0][j][k][l][m][0] = (void *)((unsigned char *)array[0][j][k][l][m-1][0] + nn*no*size);

	  for(n=1; n<nn; n++)
	    array[0][j][k][l][m][n] = (void *)((unsigned char *)array[0][j][k][l][m][n-1] + no*size);
	}  
      } 
    }
  }

  for(i=1; i<ni; i++){
    array[i][0][0][0][0][0] = (void *)((unsigned char *)array[i-1][0][0][0][0][0] + nj*nk*nl*nm*nn*no*size);

    for(n=1; n<nn; n++)
      array[i][0][0][0][0][n] = (void *)((unsigned char *)array[i][0][0][0][0][n-1] + no*size);

    for(m=1; m<nm; m++) {
      array[i][0][0][0][m][0] = (void *)((unsigned char *)array[i][0][0][0][m-1][0] + nn*no*size);

      for(n=1; n<nn; n++)
	array[i][0][0][0][m][n] = (void *)((unsigned char *)array[i][0][0][0][m][n-1] + no*size);
    }

    for(l=1; l<nl; l++) {
      array[i][0][0][l][0][0] = (void *)((unsigned char *)array[i][0][0][l-1][0][0] + nm*nn*no*size);

      for(n=1; n<nn; n++)
	array[i][0][0][l][0][n] = (void *)((unsigned char *)array[i][0][0][l][0][n-1] + no*size);

      for(m=1; m<nm; m++) {
	array[i][0][0][l][m][0] = (void *)((unsigned char *)array[i][0][0][l][m-1][0] + nn*no*size);

	for(n=1; n<nn; n++)
	  array[i][0][0][l][m][n] = (void *)((unsigned char *)array[i][0][0][l][m][n-1] + no*size);
      }
    }

    for(k=1; k<nk; k++){
      array[i][0][k][0][0][0] = (void *)((unsigned char *)array[i][0][k-1][0][0][0] + nl*nm*nn*no*size);

      for(n=1; n<nn; n++)
	array[i][0][k][0][0][n] = (void *)((unsigned char *)array[i][0][k][0][0][n-1] + no*size);

      for(m=1; m<nm; m++) {
	array[i][0][k][0][m][0] = (void *)((unsigned char *)array[i][0][k][0][m-1][0] + nn*no*size);

	for(n=1; n<nn; n++)
	  array[i][0][k][0][m][n] = (void *)((unsigned char *)array[i][0][k][0][m][n-1] + no*size);
      }

      for(l=1; l<nl; l++) {
	array[i][0][k][l][0][0] = (void *)((unsigned char *)array[i][0][k][l-1][0][0] + nm*nn*no*size);    

	for(n=1; n<nn; n++)
	  array[i][0][k][l][0][n] = (void *)((unsigned char *)array[i][0][k][l][0][n-1] + no*size);

	for(m=1; m<nm; m++) {
	  array[i][0][k][l][m][0] = (void *)((unsigned char *)array[i][0][k][l][m-1][0] + nn*no*size);  

	  for(n=1; n<nn; n++)
	    array[i][0][k][l][m][n] = (void *)((unsigned char *)array[i][0][k][l][m][n-1] + no*size);
	}
      } 
    }

    for(j=1; j<nj; j++){
      array[i][j][0][0][0][0] = (void *)((unsigned char *)array[i][j-1][0][0][0][0] + nk*nl*nm*nn*no*size);

      for(n=1; n<nn; n++)
	array[i][j][0][0][0][n] = (void *)((unsigned char *)array[i][j][0][0][0][n-1] + no*size);

      for(m=1; m<nm; m++) {
	array[i][j][0][0][m][0] = (void *)((unsigned char *)array[i][j][0][0][m-1][0] + nn*no*size);

	for(n=1; n<nn; n++)
	  array[i][j][0][0][m][n] = (void *)((unsigned char *)array[i][j][0][0][m][n-1] + no*size);
      }

      for(l=1; l<nl; l++) {
	array[i][j][0][l][0][0] = (void *)((unsigned char *)array[i][j][0][l-1][0][0] + nm*nn*no*size);
    
	for(n=1; n<nn; n++)
	  array[i][j][0][l][0][n] = (void *)((unsigned char *)array[i][j][0][l][0][n-1] + no*size);

	for(m=1; m<nm; m++) {
	  array[i][j][0][l][m][0] = (void *)((unsigned char *)array[i][j][0][l][m-1][0] + nn*no*size);

	  for(n=1; n<nn; n++)
	    array[i][j][0][l][m][n] = (void *)((unsigned char *)array[i][j][0][l][m][n-1] + no*size);
	}
      }

      for(k=1; k<nk; k++){
	array[i][j][k][0][0][0] = (void *)((unsigned char *)array[i][j][k-1][0][0][0] + nl*nm*nn*no*size);

	for(n=1; n<nn; n++)
	  array[i][j][k][0][0][n] = (void *)((unsigned char *)array[i][j][k][0][0][n-1] + no*size);

	for(m=1; m<nm; m++) {
	  array[i][j][k][0][m][0] = (void *)((unsigned char *)array[i][j][k][0][m-1][0] + nn*no*size);

	  for(n=1; n<nn; n++)
	    array[i][j][k][0][m][n] = (void *)((unsigned char *)array[i][j][k][0][m][n-1] + no*size);	
	}

	for(l=1; l<nl; l++) {
	  array[i][j][k][l][0][0] = (void *)((unsigned char *)array[i][j][k][l-1][0][0] + nm*nn*no*size);    
	  
	  for(n=1; n<nn; n++)
	    array[i][j][k][l][0][n] = (void *)((unsigned char *)array[i][j][k][l][0][n-1] + no*size);

	  for(m=1; m<nm; m++) {
	    array[i][j][k][l][m][0] = (void *)((unsigned char *)array[i][j][k][l][m-1][0] + nn*no*size);  

	    for(n=1; n<nn; n++)
	      array[i][j][k][l][m][n] = (void *)((unsigned char *)array[i][j][k][l][m][n-1] + no*size);  
	  }
	} 
      }
    }
  }



  return array;
}

/*----------------------------------------------------------------------------*/
/* free_1d_array: free memory used by 1D array  */

void free_1d_array(void *array)
{
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_2d_array: free memory used by 2D array  */

void free_2d_array(void *array)
{
  void **ta = (void **)array;

  free(ta[0]);
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_3d_array: free memory used by 3D array  */

void free_3d_array(void *array)
{
  void ***ta = (void ***)array;

  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_4d_array: free memory used by 3D array  */

void free_4d_array(void *array)
{
  void ****ta = (void ****)array;

  free(ta[0][0][0]);
  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_5d_array: free memory used by 3D array  */

void free_5d_array(void *array)
{
  void *****ta = (void *****)array;

  free(ta[0][0][0][0]);
  free(ta[0][0][0]);
  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_6d_array: free memory used by 3D array  */

void free_6d_array(void *array)
{
  void ******ta = (void ******)array;

  free(ta[0][0][0][0][0]);
  free(ta[0][0][0][0]);
  free(ta[0][0][0]);
  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

/*----------------------------------------------------------------------------*/
/* free_7d_array: free memory used by 3D array  */

void free_7d_array(void *array)
{
  void *******ta = (void *******)array;

  free(ta[0][0][0][0][0][0]);
  free(ta[0][0][0][0][0]);
  free(ta[0][0][0][0]);
  free(ta[0][0][0]);
  free(ta[0][0]);
  free(ta[0]);
  free(array);
}
