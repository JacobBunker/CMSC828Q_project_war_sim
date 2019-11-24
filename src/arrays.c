#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "arrays.h"
#define DEBUG_IND 0

array3d_int* array3d_int_init(int dim1,int dim2,int dim3){
  array3d_int *arr=malloc(sizeof(int)*(dim1*dim2*dim3)+sizeof(array3d_int));

  arr->dim1=dim1;
  arr->dim2=dim2;
  arr->dim3=dim3;
  arr->array= (int *)malloc(sizeof(int)*dim1*dim2*dim3);
  assert(arr->array!=NULL);
  int ind;
  for(int k=0;k<arr->dim3;++k){
    for(int i=0;i<arr->dim1;++i){
      for(int j=0;j<arr->dim2;++j){
	ind=array3d_int_index(arr,i,j,k);
	arr->array[ind]=0;
      }
    }
  }
  return arr;
}

int array3d_int_index(array3d_int *arr,int i,int j,int k){
  if(DEBUG_IND){
    printf("INT i %d dim1 %d j %d dim2 %d k %d dim3 %d \n",i,arr->dim1,j,arr->dim2,k,arr->dim3);}
  assert((i<arr->dim1) && (j<arr->dim2) && (k<arr->dim3));
  return i +j*arr->dim1 +k*arr->dim2*arr->dim1;
}

void array3d_int_show(array3d_int *arr){
  int ind;
  for(int k=0;k<arr->dim3;++k){
    for(int i=0;i<arr->dim1;++i){
      for(int j=0;j<arr->dim2;++j){
	ind=array3d_int_index(arr,i,j,k);
	printf(" %d ",arr->array[ind]);
      }
      printf("\n");
    }
  }
}


 array3d_double* array3d_double_init(int dim1,int dim2,int dim3){
   array3d_double *arr=malloc(sizeof(long double)*(dim1*dim2*dim3)+sizeof(array3d_double));
  arr->dim1=dim1;
  arr->dim2=dim2;
  arr->dim3=dim3;
  arr->array= (long double *)malloc(sizeof(long double)*dim1*dim2*dim3);
  for(int i=0;i<dim1*dim2*dim3;++i){
    arr->array[i]=0;
  }
  assert(arr->array!=NULL);
  return arr;
}

int array3d_double_index(array3d_double *arr,int i,int j,int k){
    if(DEBUG_IND){
    printf("Double i %d dim1 %d j %d dim2 %d k %d dim3 %d \n",i,arr->dim1,j,arr->dim2,k,arr->dim3);}
  assert(i<arr->dim1 && j<arr->dim2 && k<arr->dim3);
  return i +j*arr->dim1 +k*arr->dim2*arr->dim1;
}

void  array3d_double_show(array3d_double *arr){
  int ind;
  for(int k=0;k<arr->dim3;++k){
    for(int i=0;i<arr->dim1;++i){
      for(int j=0;j<arr->dim2;++j){
	       ind=array3d_double_index(arr,i,j,k);
	       printf(" %.2LF ",arr->array[ind]);
      }
      printf("\n");
    }
  }
}



void index_sort(float *output, int * index,int size){
  int cmp_ind(const void *i, const void *j)
  {
    const int fi = *(const int *) i;
    const int fj = *(const int *) j;
    const float fa = output[fi];
    const float fb = output[fj];
    if (isnan(fa))
      {
        if (isnan(fb))
	  {
            return 0;
	  }
        return 1;
      }
    if (isnan(fb))
      {
        return -1;
      }
    if (fa > fb) return -1;
    if (fa < fb) return 1;

    /* no more comparisons needed */
    return 0;
  }

  qsort(index,size,sizeof(int),cmp_ind);
}

void index_sort_dcr(float *output, int * index,int size){
  // Sorts index by increasing order of output[index]
  int cmp_ind(const void *i, const void *j)
  {
    const int fi = *(const int *) i;
    const int fj = *(const int *) j;
    const float fa = output[fi];
    const float fb = output[fj];
    if (isnan(fa))
      {
        if (isnan(fb))
	  {
            return 0;
	  }
        return -1;
      }
    if (isnan(fb))
      {
        return 1;
      }
    if (fa > fb) return 1;
    if (fa < fb) return -1;

    /* no more comparisons needed */
    return 0;
  }

  qsort(index,size,sizeof(int),cmp_ind);
}


void value_sort(int n, int * tbs){
  // sort tbs in increasing order.
  int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

  qsort(tbs,n,sizeof(int),cmpfunc);
}


array3d_int * array3d_int_copy(array3d_int * orig){
  array3d_int * copy=array3d_int_init(orig->dim1,orig->dim2,orig->dim3);
    for(int i=0;i<copy->dim1*copy->dim2*copy->dim3;++i){
    copy->array[i]=orig->array[i];
  }
  return copy;
}
  
void array3d_int_replace(array3d_int * destination,array3d_int *source){
  assert(destination->dim1==source->dim1 && destination->dim2==source->dim2 && destination->dim3==source->dim3) ;
  int size_dest=destination->dim1*destination->dim2*destination->dim3*sizeof(int);
  memcpy(destination->array,source->array,size_dest);
  /*   forsizede(int i=0;i<orig->dim1*orig->dim2*orig->dim3;++i){ */
  /*   orig->array[i]=source->array[i]; */
  /* } */
}


array3d_double * array3d_double_copy(array3d_double * orig){
  array3d_double * copy=array3d_double_init(orig->dim1,orig->dim2,orig->dim3);
    for(int i=0;i<copy->dim1*copy->dim2*copy->dim3;++i){
    copy->array[i]=orig->array[i];
  }
  return copy;
}

void array3d_double_replace(array3d_double * destination,array3d_double *source){
  assert(destination->dim1==source->dim1 && destination->dim2==source->dim2 && destination->dim3==source->dim3);
    int size_dest=destination->dim1*destination->dim2*destination->dim3*sizeof(long double);
    memcpy(destination->array,source->array,size_dest);
}

int array3d_double_write(FILE * file, array3d_double *arr3d){
  int size_arr=arr3d->dim1*arr3d->dim2*arr3d->dim3;
  
  if(fwrite(arr3d->array, sizeof(long double),size_arr, file) != size_arr){
    printf("File write error.\n");
  }
  return size_arr;
}
void array3d_double_read(FILE * file,array3d_double *arr3d){
  int size=arr3d->dim1*arr3d->dim2*arr3d->dim3;
  if(fread(arr3d->array, sizeof(long double), size, file)!= size) {
    if(feof(file))
      printf("Premature end of file.\n");
    else
      printf("File read error.\n");
  }
}
int array3d_int_write(FILE * file, array3d_int *arr3d){
  int size_arr=arr3d->dim1*arr3d->dim2*arr3d->dim3;
  
  if(fwrite(arr3d->array, sizeof(int),size_arr, file) != size_arr){
    printf("File write error.\n");
  }
  return size_arr;
}
void array3d_int_read(FILE * file,array3d_int *arr3d){
  int size=arr3d->dim1*arr3d->dim2*arr3d->dim3;
  if(fread(arr3d->array, sizeof(int), size, file)!= size) {
    if(feof(file))
      printf("Premature end of file.\n");
    else
      printf("File read error.\n");
  }
}

