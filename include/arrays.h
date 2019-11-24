#ifndef ARRAYS_H
#define ARRAYS_H

typedef struct array3d_int{
  int dim1,dim2,dim3;
  int *array;
}array3d_int;


array3d_int * array3d_int_init(int dim1,int dim2,int dim3);
int array3d_int_index(array3d_int *arr,int i,int j,int k);
void array3d_int_show(array3d_int *arr);
void array3d_int_replace(array3d_int * destination,array3d_int *source);


typedef struct array3d_double{
  int dim1,dim2,dim3;
  long double *array;
}array3d_double;


array3d_double * array3d_double_init(int dim1,int dim2,int dim3);
int array3d_double_index(array3d_double *arr,int i,int j,int k);
void array3d_double_show(array3d_double *arr);

void value_sort(int n, int * tbs);
array3d_int * array3d_int_copy(array3d_int * orig);
array3d_double * array3d_double_copy(array3d_double * orig);
void array3d_double_replace(array3d_double * destination,array3d_double *source);

void array3d_int_read(FILE * file,array3d_int *arr3d);
int array3d_int_write(FILE * file,array3d_int *arr3d);
void array3d_double_read(FILE * file,array3d_double *arr3d);
int array3d_double_write(FILE * file,array3d_double *arr3d);
void index_sort_inc(float *array, int * index,int size); 
void index_sort_dcr(float *output, int * index,int size);
#endif
