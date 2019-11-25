#ifndef NEURALNET_H
#define NEURALNET_H
#include "arrays.h"

typedef struct neuralnet{
  array3d_int *A;// Adjacency matrix has 1 at pos i,j if links between from  j to i. 0 if not
  array3d_double *W; // Weights matrix wij=weight form j to i
  array3d_double *a;      // Activation
  array3d_double *table_act;      // Activation table
  
  int MAX_NEURON;
  int MAX_LINKS;
  int Nneuron; //*Total number of neurons Input+output+hidden
  int Nhidden;
  int Ninput;
  int Noutput;
  int Nlinks;
  int t;
  int cur,old;
}neuralnet;

neuralnet* neuralnet_init(int Ninput,int Noutput,int MAX_NEURON,int MAX_LINKS);
neuralnet* neuralnet_full_init(int Ninput,int Noutput,int MAX_NEURON );

void neuralnet_replace(neuralnet *destination ,neuralnet *source );
void full_neuralnet_replace(neuralnet *destination ,neuralnet *source );
neuralnet * neuralnet_copy(neuralnet *nn); 

long double sigmoid(long double x);
void advance_state(neuralnet *nn);
void compute_act(neuralnet *nn);
void neuralnet_lin_computation(neuralnet *nn);
void run_tanh(neuralnet *nn);
void pass_int_input(neuralnet *nn,array3d_int * input);
void pass_float_input(neuralnet *nn,float * input);
void showact(neuralnet *nn);
long double * get_output(neuralnet *nn);
long double * forward_pass(neuralnet *nn,array3d_int* gamestate);
void float_forward_pass(neuralnet *nn, float* gamestate,float *output);
void neuralnet_reset_act(neuralnet *nn);
void full_mutate_weights(neuralnet *nn,double sigma);
void neuralnet_mutate_table(neuralnet *nn,double sigma);

void mutate_weights(neuralnet *nn,double sigma);
long double randn (double mu, double sigma);

void  float_get_output(neuralnet *nn,float * output);

void neuralnet_write(neuralnet * nn);
void full_neuralnet_write(neuralnet * nn);

void neuralnet_read(neuralnet * nn);

void full_neuralnet_free(neuralnet* nn );
void neuralnet_free(neuralnet* nn);

void full_neuralnet_read(neuralnet * nn);

void neuralnet_write2(neuralnet * nn,FILE* fW,FILE* fA,FILE* fa,FILE * ftable);
void neuralnet_read2(neuralnet * nn,FILE* fW,FILE* fA,FILE* fa,FILE * ftable);

void full_neuralnet_read2(neuralnet * nn,FILE* fW,FILE* fa,FILE * ftable);

void full_neuralnet_write2(neuralnet * nn,FILE* fW,FILE* fa,FILE * ftable);

#endif
