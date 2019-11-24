#ifndef Brain_H
#define Brain_H

#include "neuralnet.h"
#include "arrays.h"

typedef struct brain{
  int Ninput,Noutput,Ncluster, Size_cluster,Ncluster_links,Ntotal_links,br_size;
  neuralnet **cluster;
  array3d_int *A;
  array3d_double *W;
}brain;

brain * brain_init (int Ninput,int Noutput,int Ncluster,int Size_cluster,int Ncluster_links);
brain * brain_graph_init(int Ninput,int Noutput,int nx,int ny,int Size_cluster,int Ncluster_links);
void brain_show_act(brain *br);
void brain_show_table(brain *br);
void brain_outer_cluster_compute_lin(brain* br,int i,int j,int weight_ind);
void brain_forward_pass(brain *br,float* input,float* output);
void brain_show_weights(brain *br);
void brain_free(brain *br);



void brain_replace(brain *destination ,brain *source );
void brain_mutate_weights(brain *br,double sigma);
void  brain_float_get_output(brain *br,float * output);

void brain_pass_float_input(brain *br,float * input);
array3d_int * graph_builder(int nx, int ny);

void brain_reset_act(brain * br);
void brain_mutate_table(brain * br, double sigma);






void brain_write(brain * br);
void brain_read(brain *br );

#endif
