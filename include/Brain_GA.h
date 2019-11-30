#ifndef Brain_GA_H
#define Brain_GA_H

#include "arrays.h"
#include "brain.h"
#include "GA.h"

typedef struct Brain_GA{
  int n; // Pop size
  int Ncluster,Size_cluster,Ncluster_links,Ninput,Noutput;

  float *fit_array; // fitness of chrom in pop
  float *copy_fit; // copy o f fit_array for selection
  float *sigma; // Variance for ES for each chrom in pop 
  brain **pop;// population of chromosones
  brain **copy_pop; // used to copy the pop for selection
  int Nimproved_fit; // number of chrom whose fitness improved
  double c;   // sigma evolution parameter
  
}Brain_GA;
Brain_GA * Brain_GA_init(int n,int Ninput,int Noutput,int Ncluster,int Size_cluster, int Ncluster_links);
Brain_GA * Brain_GA_graph_init(int n,int Ninput,int Noutput,int nx,int ny,int Size_cluster, int Ncluster_links);


int Brain_GA_max_fit(Brain_GA *ga);

int Brain_GA_n_best(Brain_GA *ga, int n);
void Brain_GA_tournament_selection(Brain_GA* bga);
void Brain_GA_mutate_sigma(Brain_GA* bga);
void Brain_GA_mutate_weights(Brain_GA *bga,float pm);
void Brain_GA_free(Brain_GA* bga);
void Brain_GA_mutate_table(Brain_GA *bga,float pm);
  
void Brain_GA_out_fit(FILE * file, Brain_GA *bga);
void Brain_GA_out_sig(FILE * file, Brain_GA *bga);

void Brain_GA_write(Brain_GA *bga);

void Brain_GA_out_fit2(FILE * file, Brain_GA *bga);
void Brain_GA_out_sig2(FILE * file, Brain_GA *bga);
void Brain_GA_read_fit2(FILE * file, Brain_GA *bga);
void Brain_GA_read_sig2(FILE * file, Brain_GA *bga);
void Brain_GA_read(Brain_GA *bga);
#endif
