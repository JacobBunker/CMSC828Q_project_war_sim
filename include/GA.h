#ifndef GA_H
#define GA_H
#include "arrays.h"
#include "neuralnet.h"

typedef struct GA{
  int n; // Pop size
  int MAX_NEURON,MAX_LINKS,Ninput,Noutput;

  float *fit_array; // fitness of chrom in pop
  float *copy_fit; // copy o f fit_array for selection
  float  *sigma; // Variance for ES for each chrom in pop 
  neuralnet **pop;// population of chromosones
  neuralnet **copy_pop; // used to copy the pop for selection
  int Nimproved_fit; // number of chrom whose fitness improved
  double c;   // sigma evolution parameter
  
}GA;
GA * GA_init(int n,int Ninput,int Noutput,int MAX_NEURON,int MAX_LINKS);
void permute(int* per,int n); 
int max_fit(GA *ga);
void tournament_selection(GA* ga);
void mutate_sigma(GA*ga);
void GA_mutate_weights(GA *ga,float pm);
double rand_double();
void GA_free(GA* ga);
void race(GA *ga,int Ngame);
void out_fit(FILE * file, GA *ga);
void out_sig(FILE * file, GA *ga);
int n_best(GA *ga, int n);  
#endif
