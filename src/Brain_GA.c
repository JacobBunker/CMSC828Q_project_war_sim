#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "GA.h"
#include "neuralnet.h"
//#include "brain.h"
#include "Brain_GA.h"

#define elitism 1


Brain_GA * Brain_GA_init(int n,int Ninput,int Noutput,int Ncluster,int Size_cluster,int Ncluster_links){

  int sizefit=n*sizeof(float);
  int sizesig=n*sizeof(float);
  int cluster_sizeA=Size_cluster*Size_cluster*sizeof(int),
    cluster_sizeW=Size_cluster*Size_cluster*sizeof(long double),
    cluster_sizea=Size_cluster*2*sizeof(long double),
    cluster_sizetable=Size_cluster*sizeof(long double),
    sizeA=Ncluster*Ncluster*sizeof(int), 
    sizeW=Ncluster*Ncluster*sizeof(long double); 
  
  int cluster_size=sizeof(neuralnet)+cluster_sizea+cluster_sizeW+cluster_sizeA+cluster_sizetable; 
  int brain_size=sizeof(brain)+cluster_size*Ncluster+sizeA+sizeW;

  Brain_GA * bga=malloc(sizeof(Brain_GA));// +2*sizefit+sizesig+2*n*brain_size);

  bga->pop=malloc(n*brain_size);
  bga->copy_pop=malloc(n*brain_size);
  bga->n=n;
  bga->Ncluster=Ncluster;
  bga->Ncluster_links=Ncluster_links;
  bga->Size_cluster=Size_cluster;
  bga->Ninput=Ninput;
  bga->Noutput=Noutput;
  bga->fit_array=malloc(sizefit);
  bga->copy_fit=malloc(sizefit);
  bga->sigma=malloc(sizesig);
  bga->c=.99;

  for(int i=0; i<bga->n; ++i){ bga->pop[i]=brain_init(bga->Ninput,bga->Noutput,bga->Ncluster,bga->Size_cluster,bga->Ncluster_links);  bga->copy_pop[i]=brain_init(bga->Ninput,bga->Noutput,bga->Ncluster,bga->Size_cluster,bga->Ncluster_links);
    bga->fit_array[i]=0;
    bga->sigma[i]=rand()%10 + 0.01;
  }

  return bga;
}


Brain_GA * Brain_GA_graph_init(int n,int Ninput,int Noutput,int nx,int ny,int Size_cluster,int Ncluster_links){
  int Ncluster=nx*ny;
  int sizefit=n*sizeof(float);
  int sizesig=n*sizeof(float); 
  int cluster_sizeA=Size_cluster*Size_cluster*sizeof(int),
    cluster_sizeW=Size_cluster*Size_cluster*sizeof(long double),
    cluster_sizea=Size_cluster*2*sizeof(long double),
    cluster_sizetable=Size_cluster*sizeof(long double),
    sizeA=Ncluster*Ncluster*sizeof(int), 
    sizeW=Ncluster*Ncluster*sizeof(long double); 
  
  int cluster_size=sizeof(neuralnet)+cluster_sizea+cluster_sizeW+cluster_sizeA+cluster_sizetable; 
  int brain_size=sizeof(brain)+cluster_size*Ncluster+sizeA+sizeW;

  Brain_GA * bga=malloc(sizeof(Brain_GA)+2*sizefit+sizesig+2*n*brain_size);

  bga->pop=malloc(n*brain_size);
  bga->copy_pop=malloc(n*brain_size);
  bga->n=n;
  bga->Ncluster=Ncluster;
  bga->Ncluster_links=Ncluster_links;
  bga->Size_cluster=Size_cluster;
  bga->Ninput=Ninput;
  bga->Noutput=Noutput;
  bga->fit_array=malloc(sizefit);
  bga->copy_fit=malloc(sizefit);
  bga->sigma=malloc(sizesig);
  bga->c=.99;

  for(int i=0; i<bga->n; ++i){ bga->pop[i]=brain_graph_init(bga->Ninput,bga->Noutput,nx,ny,bga->Size_cluster,bga->Ncluster_links);  bga->copy_pop[i]=brain_graph_init(bga->Ninput,bga->Noutput,nx,ny,bga->Size_cluster,bga->Ncluster_links);
    bga->fit_array[i]=0;
    bga->sigma[i]=rand()%10 + 0.01;
  }

  return bga;
}






int Brain_GA_max_fit(Brain_GA *bga){
  int max=0;
  for(int i=1;i<bga->n;++i){
    if(bga->fit_array[i]>bga->fit_array[max])
      {max=i;}
  }
  return max;
}

void Brain_GA_tournament_selection(Brain_GA* bga){
  /*
    replace ga->pop with new population of same size.
    by doing n times the following:
    chose 2 parents at random, copy one with best fitness in new pop
  */
  int parents[bga->n]; //array keeping indices of chosen parents
  int par1=0,par2=0;   // random indices
  
  for(int i=0;i<bga->n;++i){
    parents[i]=i;
  }
 
  index_sort_dcr(bga->fit_array,parents,bga->n);
  for(int i=elitism;i<bga->n;++i){
    parents[i]=0;
    par1=0;par2=0;
    while(par1==par2){
      // chose 2 diffent parents at random
      par1=rand()%bga->n;
      par2=rand()%bga->n;
    }
    switch(bga->fit_array[par1] > bga->fit_array[par2]){
      //keep parent with best fitness
    case 1:
      parents[i]=par1;
      break;
    case 0:
      parents[i]=par2;
      break;
    }
  }
  value_sort(bga->n,parents); // sorts parents indices in increasing order
 
//------------------------------------------------------------------
  // copy all chrom selected in copy_pop to be copied back later in pop
 brain_replace(bga->copy_pop[parents[elitism]],bga->pop[parents[elitism]]);
  for(int i=elitism;i<bga->n;++i){
    if(parents[i]!=parents[i-1]){   // Test for unnecessary copying
      brain_replace(bga->copy_pop[parents[i]],bga->pop[parents[i]]);
    }
  }
}

void Brain_GA_mutate_sigma(Brain_GA* bga){
  //printf("sw %d\n",ga->Nimproved_fit>ga->n/4.0);
  switch(bga->Nimproved_fit>bga->n/2.0){
  case 1:
    for(int i=0;i<bga->n;++i){
      /* printf("case 1 old sig %f",ga->sigma[i]); */
      bga->sigma[i]/=fmax(bga->c,.2);
      /* printf(" new sig %f\n",ga->sigma[i]); */
    }
    break;
  case 0:
    for(int i=0;i<bga->n;++i){
      /* printf("case 0 old sig %f",ga->sigma[i]);  */
      bga->sigma[i]*=fmax(bga->c,.2);
      /* printf(" new sig %f\n",ga->sigma[i]); */
    }
    break;
  }
}    


int Brain_GA_n_best(Brain_GA *ga, int n) {
  int parents[ga->n]; //array keeping indices of chosen parents
 for(int i=0;i<ga->n;++i){
   parents[i]=i;
 }
 index_sort_dcr(ga->fit_array,parents,ga->n);
 return parents[n];
}



void Brain_GA_mutate_weights(Brain_GA *bga,float pm){
  assert(0<pm && pm<1);
  for(int i=0;i<bga->n;++i){
    if(rand_double()<pm){
      brain_mutate_weights(bga->pop[i],bga->sigma[i]);
    }
  }
}

void Brain_GA_mutate_table(Brain_GA *bga,float pm){
  assert(0<pm && pm<1);
  for(int i=0;i<bga->n;++i){
    if(rand_double()<pm){
      brain_mutate_table(bga->pop[i],bga->sigma[i]);
    }
  }
}
void Brain_GA_free(Brain_GA* bga){
    for(int i=0;i<bga->n;++i){
    brain_free(bga->pop[i]);
    brain_free(bga->copy_pop[i]);
  }
  free(bga->fit_array);
  free(bga->copy_fit);
  free(bga->sigma);
  free(bga->pop);
  free(bga->copy_pop);
  free(bga);

}

void Brain_GA_out_fit(FILE * file, Brain_GA *bga){ 
  for(int i=0;i<bga->n;++i){
    fprintf(file,"%05.2f ",bga->fit_array[i]);    
  }
  fprintf(file,"\n");
}

void Brain_GA_out_sig(FILE * file, Brain_GA *bga){
  for(int i=0;i<bga->n;++i){
    fprintf(file,"%05.2f ",bga->sigma[i]);    
  }
  fprintf(file,"\n");
}
