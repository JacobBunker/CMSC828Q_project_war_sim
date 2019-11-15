 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rigid_body_sim.h"  
#include "GA.h"
#include "neuralnet.h"

GA * GA_init(int n,int Ninput,int Noutput,int MAX_NEURON,int MAX_LINKS){
  int sizeA=MAX_NEURON*MAX_NEURON*sizeof(int);
  int sizeW=MAX_NEURON*MAX_NEURON*sizeof(long double);
  int sizea=MAX_NEURON*2*sizeof(long double);
  int sizefit=n*sizeof(long double);
  int sizesig=n*sizeof(double);
  
  GA*ga=malloc(2*n*(sizeof(neuralnet)+sizeA+sizeW+sizea)+2*sizefit+sizesig+sizeof(GA));
  if((ga = (GA *) malloc(BUFSIZ)) == NULL) {
    printf("malloc error in GA_init");
    return 0;
  }

  assert(!(n%2));
  ga->pop=malloc(n*(sizeof(neuralnet)+sizeA+sizeW+sizea));
  ga->copy_pop=malloc(n*(sizeof(neuralnet)+sizeA+sizeW+sizea));
  ga->n=n; 
  ga->MAX_NEURON=MAX_NEURON;
  ga->MAX_LINKS=MAX_LINKS;
  ga->Ninput=Ninput;
  ga->Noutput=Noutput;
  ga->fit_array=malloc(sizefit);
  ga->copy_fit=malloc(sizefit);
  ga->sigma=malloc(sizesig);
  ga->c=0.85;
  
  for(int i=0;i<n;++i){
    ga->pop[i]=neuralnet_init(Ninput,Noutput,MAX_NEURON,MAX_LINKS);
    ga->copy_pop[i]=neuralnet_init(Ninput,Noutput,MAX_NEURON,MAX_LINKS);
    ga->fit_array[i]=0;
    ga->sigma[i]=rand()%10 + 0.01;
  }
  return ga;
} 


void permute(int* per,int n){
  // return a random permutation
  //of 1,2,...,n
  int temp,rd,i;
  for(i=0;i<n;++i){
    per[i]=i;
  }
  
  for(i=n-1;i>0;--i){
    rd=rand()%i+1;
    temp=per[i];
    per[i]=per[rd];
    per[rd]=temp;
  } 
}

void mutate_sigma(GA* ga){
  //printf("sw %d\n",ga->Nimproved_fit>ga->n/4.0);
  switch(ga->Nimproved_fit>ga->n/2.0){
  case 1:
    for(int i=0;i<ga->n;++i){
      /* printf("case 1 old sig %f",ga->sigma[i]); */
      ga->sigma[i]*=1/ga->c;
      /* printf(" new sig %f\n",ga->sigma[i]); */
    }
    break;
  case 0:
    for(int i=0;i<ga->n;++i){
      /* printf("case 0 old sig %f",ga->sigma[i]);  */
      ga->sigma[i]*=ga->c;
      /* printf(" new sig %f\n",ga->sigma[i]); */
    }
    break;
  }
    
}

int max_fit(GA *ga){
  int max=0;
  for(int i=1;i<ga->n;++i){
    if(ga->fit_array[i]>ga->fit_array[max])
      {max=i;}
  }
  return max;
}
void tournament_selection(GA* ga){
  /*
    replace ga->pop with new population of same size.
    by doing n times the following:
    chose 2 parents at random, copy one with best fitness in new pop
  */
  int parents[ga->n]; //array keeping indices of chosen parents 
  int par1=0,par2=0;   // random indices
  for(int i=0;i<ga->n;++i){
    parents[i]=0;
    par1=0;par2=0;
    while(par1==par2){
      // chose 2 diffent parents at random
      par1=rand()%ga->n;
      par2=rand()%ga->n;
    }
    switch(ga->fit_array[par1]>ga->fit_array[par2]){
      //keep parent with best fitness
    case 1:
      parents[i]=par1;
      break;
    case 0:
      parents[i]=par2;
      break;
    }
  }
  value_sort(ga->n,parents); // sorts parents indices in increasing order
  //------------------------------------------------------------------
  // copy all chrom selected in copy_pop to be copied back later in pop
  neuralnet_replace(ga->copy_pop[parents[0]],ga->pop[parents[0]]);
  for(int i=1;i<ga->n;++i){
    if(parents[i]!=parents[i-1]){   // Test for unnecessary copying 
      neuralnet_replace(ga->copy_pop[parents[i]],ga->pop[parents[i]]);
    }
  }
  //------------------------------------------------------------------

  for(int i=1;i<ga->n;++i){
    neuralnet_replace(ga->pop[i],ga->copy_pop[parents[i]]); // replace pop 
    ga->fit_array[i]=ga->copy_fit[parents[i]];         // replace fitness
  }
  //  ga->fit_avg*=1/ga->n;	    
}    



void GA_mutate_weights(GA *ga,double pm){
  assert(0<pm && pm<1);
  for(int i=0;i<ga->n;++i){
    if(rand_double()<pm){
      mutate_weights(ga->pop[i],ga->sigma[i]);
    }
  }

}

double rand_double() {
  return rand()/(double)RAND_MAX;
}

void GA_free(GA* ga){
  for(int i=0;i<ga->n;++i){
    neuralnet_free(ga->pop[i]);
    neuralnet_free(ga->copy_pop[i]);
  }
  free(ga->fit_array);
  free(ga->copy_fit);
  free(ga->sigma);
  free(ga->pop);
  free(ga->copy_pop);
  free(ga);
  
}

void race(GA *ga,int Ngame){
  int sizeA=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(int);
  int sizeW=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(long double);
  int sizea=ga->MAX_NEURON*2*sizeof(long double);
  neuralnet **players=malloc(4*(sizeof(neuralnet)+sizeA+sizea+sizeW)) ;
  assert(ga->n%4==0);
  int schedule[ga->n];
  //  long double fit_change[ga->n];
  int outcome;
  ga->Nimproved_fit=0;
  for (int i=0;i<ga->n;++i){
    ga->copy_fit[i]=0;
  }
  for(int i=0;i<Ngame;++i){
    permute(schedule,ga->n);
    for(int j=0;j<ga->n/4;++j){
      for(int ii=0;ii<4;++ii){
	       players[ii]=ga->pop[schedule[2*j+ii]];
      }
      outcome=RunRigidBodySimulation(players,0);

      if(outcome!=4){
        ga->copy_fit[schedule[2*j+outcome]]+=10;
        for(int iii=0;iii<4;++iii){
          ga->copy_fit[schedule[2*j+iii]]-=1;
        }
      }
    }
  }
  for (int i=0;i<ga->n;++i){
    /* printf("old fit : %LF , newfit : %LF\n",ga->fit_array[i],ga->copy_fit[i]); */
    switch(ga->copy_fit[i]>ga->fit_array[i]){
    case 1:
      ++ga->Nimproved_fit;
      break;
    }
  }
  /* printf("Nimproved %d\n",ga->Nimproved_fit); */
  
  memcpy(ga->fit_array,ga->copy_fit,ga->n*sizeof(long double));
  /* for (int i=0;i<ga->n;++i){ */
  /*   /\* printf("Again old fit : %LF , newfit : %LF\n",ga->fit_array[i],ga->copy_fit[i]); *\/ */
  /*   } */
  free(players);
}


void out_fit(FILE * file, GA *ga){ 
  for(int i=0;i<ga->n;++i){
    fprintf(file,"%05.2LF ",ga->fit_array[i]);    
  }
  fprintf(file,"\n");
}

void out_sig(FILE * file, GA *ga){
  for(int i=0;i<ga->n;++i){
    fprintf(file,"%05.2f ",ga->sigma[i]);    
  }
  fprintf(file,"\n");
}
