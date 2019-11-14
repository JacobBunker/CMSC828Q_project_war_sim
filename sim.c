#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "neuralnet.h"
#include "GA.h"
#include "rigid_body_sim.h"


int main( int argc, const char* argv[] )
{
  srand(time(NULL));
  int popsize=40;
  int Ninput=9,Noutput=2,MAX_NEURON=50,MAX_LINKS=2000;
  int Ngame=10;
  int  Learning_time=20;
  GA *ga=GA_init(popsize,Ninput,Noutput,MAX_NEURON,MAX_LINKS);
    FILE * ffit,* fsig;
  ffit=fopen("fit.txt","w");
  fsig=fopen("sig.txt","w");
  
  void race_learning(GA * ga,int Learning_time){
    for(int i=0;i<Learning_time;++i){
      out_fit(ffit,ga);
      out_sig(fsig,ga);
      printf("time: %d\n ",i);
      race(ga,Ngame);
      mutate_sigma(ga);
      tournament_selection(ga);
      GA_mutate_weights(ga,.8);
    }
  }
  printf("Training done \n\n");
  race_learning(ga,Learning_time);
  neuralnet_write(ga->pop[max_fit(ga)]);

  
  int sizeA=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(int);
  int sizeW=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(long double);
  int sizea=ga->MAX_NEURON*2*sizeof(long double);
  neuralnet **players=malloc(4*(sizeof(neuralnet)+sizeA+sizea+sizeW)) ;
  for(int i=0;i<4;++i){
    players[i]=neuralnet_init(Ninput,Noutput,MAX_NEURON,MAX_LINKS);
    neuralnet_replace(players[i],ga->pop[max_fit(ga)]);
  }
  printf("here\n");
  RunRigidBodySimulation(players,1);
  
  for(int i=0;i<4;++i){
    printf("Player %d\n",i);
    array3d_double_show(players[i]->W);
    printf("\n");
  }

  for(int i=0;i<4;++i){

    printf("Player %d\n",i);
    array3d_double_show(players[i]->a);
    printf("\n");

  }

    
  for(int i=0;i<4;++i){

    printf("Player %d\n",i);
    array3d_int_show(players[i]->A);
    printf("\n");

  }
  for(int i=0;i<4;++i){
      neuralnet_free(players[i]);
  }
  free(players);
    
  out_fit(ffit,ga);
  out_sig(fsig,ga);
  GA_free(ga);
  fclose(ffit);
  fclose(fsig);


}
