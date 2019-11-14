#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "arrays.h"
#include "neuralnet.h"
#include "tictactoe.h"
#include "GA.h"

/*
FIX ES
change fit to double
input treatment
use memcpy
Some neurons may not have weights :s
Implementer maxlinks
nutate weights fix non opti
fix index sort (not general to number of input)
Tournament selection rethink way of doing things.
Change get output/input to pointers
*/


int main(){
   srand(time(NULL));
   int popsize=10;
  int Ninput=9,Noutput=9,MAX_NEURON=24,MAX_LINKS=400;
  int Ngame=20;
  int  Learning_time=10;
  GA *ga=GA_init(popsize,Ninput,Noutput,MAX_NEURON,MAX_LINKS);
  FILE * ffit,* fsig;
  ffit=fopen("fit.txt","w");
  fsig=fopen("sig.txt","w");
  void learning(GA * ga,int Learning_time){
    for(int i=0;i<Learning_time;++i){
      out_fit(ffit,ga);
      out_sig(fsig,ga);
      printf("time: %d\n ",i);
      match(ga,Ngame);
      mutate_sigma(ga);
      tournament_selection(ga);
      GA_mutate_weights(ga,.8);
    }
  }
  
  learning(ga,Learning_time);
  show_game_of_ttt(ga->pop[max_fit(ga)],ga->pop[max_fit(ga)]);
  human_vs_nn(ga->pop[max_fit(ga)]);
    
  out_fit(ffit,ga);
  out_sig(fsig,ga);
  GA_free(ga);
  fclose(ffit);
  fclose(fsig);



  return 0;
}
