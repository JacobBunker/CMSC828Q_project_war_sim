#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>

#include "neuralnet.h"
#include "GA.h"
#include "Brain_GA.h"
#include "brain.h"

int popsize=10,
  Ninput=10,
  Noutput=10,
  nx=10,
  ny=5,
  Size_cluster=30,
  Ncluster_links=100,
  Ngame=4,
  Learning_time=40;
int main(){
  Brain_GA *bga1=Brain_GA_graph_init(popsize,Ninput,Noutput,nx,ny,Size_cluster,Ncluster_links),
    *bga2=Brain_GA_graph_init(popsize,Ninput,Noutput,nx,ny,Size_cluster,Ncluster_links);
  float input[10]={0,1,2,3,4,5,6,7,8,9};
  brain_pass_float_input(bga1->pop[0],input);
  Brain_GA_write(bga1);
  Brain_GA_read(bga2);
 
  Brain_GA_free(bga1);
  Brain_GA_free(bga2);

  
 return 0;
}
 
