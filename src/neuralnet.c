
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "neuralnet.h"

#define SPREAD 1

neuralnet* neuralnet_init(int Ninput,int Noutput,int MAX_NEURON,int MAX_LINKS){
  /*This functions instanciates the struct neuralnet.It defines the numberof inputs and outputs and  initializes the weight matrix (W), adjancy matrix(A) and  inital activations (a) at random.*/
  int indA,indW,Notyet;
  assert(Ninput+Noutput+1<MAX_NEURON && MAX_LINKS<MAX_NEURON*MAX_NEURON+1 &&MAX_LINKS>Noutput); 
  
  int Nhidden=rand()%(MAX_NEURON-1 -Ninput-Noutput)+1; // Nhidden>5.
  /*  const in size= sizeof(neuralnet)+sizeof(long double) *MAX_NEURON*(MAX_NEURON*2+2)+sizeof(int) *MAX_NEURON*MAX_NEURON*2;*/
  /* int sizeA=MAX_NEURON*MAX_NEURON*sizeof(int); */
  /* int   sizeW=MAX_NEURON*MAX_NEURON*sizeof(long double); */
  /* int   sizea=MAX_NEURON*2*sizeof(long double); */
  /* int   sizetable=MAX_NEURON*sizeof(long double); */
    
  neuralnet* nn=malloc(sizeof(neuralnet));//+sizeA+sizeW+sizea+sizetable);
  

  nn->A=array3d_int_init(MAX_NEURON,MAX_NEURON,1);
  nn->W=array3d_double_init(MAX_NEURON,MAX_NEURON,1);
  nn->a=array3d_double_init(MAX_NEURON,2,1);
  nn->table_act=array3d_double_init(MAX_NEURON,1,1);
  
  nn->t=0;
  nn->old=0;
  nn->cur=1;

  nn->MAX_NEURON=MAX_NEURON;
  nn->MAX_LINKS=MAX_LINKS;
  nn->Nhidden=Nhidden;
  nn->Ninput=Ninput;
  nn->Noutput=Noutput;
  nn->Nneuron=Ninput+Noutput+Nhidden;
  nn->Nlinks=0;

  for(int i=0;i<nn->Nneuron;++i){
    nn->table_act->array[i]=(long double) (2*random())/ INT_MAX -1;
  }
  
  int a=0,b=0;
  /* Next loop initialises A and W
     Nlinks are chosen at random between the neurons
     weights are chosen uniformly in [-10,10]
  */

  // Link_control is the max number of links created at init
  // it is random, but has to allow for enough bias links
  
  int Link_control=fmin(nn->MAX_LINKS,rand()%((nn->Nneuron*nn->Nneuron)-nn->Nneuron));
  Link_control=fmax(Link_control,nn->Noutput);
  long double temp;
  
  Notyet=1;
  for(int i=nn->Nneuron-nn->Noutput;i<nn->Nneuron;++i){
    // Loop for bias weights
    indA=array3d_int_index(nn->A,i,nn->Ninput,0);
    nn->A->array[indA]=1;
    indW=array3d_double_index(nn->W,i,nn->Ninput,0);
    temp=(long double) (SPREAD*random())/ INT_MAX -SPREAD/2;
    nn->W->array[indW]=temp;
    ++nn->Nlinks;
  }

  while(nn->Nlinks<Link_control){
    Notyet=1;
    while(Notyet){
      // Loop that initiate weights at random.
      a=0,b=0;
      while(a==b){
	// this loop make sure a!=b so that a neuron is not linked to itself
	a=rand() %nn->Nneuron;
	b=rand() %nn->Nneuron;
      }
      indA=array3d_int_index(nn->A,a,b,0);
      if(nn->A->array[indA]==0){
	nn->A->array[indA]=1;
	indW=array3d_double_index(nn->W,a,b,0);
	temp=(long double) (SPREAD*random())/ INT_MAX -SPREAD/2;
	nn->W->array[indW]=temp;
	++nn->Nlinks;
	Notyet=0;
      }
    }
  }
  
  return nn;
  
}



neuralnet* neuralnet_full_init(int Ninput,int Noutput,int MAX_NEURON){
  /*This functions instanciates the struct neuralnet.It defines the numberof inputs and outputs and  initializes the weight matrix (W), adjancy matrix(A) and  inital activations (a) at random.*/
  int indW;
    assert(Ninput+Noutput+1<MAX_NEURON); 
  
    int Nhidden=MAX_NEURON-Ninput-Noutput;
  neuralnet* nn=malloc(sizeof(neuralnet));
  

  //  nn->A=array3d_int_init(MAX_NEURON,MAX_NEURON,1);
  nn->W=array3d_double_init(MAX_NEURON,MAX_NEURON,1);
  nn->a=array3d_double_init(MAX_NEURON,2,1);
  nn->table_act=array3d_double_init(MAX_NEURON,1,1);
  
  nn->t=0;
  nn->old=0;
  nn->cur=1;

  nn->MAX_NEURON=MAX_NEURON;
  nn->MAX_LINKS=MAX_NEURON*MAX_NEURON-MAX_NEURON;
  nn->Nhidden=Nhidden;
  nn->Ninput=Ninput;
  nn->Noutput=Noutput;
  nn->Nneuron=MAX_NEURON;
  nn->Nlinks=0;


  for(int i=0;i<nn->Nneuron;++i){
    nn->table_act->array[i]=(long double) (2*random())/ INT_MAX -1;
  }
  /* Next loop initialises A and W

     weights are chosen uniformly in [-10,10]*/

 int inda= array3d_double_index(nn->a,Ninput,nn->cur,0);
 nn->a->array[inda]=1;
  for(int i=0;i<nn->Nneuron;++i){
    for(int j=0;j<nn->Nneuron;++j){
      if(i!=j){
	indW=array3d_double_index(nn->W,i,j,0);
	nn->W->array[indW]=(long double) (SPREAD*random())/ INT_MAX -SPREAD/2;
	++nn->Nlinks;
      }
    }
  }

      
  return nn;
  
}




  void advance_state(neuralnet *nn){
    int tmp;                        /* temporary storage */
   int inda;
    ++(nn->t);                            /* increment time */
    tmp=nn->cur; nn->cur=nn->old; nn->old=tmp;     /* swap state location */
    for(int i=0;i<nn->Nneuron;++i){
      /*Sets current activation at 0 before they are computed*/
      /*from the old activation*/
      inda=array3d_double_index(nn->a,i,nn->cur,0);
      nn->a->array[inda]=0;
    }
  } /* end advancestate */


long double sigmoid (long double x){
  return 1/(1+exp(-x));
}

void compute_act(neuralnet *nn){
  int indW,indacur,indaold,inda;
  long double temp=0;
  inda=array3d_double_index(nn->a,nn->Ninput,nn->old,0);
  nn->a->array[inda]=1;
  for(int i=nn->Ninput;i<nn->Nneuron ;++i){
    indacur=array3d_double_index(nn->a,i,nn->cur,0);
    temp=0;
    for(int j=0;j<nn->Nneuron;++j){
      indW=array3d_double_index(nn->W,i,j,0);
      indaold=array3d_double_index(nn->a,j,nn->old,0);
      temp+=(nn->W->array[indW])*(nn->a->array[indaold]);
    }
    
    nn->a->array[indacur]=tanh(0.1*temp);
     
  }
}

void neuralnet_lin_computation(neuralnet *nn){
  int indW,indacur,indaold;
  long double temp=0;
  
  int inda= array3d_double_index(nn->a,nn->Ninput,nn->old,0);
  nn->a->array[inda]=1;// Bias reset
  
  for(int i=nn->Ninput;i<nn->Nneuron ;++i){
    indacur=array3d_double_index(nn->a,i,nn->cur,0);
    temp=0;
    for(int j=0;j<nn->Nneuron;++j){
      indW=array3d_double_index(nn->W,i,j,0);
      indaold=array3d_double_index(nn->a,j,nn->old,0);
      temp+=(nn->W->array[indW])*(nn->a->array[indaold]);
    }
    nn->a->array[indacur]=temp;
  }

}

  void run_tanh(neuralnet *nn){
    int inda;
    for (int i=nn->Ninput;i<nn->Nneuron;++i){
      inda=array3d_double_index(nn->a,i,nn->cur,0);
      nn->a->array[inda]=tanh(.1*(nn->a->array[inda]));
    }
  }

void pass_int_input(neuralnet* nn,array3d_int* input){
  assert(input->dim1==nn->Ninput && input->dim2==1 && input->dim3==1);
  int inda;
  for (int i=0;i<nn->Ninput;++i){
    inda=array3d_double_index(nn->a,i,nn->cur,0);
    nn->a->array[inda]=input->array[i];
  }
}
void pass_float_input(neuralnet* nn, float* input){
  int inda;
  for (int i=0;i<nn->Ninput;++i){
    inda=array3d_double_index(nn->a,i,nn->cur,0);
    nn->a->array[inda]=input[i];
  }
}

void showact(neuralnet *nn){
  int ind;
  for(int i=0;i<nn->a->dim1;++i)
    {
      ind=array3d_double_index(nn->a,i,nn->cur,0);
      printf("%LF \n",nn->a->array[ind]);
    }
      
}


void neuralnet_reset_act(neuralnet *nn){
  int inda,size_of_act;
  size_of_act=sizeof(long double)*nn->MAX_NEURON;
  inda=array3d_double_index(nn->a,0,nn->cur,0);
  memcpy(&(nn->a->array[inda]),&(nn->table_act->array[0]),size_of_act);
}




long double * get_output(neuralnet *nn){
  long double * output;
  int inda;
  output=malloc(sizeof(long double)*nn->Noutput);
  for (int i=0;i<nn->Noutput;++i){
    inda=array3d_double_index(nn->a,nn->Ninput+nn->Nhidden+i,nn->cur,0);
    output[i]=nn->a->array[inda];
  }
  return output;
}

void float_get_output(neuralnet *nn,float * output){

  int inda;
  for (int i=0;i<nn->Noutput;++i){
    inda=array3d_double_index(nn->a,nn->Ninput+nn->Nhidden+i,nn->cur,0);
    output[i]=nn->a->array[inda];
  }
  /* printf("OUT :%f \n",output[0]); */
  /* printf("OUT :%f \n",output[1]);  */
}

long double* forward_pass(neuralnet *nn, array3d_int * gamestate){

  pass_int_input(nn,gamestate);
  advance_state(nn);
  compute_act(nn);
  long double * output=get_output(nn);
  return output;
}

void float_forward_pass(neuralnet *nn, float* gamestate,float *output){

  pass_float_input(nn,gamestate);
  advance_state(nn);
  compute_act(nn);
  float_get_output(nn,output);
  /* printf("OUT2 :%f \n",output[0]); */
  /* printf("OUT2 :%f \n",output[1]); */
}


void mutate_weights(neuralnet *nn,double sigma){
    int indA,indW;
    for(int i=0;i<nn->MAX_NEURON;++i){
      for(int j=0;j<nn->MAX_NEURON;++j){
	indA=array3d_int_index(nn->A,i,j,0);
	indW=array3d_double_index(nn->W,i,j,0);
	if(nn->A->array[indA]){
	nn->W->array[indW]=nn->W->array[indW]+randn(0,sigma);
	}
      }
    }
  }

void full_mutate_weights(neuralnet *nn,double sigma){ 
  int indW;
  for(int i=0;i<nn->MAX_NEURON;++i){
    for(int j=0;j<nn->MAX_NEURON;++j){
      indW=array3d_double_index(nn->W,i,j,0);
      nn->W->array[indW]=nn->W->array[indW]+randn(0,sigma);
    }
  }
}


void neuralnet_mutate_table(neuralnet *nn,double sigma){
  for(int i=0;i<nn->Nneuron;++i){
    nn->table_act->array[i]+=randn(0,sigma);
  }
}



long double randn (double mu, double sigma)
{
  long double U1, U2, W, mult;
  static long double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (long double) X2);
    }
 
  do
    {
      U1 = -1 + ((long double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((long double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
   
  call = !call;
 
  return (mu + sigma * (long double) X1);
}
void full_neuralnet_free(neuralnet* nn ){
  free(nn->a->array);
  nn->a->array=NULL;
  free(nn->a);

  free(nn->W->array);
  nn->W->array=NULL;
  free(nn->W);
  
  free(nn->table_act->array);
  nn->table_act->array=NULL;
  free(nn->table_act);
  free(nn);

}

void neuralnet_free(neuralnet* nn ){
  free(nn->a->array);
  nn->a->array=NULL;
  free(nn->a);

  free(nn->W->array);
  nn->W->array=NULL;
  free(nn->W);

  free(nn->A->array);
  nn->A->array=NULL;
  free(nn->A);


  free(nn->table_act->array);
  nn->table_act->array=NULL;
  free(nn->table_act);
  
  free(nn);

}

  neuralnet * neuralnet_copy(neuralnet *nn){
    int sizeA=nn->MAX_NEURON*nn->MAX_NEURON*sizeof(int);
    int   sizeW=nn->MAX_NEURON*nn->MAX_NEURON*sizeof(long double);
    int   sizea=nn->MAX_NEURON*2*sizeof(long double);
    neuralnet* cop=malloc(sizeof(neuralnet)+sizeA+sizeW+sizea);
    if((cop = (neuralnet *) malloc(BUFSIZ)) == NULL) {
      printf("malloc error in neuralnet_copy");
      return 0;
    }

    cop->t=nn->t;
    cop->old=nn->old;
    cop->cur=nn->cur;
    cop->MAX_NEURON=nn->MAX_NEURON;
    cop->MAX_LINKS=nn->MAX_LINKS;
    cop->Nhidden=nn->Nhidden;
    cop->Ninput=nn->Ninput;
    cop->Noutput=nn->Noutput;
    cop->Nneuron=nn->Nneuron;
    cop->Nlinks=nn->Nlinks;

    cop->A=array3d_int_copy(nn->A);
    cop->W=array3d_double_copy(nn->W);
    cop->a=array3d_double_copy(nn->a);
    return cop;
  }
void neuralnet_replace(neuralnet *destination ,neuralnet *source ){
  //printf("nn_repl for 0 \n");
  assert(destination->MAX_NEURON==source->MAX_NEURON && destination->MAX_LINKS==source->MAX_LINKS);
  //  printf("nn_repl for 1 \n");
  destination->t=source->t;
  destination->old=source->old;
  destination->cur=source->cur;
  destination->MAX_NEURON=source->MAX_NEURON;
  destination->MAX_LINKS=source->MAX_LINKS;
  destination->Nhidden=source->Nhidden;
  destination->Ninput=source->Ninput;
  destination->Noutput=source->Noutput;
  destination->Nneuron=source->Nneuron; 
  destination->Nlinks=source->Nlinks;
  //  printf("nn_repl for 2 \n");
  array3d_int_replace(destination->A,source->A);
  //printf("nn_repl for 3 \n");
  array3d_double_replace(destination->W,source->W);
  //printf("nn_repl for 4 \n");
  array3d_double_replace(destination->a,source->a);
  //printf("nn_repl for 5 \n");
  array3d_double_replace(destination->table_act,source->table_act);
  //printf("nn_repl for out \n");
  
}

void full_neuralnet_replace(neuralnet *destination ,neuralnet *source ){
  //printf("nn_repl for 0 \n");
  assert(destination->MAX_NEURON==source->MAX_NEURON && destination->MAX_LINKS==source->MAX_LINKS);
  //  printf("nn_repl for 1 \n");
  destination->t=source->t;
  destination->old=source->old;
  destination->cur=source->cur;
  destination->MAX_NEURON=source->MAX_NEURON;
  destination->MAX_LINKS=source->MAX_LINKS;
  destination->Nhidden=source->Nhidden;
  destination->Ninput=source->Ninput;
  destination->Noutput=source->Noutput;
  destination->Nneuron=source->Nneuron; 
  destination->Nlinks=source->Nlinks;
  //  printf("nn_repl for 2 \n");
  //printf("nn_repl for 3 \n");
  array3d_double_replace(destination->W,source->W);
  //printf("nn_repl for 4 \n");
  array3d_double_replace(destination->a,source->a);
  //printf("nn_repl for 5 \n");
  array3d_double_replace(destination->table_act,source->table_act);
  //printf("nn_repl for out \n");
  
}

void neuralnet_write(neuralnet * nn){
  FILE* fW,*fA,*fa,*ftable;
    assert((fW=fopen("./txt/W.txt","wb"))!=NULL);
    assert((fA=fopen("./txt/A.txt","wb"))!=NULL);
    assert((fa=fopen("./txt/a.txt","wb"))!=NULL);
    assert((ftable=fopen("./txt/table_act.txt","wb"))!=NULL);
    
    array3d_double_write(fW,nn->W);
    array3d_int_write(fA,nn->A);
    array3d_double_write(fa,nn->a);
    array3d_double_write(ftable,nn->table_act);

    fclose(fW);
    fclose(fA);
    fclose(fa);
    fclose(ftable);
  }


void full_neuralnet_write(neuralnet * nn){
  FILE* fW,*fa,*ftable;
    assert((fW=fopen("./txt/W.txt","wb"))!=NULL);
    assert((fa=fopen("./txt/a.txt","wb"))!=NULL);
    assert((ftable=fopen("./txt/table_act.txt","wb"))!=NULL);
    
    array3d_double_write(fW,nn->W);
    array3d_double_write(fa,nn->a);
    array3d_double_write(ftable,nn->table_act);

    fclose(fW);
    fclose(fa);
    fclose(ftable);
  }


  void neuralnet_read(neuralnet * nn){
    FILE* fW,*fA,*fa,*ftable;
    assert((fW=fopen("./txt/W.txt","rb"))!=NULL);
    assert((fA=fopen("./txt/A.txt","rb"))!=NULL);
    assert((fa=fopen("./txt/a.txt","rb"))!=NULL);
    assert((ftable=fopen("./txt/table_act.txt","rb"))!=NULL);
	
    array3d_double_read(fW,nn->W);
    array3d_int_read(fA,nn->A);
    array3d_double_read(fa,nn->a);
    array3d_double_read(ftable,nn->table_act);
	
    fclose(fW);
    fclose(fA);
    fclose(fa);
    fclose(ftable);
    
  }
  void full_neuralnet_read(neuralnet * nn){
    FILE* fW,*fa,*ftable;
    assert((fW=fopen("./txt/W.txt","rb"))!=NULL);
    assert((fa=fopen("./txt/a.txt","rb"))!=NULL);
    assert((ftable=fopen("./txt/table_act.txt","rb"))!=NULL);
    
    array3d_double_read(fW,nn->W);
    array3d_double_read(fa,nn->a);
    array3d_double_read(ftable,nn->table_act);
	
    fclose(fW);
    fclose(fa);
    fclose(ftable);
    
  }

void neuralnet_write2(neuralnet * nn,FILE* fW,FILE* fA,FILE* fa,FILE * ftable){
    array3d_double_write(fW,nn->W);
    array3d_int_write(fA,nn->A);
    array3d_double_write(fa,nn->a);
    array3d_double_write(ftable,nn->table_act);
  }
void neuralnet_read2(neuralnet * nn,FILE* fW,FILE* fA,FILE* fa,FILE * ftable){
    array3d_double_read(fW,nn->W);
    array3d_int_read(fA,nn->A);
    array3d_double_read(fa,nn->a);
    array3d_double_read(ftable,nn->table_act);
  }


void full_neuralnet_write2(neuralnet * nn,FILE* fW,FILE* fa,FILE * ftable){
  array3d_double_write(fW,nn->W);
  array3d_double_write(fa,nn->a);
  array3d_double_write(ftable,nn->table_act);
}
void full_neuralnet_read2(neuralnet * nn,FILE* fW,FILE* fa,FILE * ftable){
  array3d_double_read(fW,nn->W);
  array3d_double_read(fa,nn->a);
  array3d_double_read(ftable,nn->table_act);
}
