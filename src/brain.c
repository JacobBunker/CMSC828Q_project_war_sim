#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "brain.h"

#define SPREAD 1
#define DEBUG_HARD 0


brain * brain_init(int Ninput,int Noutput,int Ncluster,int Size_cluster,int Ncluster_links){
  assert(Ninput + Noutput + 1<Size_cluster) ;
  if(DEBUG_HARD){
    printf("init\n");
  }
  // add assert number of links
  
  int intra_sizeA=Size_cluster*Size_cluster*sizeof(int),
    intra_sizeW=Size_cluster*Size_cluster*sizeof(long double),
    intra_sizea=Size_cluster*2*sizeof(long double),
    intra_sizetable=Size_cluster*sizeof(long double),
    sizeA=Ncluster*Ncluster*sizeof(int), 
    sizeW=Ncluster*Ncluster*sizeof(long double); 
  
  int intrann_size=sizeof(neuralnet)+intra_sizea+intra_sizeW+intra_sizeA+intra_sizetable; // 
 
  brain * br=malloc(sizeof(brain)+intrann_size*Ncluster+sizeA+sizeW);
  
  br->Ncluster=Ncluster;
  br->Noutput=Noutput;
  br->Ninput=Ninput;
  br->Size_cluster=Size_cluster;
  br->Ncluster_links=Ncluster_links;
  br->cluster=malloc(Ncluster*intrann_size);
  
  br->cluster[0]=neuralnet_full_init(Ninput,0,Size_cluster); //input nn
  br->cluster[Ncluster-1]=neuralnet_full_init(0,Noutput,Size_cluster); //output nn
  
  for(int i=1;i<Ncluster-1;++i){ 
    br->cluster[i]=neuralnet_full_init(0,0,Size_cluster);
  }
  br->A=array3d_int_init(Ncluster,Ncluster,1);
  br->W=array3d_double_init(Ncluster_links,Size_cluster,1);


  
  long double temp;
  
  int a,b,Notyet,indA,indW;
  int Nlinks=0;
  assert(Ncluster_links< br->Ncluster * br->Ncluster - br->Ncluster);
  while(Nlinks<Ncluster_links){
    Notyet=1;
    while(Notyet){
      // Loop that initiate weights at random.
      a=0,b=0;
      while(a==b){
	// this loop make sure a!=b so that a neuron is not linked to itself
	a=rand() %br->Ncluster;
	b=rand() %br->Ncluster;
      }
      indA=array3d_int_index(br->A,a,b,0);
      if(br->A->array[indA]==0){
	br->A->array[indA]=1;
	for(int i =0;i<Size_cluster;++i){
	  indW=array3d_double_index(br->W,Nlinks,i,0);
	  temp=(long double) (SPREAD*random())/ INT_MAX -SPREAD/2;
	  br->W->array[indW]=temp;
	}
	++Nlinks;
	Notyet=0;
      }
    }
  }
  return br;
}


brain * brain_graph_init(int Ninput,int Noutput,int nx,int ny,int Size_cluster,int Ncluster_links){
  if(DEBUG_HARD){
    printf("graph_init\n");
  }
  assert(Ninput + Noutput + 1<Size_cluster &&
	 Ncluster_links+(nx-1)*ny+(ny-1)*nx<nx*nx*ny*ny  ) ; 
  // add assert number of links
  int Ncluster=nx*ny;
  
  int intra_sizeA=Size_cluster*Size_cluster*sizeof(int),
    intra_sizeW=Size_cluster*Size_cluster*sizeof(long double),
    intra_sizea=Size_cluster*2*sizeof(long double),
    intra_sizetable=Size_cluster*sizeof(long double),
    sizeA=Ncluster*Ncluster*sizeof(int), 
    sizeW=Ncluster*Ncluster*sizeof(long double); 
  
  int intrann_size=sizeof(neuralnet)+intra_sizea+intra_sizeW+intra_sizeA+intra_sizetable; //  
  int br_size=sizeof(brain)+intrann_size*Ncluster+sizeA+sizeW;
  brain * br=malloc(br_size);

  br->br_size=br_size;
  br->Ncluster=Ncluster;
  br->Noutput=Noutput;
  br->Ninput=Ninput;
  br->Size_cluster=Size_cluster;
  br->Ncluster_links=Ncluster_links;
  br->cluster=malloc(Ncluster*intrann_size);
  br->Ntotal_links=2*((nx-1)*ny+(ny-1)*nx)+Ncluster_links;
  
  br->cluster[0]=neuralnet_full_init(Ninput,0,Size_cluster); //input nn
  br->cluster[Ncluster-1]=neuralnet_full_init(0,Noutput,Size_cluster); //output nn
  
  for(int i=1;i<Ncluster-1;++i){ 
    br->cluster[i]=neuralnet_full_init(1,1,Size_cluster);
  }
  br->A=graph_builder(nx,ny);
  br->W=array3d_double_init(br->Ntotal_links,Size_cluster,1);
  
  long double temp;
  
  int a,b,Notyet,indA,indW;
  int Nlinks=0;
  assert(br->Ntotal_links< br->Ncluster * br->Ncluster - br->Ncluster);
  while(Nlinks<Ncluster_links){
    Notyet=1;
    while(Notyet){
      // Loop that initiate weights at random.
      a=0,b=0;
      while(a==b){
	// this loop make sure a!=b so that a neuron is not linked to itself
	a=rand() %br->Ncluster;
	b=rand() %br->Ncluster;
      }
      indA=array3d_int_index(br->A,a,b,0);
      if(br->A->array[indA]==0){
	br->A->array[indA]=1;
	for(int i =0;i<Size_cluster;++i){
	  indW=array3d_double_index(br->W,Nlinks,i,0);
	  temp=(long double) (SPREAD*random())/ INT_MAX -SPREAD/2;
	  br->W->array[indW]=temp;
	}
	++Nlinks;
	Notyet=0;
      }
    }
  }
  int act=0;
  for(int i=0;i<Ncluster*Ncluster;++i){
    if(br->A->array[i]==1){
      ++act;
    }
  }


  return br;
}

   
void brain_show_act(brain *br){
  for(int i=0;i<br->Ncluster;++i){
    printf("\n Network %d\n",i);
    array3d_double_show(br->cluster[i]->a);
  }
}

void brain_show_table(brain *br){
  for(int i=0;i<br->Ncluster;++i){
    printf("\n Network %d\n",i);
    array3d_double_show(br->cluster[i]->table_act);
  }
}

void brain_outer_cluster_compute_lin(brain* br,int i,int j,int weight_ind){
  if(DEBUG_HARD){
    printf("outer_cluster\n");
  }
  /* Compute the lineaire addition created by the link */
  /* between cluster i to cluster j */
  assert(i!=j && i<br->Ncluster && j<br->Ncluster);
  int indai,indaj,indW;

  //  indai=array3d_double_index(br->cluster[i]->a,br->cluster[i]->Ninput,br->cluster[i]->old,0);
  //br->cluster[i]->a->array[indai]=1; 
  
  for(int ii=0;ii<br->Size_cluster;++ii){
    indai=array3d_double_index(br->cluster[i]->a,ii,br->cluster[i]->old,0);
    indaj=array3d_double_index(br->cluster[j]->a,ii,br->cluster[j]->cur,0);
    indW=array3d_double_index(br->W,weight_ind,ii,0);
    br->cluster[j]->a->array[indaj]+=br->W->array[indW]*br->cluster[i]->a->array[indai];
  }
}
void brain_forward_pass(brain *br,float* input,float* output){
  if(DEBUG_HARD){
    printf("br_forward_pass\n");
  }
  brain_pass_float_input(br,input);
  int indA=0,
    weight_ind=0;
  for(int i=0;i<br->Ncluster;++i){
    advance_state(br->cluster[i]);
    neuralnet_lin_computation(br->cluster[i]);
    for(int j=0;j<br->Ncluster;++j){
      indA=array3d_int_index(br->A,j,i,0);
      if(br->A->array[indA]){       
	brain_outer_cluster_compute_lin(br,j,i,weight_ind);
	++weight_ind;
      }
    }
    run_tanh(br->cluster[i]);
  }
  brain_float_get_output(br,output);
}


void brain_show_weights(brain *br){
  for(int i=0;i<br->Ncluster;++i){
    printf("\n Network %d\n",i);
    array3d_double_show(br->cluster[i]->W);
  }
}


void brain_free(brain *br){
  for(int i=0;i<br->Ncluster;++i){
    neuralnet_free(br->cluster[i]);
  }


  free(br->W->array);
  free(br->W);

  free(br->A->array);
  free(br->A);
  free(br->cluster);
  free(br);
}

void brain_replace(brain *destination ,brain *source ){
  //  printf("br_repl 0 \n");
  if(DEBUG_HARD){
    printf("br replace\n");
  }
  assert(destination->Ncluster==source->Ncluster &&
	 destination->Size_cluster==source->Size_cluster &&
	 destination->Ncluster_links==source->Ncluster_links &&
	 destination->Noutput==source->Noutput &&
	 destination->Ninput==source->Ninput);

  for(int i=0;i<destination->Ncluster;++i){
    //printf("br_repl for in %d \n",i);
    neuralnet_replace(destination->cluster[i],source->cluster[i]);
    //printf("br_repl for out %d \n",i);
  }
  // printf("br_repl 1 \n");
  array3d_double_replace(destination->W,source->W);
  //printf("br_repl 2 \n");
  array3d_int_replace(destination->A,source->A);
  //printf("br_repl out \n");
}


void brain_mutate_weights(brain *br,double sigma){
  int indW;
  if(DEBUG_HARD){
    printf("br mut we\n");
  }
  for(int i=0;i<br->Ncluster;++i){
    mutate_weights(br->cluster[i],sigma);
  }
    
  for(int i =0 ;i<br->Ntotal_links;++i){
    for(int j =0 ;j<br->Size_cluster;++j){
      indW=array3d_double_index(br->W,i,j,0);
      br->W->array[indW]+=randn(0,sigma);	
    }
  }
}

void brain_float_get_output(brain * br,float * output){
  if(DEBUG_HARD){
    printf("br get out\n");
  }
  int inda;
  for (int i=0;i<br->Noutput;++i){
    inda=array3d_double_index(br->cluster[br->Ncluster-1]->a,
			      br->cluster[br->Ncluster-1]->Ninput+br->cluster[br->Ncluster-1]->Nhidden+i,
			      br->cluster[br->Ncluster-1]->cur,0);
    output[i]=br->cluster[br->Ncluster-1]->a->array[inda];
  }
}

void brain_pass_float_input(brain* br, float* input){
  if(DEBUG_HARD){
    printf("br pass in\n");
  }
  int inda;
  for (int i=0;i<br->Ninput;++i){
    
    inda=array3d_double_index(br->cluster[0]->a
			      ,i,br->cluster[0]->cur,0);
    br->cluster[0]->a->array[inda]=input[i];
  }
}
 
array3d_int * graph_builder(int nx, int ny){
  if(DEBUG_HARD){
    printf("graph build\n");
  }
  array3d_int *  A=array3d_int_init(nx*ny,nx*ny,1);
  int indA;
  int notacorner(int i){
    return i!=0 && i!=(ny-1) && i!=ny*(nx-1) && i!=nx*ny-1;	
  }
  for(int i=0;i<nx*ny;++i){
    
    if(i%ny==0 && notacorner(i)){
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
    }
    else if (0<=i&& i<ny && notacorner(i)){
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
    }

    else if ((i-ny+1)%ny==0 && notacorner(i)){
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
    }

    else if( ny*(nx-1)<=i && i<nx*ny  && notacorner(i)){
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
    }
    
    else if(i==0){
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
    }
    else if(i==ny-1){
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
    }
    else if(i==ny*(nx-1)){
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
    }
    else if(i==nx*ny-1){
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
    }

    else{
      indA=array3d_int_index(A,i+ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-1,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i-ny,i,0);
      A->array[indA]=1;
      indA=array3d_int_index(A,i+1,i,0);
      A->array[indA]=1;
    }
  }
  return A;
}




void brain_reset_act(brain * br){
  if(DEBUG_HARD){
    printf("reset act\n");
  }
  for(int i=0;i<br->Ncluster;++i){
    neuralnet_reset_act(br->cluster[i]);
  }
}

void brain_mutate_table(brain * br, double sigma){
  if(DEBUG_HARD){
    printf("br mut table\n");
  }
  for(int i=0;i<br->Ncluster;++i){
    neuralnet_mutate_table(br->cluster[i],sigma);
  }  
}

void brain_write(brain * br){
  if(DEBUG_HARD){
    printf("br writer\n");
  }
  FILE *fcW,*fcA,*fca,*fW,*fA,*fctable;
  assert((fW=fopen("./txt/W.txt","wb"))!=NULL);
  assert((fA=fopen("./txt/A.txt","wb"))!=NULL);
  assert((fcW=fopen("./txt/fcW.txt","wb"))!=NULL);
  assert((fcA=fopen("./txt/fcA.txt","wb"))!=NULL);
  assert((fca=fopen("./txt/fca.txt","wb"))!=NULL);
  assert((fctable=fopen("./txt/fctable.txt","wb"))!=NULL);
  
  for(int i=0;i<br->Ncluster;++i){
    neuralnet_write2(br->cluster[i],fcW,fcA,fca,fctable); 
  }

  array3d_double_write(fW,br->W);
  array3d_int_write(fA,br->A);

  fclose(fW);
  fclose(fA);
  fclose(fca);
  fclose(fcW);
  fclose(fcA);
  fclose(fctable);
}




void brain_read(brain * br){
  FILE *fcW,*fcA,*fca,*fW,*fA,*fctable;
  assert((fW=fopen("./txt/W.txt","rb"))!=NULL);
  assert((fA=fopen("./txt/A.txt","rb"))!=NULL);
  assert((fcW=fopen("./txt/fcW.txt","rb"))!=NULL);
  assert((fcA=fopen("./txt/fcA.txt","rb"))!=NULL);
  assert((fca=fopen("./txt/fca.txt","rb"))!=NULL);
  assert((fctable=fopen("./txt/fctable.txt","rb"))!=NULL);
  
  for(int i=0;i<br->Ncluster;++i){
    neuralnet_read2(br->cluster[i],fcW,fcA,fca,fctable); 
  }

  array3d_double_read(fW,br->W);
  array3d_int_read(fA,br->A);

  fclose(fW);
  fclose(fA);
  fclose(fca);
  fclose(fcW);
  fclose(fcA);
  fclose(fctable);
}
