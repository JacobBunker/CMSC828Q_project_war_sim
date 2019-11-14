#include "tictactoe.h"

_Bool check_win(int *game_state   ){
  // check if game is won
  // vertical win
  if(game_state[0]==game_state[3] && game_state[0]==game_state[6] && game_state[0]!=0){return game_state[0];}
  if((game_state[1]==game_state[4]) && game_state[1]==game_state[7] && (game_state[1]!=0)){return game_state[1];}
  if(game_state[2]==game_state[5] && game_state[2]==game_state[8] && game_state[2]!=0){return game_state[2];}
  // horizontal win
  if(game_state[0]==game_state[1] && game_state[0]==game_state[2] && game_state[0]!=0){return game_state[0];}
  if(game_state[3]==game_state[4] && game_state[3]==game_state[5] && game_state[4]!=0){return game_state[3];}
  if(game_state[6]==game_state[7] && game_state[6]==game_state[8] && game_state[6]!=0){return game_state[6];}
  // diagonal win
  if(game_state[0]==game_state[4] && game_state[0]==game_state[8] && game_state[0]!=0){return game_state[0];}
  if(game_state[2]==game_state[4] && game_state[2]==game_state[6] && game_state[2]!=0){return game_state[2];}  
  return 0;
}

_Bool check_full(int *game_state){
  for(int i=0;i<9;++i){
    if(game_state[i]==0){
      return 0;
    }
  }
  return 1;
}



void PrintGrid(int *game_state){
  //User interface
  char game_symb[10];
  for(int i=0;i<9;++i){
    if(game_state[i]==0){game_symb[i]='.';}
    if(game_state[i]==1){game_symb[i]='X';}
    if(game_state[i]==2){game_symb[i]='O';}
  }
  game_symb[9]='\0';
  
  printf("Player 1: X\tvs\tPlayer 2: O\n");

  printf("     |     |     \n");
  printf("  %c  |  %c  |  %c  \n", game_symb[0], game_symb[1], game_symb[2]);
  printf("_____|_____|_____\n");
  printf("     |     |     \n");
  printf("  %c  |  %c  |  %c \n", game_symb[3], game_symb[4], game_symb[5]);
  printf("_____|_____|_____\n");
  printf("     |     |     \n");
  printf("  %c  |  %c  |  %c \n", game_symb[6], game_symb[7], game_symb[8]);
  printf("     |     |     \n\n");
}


void ttt_state_mod(long double* output ,array3d_int * gamestate,int player){

  int index[]={0,1,2,3,4,5,6,7,8,9};
  index_sort(output,index);
  /* printf("index sorted :"); */
  /* for (int i=0;i<9;++i){ */
  /*   printf("%d",index[i]); */
  /* } */
  /* printf(" \n "); */

  /* printf("output sorted :"); */
  /* for (int i=0;i<9;++i){ */
  /*   printf("% LF ",output[i]); */
  /* } */
  /* printf("\n"); */
    int i=0,notyet=1;

  while(notyet){
    if(gamestate->array[index[i]]==0){
      gamestate->array[index[i]]=player;
      notyet=0;
    }
    else{
      ++i;
    }

  }
  free(output);
}


  
  int game_of_ttt(neuralnet * P1,neuralnet * P2){
    neuralnet * players[]={P1,P2};
    int cross= rand()%2;
    int circle= (cross+1)%2;
    array3d_int *gamestate=array3d_int_init(9,1,1);
    array3d_int *formated_gamestate=array3d_int_init(9,1,1);
	
    while(1){
      ttt_state_mod(forward_pass(players[cross],gamestate),gamestate,1);
      //PrintGrid(gamestate->array);
      
      if(check_win(gamestate->array)){
	//	PrintGrid(gamestate->array);
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	return cross;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	//printf("full1\n\n\n\n\n");
	return 2;}
      
      for(int i=0; i<9;++i){
	switch(gamestate->array[i]){
	case 1:
	  formated_gamestate->array[i]=2;
	  break;
	case 2:
	  formated_gamestate->array[i]=1;
	  break;
	}
      }
       ttt_state_mod(forward_pass(players[circle],formated_gamestate),gamestate,2);
       //PrintGrid(gamestate->array);
      if(check_win(gamestate->array)){
	//	PrintGrid(gamestate->array);
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	/* printf("win2\n") */;
	return circle;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	//	printf("full2\n\n\n\n\n\n\n\n");
	return 2;} 
    }      
  }



  int show_game_of_ttt(neuralnet * P1,neuralnet * P2){
    neuralnet * players[]={P1,P2};
    int cross= rand()%2;
    int circle= (cross+1)%2;
    array3d_int *gamestate=array3d_int_init(9,1,1);
    array3d_int *formated_gamestate=array3d_int_init(9,1,1);
	
    while(1){
      ttt_state_mod(forward_pass(players[cross],gamestate),gamestate,1);
            PrintGrid(gamestate->array);
      
      if(check_win(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	return cross;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	/* printf("full1\n") */
	return 3;}
      
      for(int i=0; i<9;++i){
	switch(gamestate->array[i]){
	case 1:
	  formated_gamestate->array[i]=2;
	  break;
	case 2:
	  formated_gamestate->array[i]=1;
	  break;
	}
      }
       ttt_state_mod(forward_pass(players[circle],formated_gamestate),gamestate,2);
       PrintGrid(gamestate->array);
      if(check_win(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	/* printf("win2\n") */;
	return circle;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	free(formated_gamestate->array);
	free(formated_gamestate);
	/*printf("full2\n")*/;
	return 3;} 
    }      
  }


    
void human_vs_nn(neuralnet *nn){ 
  array3d_int *gamestate;
  int human,again=1;
  while(again){
    gamestate=array3d_int_init(9,1,1);
    while(1){
      ttt_state_mod(forward_pass(nn,gamestate),gamestate,1);
      PrintGrid(gamestate->array);
      
      if(check_win(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	printf("win cross\n"); 
        break;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	/* printf("full1\n") */
	break;}

      
      printf("play :\n");
      scanf("%d",&human);

      
      gamestate->array[human]=2;
      PrintGrid(gamestate->array);
      if(check_win(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	printf("win circle\n");
        break;}
      if(check_full(gamestate->array)){
	free(gamestate->array);
	free(gamestate);
	/* printf("full1\n") */
        break;}    

    }
    printf("other one ? 1 or 0 \n");
    scanf("%d",&again);
  }
}
