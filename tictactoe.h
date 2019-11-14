#ifndef TICTACTOE_H
#define TICTACTOE_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "arrays.h"
#include "neuralnet.h"

_Bool check_win(int *game_state);
_Bool check_full(int *game_state);
void PrintGrid(int *game_state);
void ttt_state_mod(long double *output,array3d_int * gamestate,int player);
#endif
  
int game_of_ttt(neuralnet * P1,neuralnet * P2);
int show_game_of_ttt(neuralnet * P1,neuralnet * P2);
void human_vs_nn(neuralnet *nn);
