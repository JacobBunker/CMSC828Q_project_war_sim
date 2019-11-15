#include <stdio.h>
#include <chipmunk.h>
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

#define NUM_PLAYERS 4
#define DEBUG_PRINT 0

int main(void) {
  // cpVect is a 2D vector and cpv() is a shortcut for initializing them.
	//cpVect gravity = cpv(0, -100);
	int arena_type = 0;
	int obstacle_type = 1;
	int team_one_type = 2;

	cpGroup player_group = team_one_type;
	cpShapeFilter player_filter = cpShapeFilterNew(player_group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);

	// Create an empty space.
	cpSpace *space = cpSpaceNew();
	//cpSpaceSetGravity(space, gravity);

	cpFloat radius = 1;
	cpFloat mass = 1;

	// The moment of inertia is like mass for rotation
	// Use the cpMomentFor*(+) functions to help you approximate it.
	cpFloat moment = cpMomentForCircle(mass, 0, radius, cpvzero);

	cpBody *player_bodies[NUM_PLAYERS];
	cpShape *player_shapes[NUM_PLAYERS];
	for(int i = 0; i < NUM_PLAYERS; ++i) {
		player_bodies[i] = cpSpaceAddBody(space, cpBodyNew(mass, moment));
		player_shapes[i] = cpSpaceAddShape(space, cpCircleShapeNew(player_bodies[i], radius, cpvzero));
		cpShapeSetFriction(player_shapes[i], 0.7);
		cpShapeSetFilter(player_shapes[i], player_filter);
		cpShapeSetUserData(player_shapes[i], &team_one_type);
	}

	int obstacle_number = 14;
	int obstacle_vertices[14] = {8,4,8,4,4,4,8,4,6,4,4,4,4,4};
	cpVect *obstacles[obstacle_number];

	obstacles[0] = malloc(sizeof(cpVect) * obstacle_vertices[0]);
	obstacles[0][0] = cpv(5, 18);
	obstacles[0][1] = cpv(5, 20);
	obstacles[0][2] = cpv(15, 20);
	obstacles[0][3] = cpv(18, 13);
	obstacles[0][4] = cpv(18, 5);
	obstacles[0][5] = cpv(16, 5);
	obstacles[0][6] = cpv(16, 13);
	obstacles[0][7] = cpv(14, 18);

	obstacles[1] = malloc(sizeof(cpVect) * obstacle_vertices[1]);
	obstacles[1][0] = cpv(23, 0);
	obstacles[1][1] = cpv(23, 6);
	obstacles[1][2] = cpv(25, 6);
	obstacles[1][3] = cpv(25, 0);

	obstacles[2] = malloc(sizeof(cpVect) * obstacle_vertices[2]);
	obstacles[2][0] = cpv(22, 17);
	obstacles[2][1] = cpv(22, 29);
	obstacles[2][2] = cpv(25, 29);
	obstacles[2][3] = cpv(25, 23);
	obstacles[2][4] = cpv(27, 23);
	obstacles[2][5] = cpv(27, 20);
	obstacles[2][6] = cpv(25, 20);
	obstacles[2][7] = cpv(25, 17);

	obstacles[3] = malloc(sizeof(cpVect) * obstacle_vertices[3]);
	obstacles[3][0] = cpv(31, 14);
	obstacles[3][1] = cpv(31, 19);
	obstacles[3][2] = cpv(36, 19);
	obstacles[3][3] = cpv(36, 14);

	obstacles[4] = malloc(sizeof(cpVect) * obstacle_vertices[4]);
	obstacles[4][0] = cpv(41, 5);
	obstacles[4][1] = cpv(41, 30);
	obstacles[4][2] = cpv(42, 30);
	obstacles[4][3] = cpv(42, 5);

	obstacles[5] = malloc(sizeof(cpVect) * obstacle_vertices[5]);
	obstacles[5][0] = cpv(32, 23);
	obstacles[5][1] = cpv(32, 24);
	obstacles[5][2] = cpv(34, 24);
	obstacles[5][3] = cpv(34, 23);

	obstacles[6] = malloc(sizeof(cpVect) * obstacle_vertices[6]);
	obstacles[6][0] = cpv(44, 27);
	obstacles[6][1] = cpv(36, 27);
	obstacles[6][2] = cpv(31, 35);
	obstacles[6][3] = cpv(31, 40);
	obstacles[6][4] = cpv(32, 40);
	obstacles[6][5] = cpv(32, 35);
	obstacles[6][6] = cpv(36, 31);
	obstacles[6][7] = cpv(44, 28);

	obstacles[7] = malloc(sizeof(cpVect) * obstacle_vertices[7]);
	obstacles[7][0] = cpv(23, 36);
	obstacles[7][1] = cpv(19, 36);
	obstacles[7][2] = cpv(19, 39);
	obstacles[7][3] = cpv(23, 41);

	obstacles[8] = malloc(sizeof(cpVect) * obstacle_vertices[8]);
	obstacles[8][0] = cpv(16, 28);
	obstacles[8][1] = cpv(6, 28);
	obstacles[8][2] = cpv(6, 31);
	obstacles[8][3] = cpv(14, 31);
	obstacles[8][4] = cpv(14, 40);
	obstacles[8][5] = cpv(16, 40);

	obstacles[9] = malloc(sizeof(cpVect) * obstacle_vertices[9]);
	obstacles[9][0] = cpv(7, 37);
	obstacles[9][1] = cpv(4, 37);
	obstacles[9][2] = cpv(4, 45);
	obstacles[9][3] = cpv(7, 45);


	obstacles[10] = malloc(sizeof(cpVect) * obstacle_vertices[9]);
	obstacles[10][0] = cpv(-5, -5);
	obstacles[10][1] = cpv(-5, 0);
	obstacles[10][2] = cpv(55, 0);
	obstacles[10][3] = cpv(55, -5);
	
	obstacles[11] = malloc(sizeof(cpVect) * obstacle_vertices[9]);
	obstacles[11][0] = cpv(55, -5);
	obstacles[11][1] = cpv(50, -5);
	obstacles[11][2] = cpv(50, 50);
	obstacles[11][3] = cpv(55, 50);

	obstacles[12] = malloc(sizeof(cpVect) * obstacle_vertices[9]);
	obstacles[12][0] = cpv(55, 45);
	obstacles[12][1] = cpv(-5, 45);
	obstacles[12][2] = cpv(-5, 50);
	obstacles[12][3] = cpv(55, 50);

	obstacles[13] = malloc(sizeof(cpVect) * obstacle_vertices[9]);
	obstacles[13][0] = cpv(0, -5);
	obstacles[13][1] = cpv(-5, -5);
	obstacles[13][2] = cpv(-5, 50);
	obstacles[13][3] = cpv(0, 50);


  //add a polygon body+shape
	cpBody **obstacle_bodies = malloc(sizeof(cpBody *) * obstacle_number);
	cpShape **obstacle_shapes = malloc(sizeof(cpShape *) * obstacle_number);

	radius = 1;

	for(int i=0; i < obstacle_number; ++i) {
		obstacle_bodies[i] = cpSpaceGetStaticBody(space);
		obstacle_shapes[i] = cpPolyShapeNew(obstacle_bodies[i], obstacle_vertices[i], obstacles[i], cpTransformIdentity, radius);
		cpShapeSetFriction(obstacle_shapes[i], 1.0f);
		cpSpaceAddShape(space, obstacle_shapes[i]);
		cpShapeSetUserData(obstacle_shapes[i], &obstacle_type);
	}



	cpFloat timeStep = 1.0/20.0;
	cpFloat arena_time_max = 20.0;
	cpVect target = cpv(40.0, 40.0);
	cpVect spawn = cpv(10.0, 10.0);
  /* NEURAL NETWORK INITIALIZATION */

	srand(time(NULL));
	int popsize=4,
	Ninput=4,
	Noutput=2,
	MAX_NEURON=100,
	MAX_LINKS=5000,
	Ngame=5,
	Learning_time=2;

	GA *ga = GA_init(popsize,Ninput,Noutput,MAX_NEURON,MAX_LINKS);
	FILE * ffit,* fsig, *fren;  
	ffit=fopen("fit.txt","w");
	fsig=fopen("sig.txt","w");
	fren=fopen("ren.txt","w");

	int sizeA=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(int);
	int sizeW=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(long double);  
	int sizea=ga->MAX_NEURON*2*sizeof(long double);  
	neuralnet **players=malloc(NUM_PLAYERS*(sizeof(neuralnet)+sizeA+sizea+sizeW));
	assert(ga->n%NUM_PLAYERS==0);
	int schedule[ga->n];
	float outcome[NUM_PLAYERS];
	cpVect pos;
	cpVect vel;
	float input[Ninput];
	float output[Noutput];
	for(int i=0;i<Ninput;++i){input[i]=0.0;}
	for(int i=0;i<Noutput;++i){output[i]=32.0;}

	for (int i=0;i<Learning_time;++i) {  
		printf("round: %d\n",i);

		/*RACE FUNCTION (GA *ga,int Ngame)*/
		ga->Nimproved_fit=0;
		for (int i=0;i<ga->n;++i){
			ga->copy_fit[i]=0;
		}
		for(int i=0;i<Ngame;++i){
			permute(schedule,ga->n);
			for(int j=0;j<(ga->n)/NUM_PLAYERS;++j){
				for(int ii=0;ii<NUM_PLAYERS;++ii){
					players[ii]=ga->pop[schedule[2*j+ii]];
				}

				//printf("arena match %d\n", j);
				//RESET THE ARENA
				for(int p_i = 0; p_i < NUM_PLAYERS; ++p_i) {
					cpBodySetPosition(player_bodies[p_i], spawn);
					outcome[p_i] = 0;
				}

	    		/* RUN THE SIMULATION */
				for(cpFloat time = 0; time < arena_time_max; time += timeStep){
					for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
						pos = cpBodyGetPosition(player_bodies[current_player]);
						vel = cpBodyGetVelocity(player_bodies[current_player]);

						//outcome[current_player] += -sqrt(((pos.x - target.x)*(pos.x - target.x)) +
						//								((pos.y - target.y)*(pos.y - target.y)));

						if(DEBUG_PRINT) {
							printf("Time is %5.2f. player %d body is at (%5.2f, %5.2f). It's velocity is (%5.2f, %5.2f)\n",
							time, current_player, pos.x, pos.y, vel.x, vel.y);
						}

						input[0] = pos.x;
						input[1] = pos.y;
						input[2] = vel.x;
						input[3] = vel.y;

						float_forward_pass(players[current_player],input,output);
						//printf("outputs: %.2f, %.2f\n", output[0], output[1]);
						vel.x = 10*output[0];
						vel.y = 10*output[1];

						cpBodySetVelocity(player_bodies[current_player], vel);
					}
					cpSpaceStep(space, timeStep);
				}

	    		/* END OF SIMULATION */
	    		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
	    			pos = cpBodyGetPosition(player_bodies[current_player]);
	    			outcome[current_player] = -sqrt(((pos.x - target.x)*(pos.x - target.x)) +
														((pos.y - target.y)*(pos.y - target.y)));
	    			ga->copy_fit[schedule[2*j+current_player]]+=outcome[current_player];
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
		memcpy(ga->fit_array,ga->copy_fit,ga->n*sizeof(long double));
		/*END OF RACE FUNCTION*/

		mutate_sigma(ga);
		tournament_selection(ga);
		GA_mutate_weights(ga,.8);
		out_fit(ffit,ga);  
		out_sig(fsig,ga);
	}


  /*TRAINING DONE*/


	printf("Training done \n\n");
	neuralnet_write(ga->pop[max_fit(ga)]);

	//RESET THE ARENA, GET BEST PLAYERS
	for(int i = 0; i < NUM_PLAYERS; ++i) {
		//printf("MAX FIT: %d\n", max_fit(ga));
		//array3d_double_show(ga->pop[max_fit(ga)]->W);
		players[i] = neuralnet_init(Ninput,Noutput,MAX_NEURON,MAX_LINKS);
		//printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
		//array3d_double_show(ga->pop[max_fit(ga)]->W);

		neuralnet_replace(players[i],ga->pop[max_fit(ga)]);
		cpBodySetPosition(player_bodies[i], spawn);
	}

	/* RUN THE SIMULATION */
	for(cpFloat time = 0; time < arena_time_max; time += timeStep){
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			pos = cpBodyGetPosition(player_bodies[current_player]);
			vel = cpBodyGetVelocity(player_bodies[current_player]);

			input[0] = pos.x;
			input[1] = pos.y;
			input[2] = vel.x;
			input[3] = vel.y;

			float_forward_pass(players[current_player],input,output);

			vel.x = 10*output[0];
			vel.y = 10*output[1];

			cpBodySetVelocity(player_bodies[current_player], vel);
		}

		/* VISUALIZATION CODE HERE */
		fprintf(fren, "%.2f distance: %f, position: %f, %f\n", time, 
			sqrt(((pos.x - target.x)*(pos.x - target.x)) + ((pos.y - target.y)*(pos.y - target.y))),
			pos.x, pos.y);
		cpSpaceStep(space, timeStep);
	}

	/* END OF SIMULATION */

	for(int i=0;i<4;++i){
		neuralnet_free(players[i]);
	}

	free(players);
	out_fit(ffit,ga);
	out_sig(fsig,ga);
	GA_free(ga);
	fclose(ffit);
	fclose(fsig);
	fclose(fren);
  // Clean up our ARENA objects and exit!

	for(int i=obstacle_number-1; i >= 0; --i) {
		cpShapeFree(obstacle_shapes[i]);
		free(obstacles[i]);
	}
	for(int i = NUM_PLAYERS - 1; i >= 0; --i) {
		cpShapeFree(player_shapes[i]);
		cpBodyFree(player_bodies[i]);
	}
	cpSpaceFree(space);

	return 0;
}