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

#include <chipmunk.h>

#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>


#include "neuralnet.h"
#include "GA.h"
#include "maze_simulator.h"

#define pi 3.142857 
#define NUM_PLAYERS 4
#define DEBUG_PRINT 0
#define N_INPUT		2
#define N_OUTPUT	4

SimulationState sim;
cpVect pos;
cpVect vel;

char *name;

void InitSim(int player_net_param) {
	sim.timeStep = 1.0/20.0;
	sim.arena_time_max = 20.0;
	sim.target = cpv(40.0, 40.0);
	sim.spawn = cpv(10.0, 10.0);

	int obstacle_type = 0;
	int team_one_type = 1;
	//int team_two_type = 2;

	//cpGroup player_group = team_one_type;
	//cpShapeFilter player_filter = cpShapeFilterNew(player_group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);

	sim.space = cpSpaceNew();

	cpFloat radius = 1;
	cpFloat mass = 1;

	// The moment of inertia is like mass for rotation
	// Use the cpMomentFor*(+) functions to help you approximate it.
	cpFloat moment = cpMomentForCircle(mass, 0, radius, cpvzero);

	sim.player_bodies = malloc(sizeof(cpBody *) * NUM_PLAYERS);
	sim.player_shapes = malloc(sizeof(cpShape *) * NUM_PLAYERS);

	for(int i = 0; i < NUM_PLAYERS; ++i) {
		sim.player_bodies[i] = cpSpaceAddBody(sim.space, cpBodyNew(mass, moment));
		sim.player_shapes[i] = cpSpaceAddShape(sim.space, cpCircleShapeNew(sim.player_bodies[i], radius, cpvzero));
		cpShapeSetFriction(sim.player_shapes[i], 0.7);
		//cpShapeSetFilter(sim.player_shapes[i], player_filter);
		cpShapeSetUserData(sim.player_shapes[i], &team_one_type);
	}

	sim.obstacle_number = 14;
	int raw_obstacle_vertices_count[14] = {8,4,8,4,4,4,8,4,6,4,4,4,4,4};
	sim.obstacle_vertices_count = malloc(sizeof(int) * sim.obstacle_number);
	for(int i=0; i < sim.obstacle_number; ++i) {
		sim.obstacle_vertices_count[i] = raw_obstacle_vertices_count[i];
	}
	float raw_obstacle_vertices[(8+4+8+4+4+4+8+4+6+4+4+4+4+4)*2] = {
		5, 18,   5, 20,   15, 20,  18, 13,  18, 5,   16, 5,   16, 13,  14, 18,
		23, 0,   23, 6,   25, 6,   25, 0,
		22, 17,  22, 29,  25, 29,  25, 23,  27, 23,  27, 20,  25, 20,  25, 17,
		31, 14,  31, 19,  36, 19,  36, 14,
		41, 5,   41, 30,  42, 30,  42, 5,
		32, 23,  32, 24,  34, 24,  34, 23,
		44, 27,  36, 27,  31, 35,  31, 40,  32, 40,  32, 35,  36, 31,  44, 28,
		23, 36,  19, 36,  19, 39,  23, 41,
		16, 28,  6, 28,   6, 31,   14, 31,  14, 40,  16, 40,
		7, 37,   4, 37,   4, 45,   7, 45,
		-5, -5,  -5, 0,   55, 0,   55, -5,
		55, -5,  50, -5,  50, 50,  55, 50,
		55, 45,  -5, 45,  -5, 50,  55, 50,
		0, -5,   -5, -5,  -5, 50,  0, 50 };


	sim.obstacles = malloc(sizeof(cpVect *) * sim.obstacle_number);
	int c = 0;
	for(int i=0; i < sim.obstacle_number; ++i) {
		sim.obstacles[i] = malloc(sizeof(cpVect) * sim.obstacle_vertices_count[i]);
		for(int ii=0; ii < sim.obstacle_vertices_count[i]; ++ii) {
			sim.obstacles[i][ii] = cpv(raw_obstacle_vertices[c+(ii*2)],raw_obstacle_vertices[c+(ii*2)+1]);
		}
		c += sim.obstacle_vertices_count[i]*2;
	}

  	//add a polygon body+shape
	sim.obstacle_bodies = malloc(sizeof(cpBody *) * sim.obstacle_number);
	sim.obstacle_shapes = malloc(sizeof(cpShape *) * sim.obstacle_number);

	for(int i=0; i < sim.obstacle_number; ++i) {
		sim.obstacle_bodies[i] = cpSpaceGetStaticBody(sim.space);
		sim.obstacle_shapes[i] = cpPolyShapeNew(sim.obstacle_bodies[i], sim.obstacle_vertices_count[i], sim.obstacles[i], cpTransformIdentity, radius);
		cpShapeSetFriction(sim.obstacle_shapes[i], 1.0f);
		cpSpaceAddShape(sim.space, sim.obstacle_shapes[i]);
		cpShapeSetUserData(sim.obstacle_shapes[i], &obstacle_type);
	}

	sim.players=malloc(NUM_PLAYERS*(sizeof(neuralnet)+player_net_param));
	sim.render = 1;
	sim.time = 0.0;
	sim.input = malloc(sizeof(float)*N_INPUT);
	sim.output = malloc(sizeof(float)*N_OUTPUT);

	for(int i=0;i<N_INPUT;++i){sim.input[i]=0.0;}
	for(int i=0;i<N_OUTPUT;++i){sim.output[i]=32.0;}
}

void FreeSim(SimulationState sim) {

	free(sim.output);
	free(sim.input);

	for(int i=0;i<NUM_PLAYERS;++i){
		neuralnet_free(sim.players[i]);
	}

	free(sim.players);

	for(int i=sim.obstacle_number-1; i >= 0; --i) {
		cpShapeFree(sim.obstacle_shapes[i]);
		free(sim.obstacles[i]);
	}
	free(sim.obstacle_shapes);
	free(sim.obstacle_bodies);
	free(sim.obstacles);
	for(int i = NUM_PLAYERS - 1; i >= 0; --i) {
		cpShapeFree(sim.player_shapes[i]);
		cpBodyFree(sim.player_bodies[i]);
	}
	free(sim.player_shapes);
	free(sim.player_bodies);
	cpSpaceFree(sim.space);
}



void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glColor3f(0.0, 0.0, 1.0);
	glMatrixMode(GL_PROJECTION);	
	glLoadIdentity();
	glOrtho(-10.0, 60.0, -10.0, 60.0, -1.0, 1.0);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glPointSize(40.0);
}
  
void timer( int value )
{
    glutTimerFunc( 16, timer, 0 );
    glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	if(sim.time < sim.arena_time_max) {
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			pos = cpBodyGetPosition(sim.player_bodies[current_player]);
			vel = cpBodyGetVelocity(sim.player_bodies[current_player]);

			sim.input[0] = pos.x;
			sim.input[1] = pos.y;
			sim.input[2] = vel.x;
			sim.input[3] = vel.y;

			float_forward_pass(sim.players[current_player],sim.input,sim.output);

			vel.x = 10*sim.output[0];
			vel.y = 10*sim.output[1];
			//cpBodySetVelocity(sim.player_bodies[current_player], vel);
			cpBodySetForce(sim.player_bodies[current_player], vel);
		}
		cpSpaceStep(sim.space, sim.timeStep);
		sim.time += sim.timeStep;

		//sleep(sim.timeStep);
		char buffer[64];
		snprintf(buffer, sizeof buffer, "%f", sim.time);
		glutSetWindowTitle(buffer);


		//RENDERING
		glColor3f(0.0f, 0.0f, 1.0f);
		for(int i=0; i < sim.obstacle_number; ++i) {
			glBegin(GL_POLYGON);
			for(int ii=0; ii < sim.obstacle_vertices_count[i]; ++ii) {
				glVertex3f(sim.obstacles[i][ii].x, sim.obstacles[i][ii].y, 0.0f);
			}
			glEnd();
		}
		glColor3f(0.0f, 1.0f, 0.0f);

		glBegin(GL_POINTS);
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			pos = cpBodyGetPosition(sim.player_bodies[current_player]);
			glVertex3f(pos.x, pos.y, 0.0f);
		}
		glEnd(); 
		glutSwapBuffers();

	} else if(sim.render == 1){
		char *name = "SIMULATION ENDED";
		glutSetWindowTitle(name);
		// Clean up our ARENA objects
		FreeSim(sim);
		sim.render = 0;
	}  
}

int main(int argc, char** argv) {

  	/* NEURAL NETWORK INITIALIZATION */

	srand(time(NULL));
	int popsize=40,
	MAX_NEURON=20,
	MAX_LINKS=200,
	Ngame=5,
	Learning_time=20;

	GA *ga = GA_init(popsize,N_INPUT,N_OUTPUT,MAX_NEURON,MAX_LINKS);
	FILE * ffit,* fsig;  
	ffit=fopen("fit.txt","w");
	fsig=fopen("sig.txt","w");

	int sizeA=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(int);
	int sizeW=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(long double);  
	int sizea=ga->MAX_NEURON*2*sizeof(long double); 
	//INITIALIZING THE SIMULATION HERE! 
	InitSim(sizeA+sizea+sizeW);
	//
	assert(ga->n%NUM_PLAYERS==0);
	int schedule[ga->n];
	float outcome[NUM_PLAYERS];

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
					sim.players[ii]=ga->pop[schedule[2*j+ii]];
				}

				//printf("arena match %d\n", j);
				//RESET THE ARENA
				for(int p_i = 0; p_i < NUM_PLAYERS; ++p_i) {
					cpBodySetPosition(sim.player_bodies[p_i], sim.spawn);
					outcome[p_i] = 0;
				}

	    		/* RUN THE SIMULATION */
				for(cpFloat time = 0; time < sim.arena_time_max; time += sim.timeStep){
					for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
						pos = cpBodyGetPosition(sim.player_bodies[current_player]);
						vel = cpBodyGetVelocity(sim.player_bodies[current_player]);

						//outcome[current_player] += -sqrt(((pos.x - target.x)*(pos.x - target.x)) +
						//								((pos.y - target.y)*(pos.y - target.y)));

						if(DEBUG_PRINT) {
							printf("Time is %5.2f. player %d body is at (%5.2f, %5.2f). It's velocity is (%5.2f, %5.2f)\n",
							time, current_player, pos.x, pos.y, vel.x, vel.y);
						}

						sim.input[0] = pos.x;
						sim.input[1] = pos.y;
						sim.input[2] = vel.x;
						sim.input[3] = vel.y;

						float_forward_pass(sim.players[current_player],sim.input,sim.output);
						//printf("outputs: %.2f, %.2f\n", output[0], output[1]);
						vel.x = 10*sim.output[0];
						vel.y = 10*sim.output[1];

						//cpBodySetVelocity(sim.player_bodies[current_player], vel);
						cpBodySetForce(sim.player_bodies[current_player], vel);
					}
					cpSpaceStep(sim.space, sim.timeStep);
				}

	    		/* END OF SIMULATION */
	    		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
	    			pos = cpBodyGetPosition(sim.player_bodies[current_player]);
	    			outcome[current_player] = -sqrt(((pos.x - sim.target.x)*(pos.x - sim.target.x)) +
														((pos.y - sim.target.y)*(pos.y - sim.target.y)));
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
		sim.players[i] = neuralnet_init(N_INPUT,N_OUTPUT,MAX_NEURON,MAX_LINKS);
		neuralnet_replace(sim.players[i],ga->pop[max_fit(ga)]);
		cpBodySetPosition(sim.player_bodies[i], sim.spawn);
	}

	/*clean the training equipment*/
	out_fit(ffit,ga);
	out_sig(fsig,ga);
	GA_free(ga);


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Simulator");
	init();
	glutTimerFunc(0, timer, 0);
	glutDisplayFunc(display);
	glutMainLoop();

	return 0;
}
