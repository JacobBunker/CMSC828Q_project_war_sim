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

#include <pthread.h>


#include "neuralnet.h"
#include "GA.h"
#include "maze_simulator.h"

#define pi 3.142857 
#define NUM_PLAYERS 4
#define PLAYERS_PER_TEAM 2
#define DEBUG_PRINT 0
#define N_INPUT		4
#define N_OUTPUT	2

#define CHANCE_MUTATION 0.8

#define MAX_THREADS 20
#define DEBUG_THREADS 0

SimulationState render_sim;
cpVect render_pos;
cpVect render_vel;

char *name;

int popsize=200,
MAX_NEURON=100,
MAX_LINKS=2500,
Ngame=10,
Learning_time=100;

void InitSim(SimulationState *sim, int player_net_param) {
	sim->timeStep = 1.0/20.0;
	sim->arena_time_max = 10.0;
	sim->target = cpv(27.0,22.0);//cpv(40.0, 40.0);
	sim->spawn = cpv(10.0, 10.0);
	sim->force_multiplier = 10.0;

	int obstacle_type = 0;
	int team_one_type = 1;
	int team_two_type = 2;

	cpGroup team_one_group = team_one_type;
	cpGroup team_two_group = team_two_type;

	cpShapeFilter team_one_filter = cpShapeFilterNew(team_one_group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);
	cpShapeFilter team_two_filter = cpShapeFilterNew(team_two_type, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);

	sim->space = cpSpaceNew();

	cpFloat radius = 1;
	cpFloat mass = 1;

	// The moment of inertia is like mass for rotation
	// Use the cpMomentFor*(+) functions to help you approximate it.
	cpFloat moment = cpMomentForCircle(mass, 0, radius, cpvzero);

	sim->player_bodies = malloc(sizeof(cpBody *) * NUM_PLAYERS);
	sim->player_shapes = malloc(sizeof(cpShape *) * NUM_PLAYERS);

	for(int i = 0; i < NUM_PLAYERS; ++i) {
		sim->player_bodies[i] = cpSpaceAddBody(sim->space, cpBodyNew(mass, moment));
		sim->player_shapes[i] = cpSpaceAddShape(sim->space, cpCircleShapeNew(sim->player_bodies[i], radius, cpvzero));
		cpShapeSetFriction(sim->player_shapes[i], 0.7);
		if(i >= PLAYERS_PER_TEAM) {
			cpShapeSetFilter(sim->player_shapes[i], team_two_filter);
			cpShapeSetUserData(sim->player_shapes[i], &team_two_type);
		} else {
			cpShapeSetFilter(sim->player_shapes[i], team_one_filter);
			cpShapeSetUserData(sim->player_shapes[i], &team_one_type);
		}
	}

	sim->obstacle_number = 14;
	int raw_obstacle_vertices_count[14] = {8,4,8,4,4,4,8,4,6,4,4,4,4,4};
	sim->obstacle_vertices_count = malloc(sizeof(int) * sim->obstacle_number);
	for(int i=0; i < sim->obstacle_number; ++i) {
		sim->obstacle_vertices_count[i] = raw_obstacle_vertices_count[i];
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


	sim->obstacles = malloc(sizeof(cpVect *) * sim->obstacle_number);
	int c = 0;
	for(int i=0; i < sim->obstacle_number; ++i) {
		sim->obstacles[i] = malloc(sizeof(cpVect) * sim->obstacle_vertices_count[i]);
		for(int ii=0; ii < sim->obstacle_vertices_count[i]; ++ii) {
			sim->obstacles[i][ii] = cpv(raw_obstacle_vertices[c+(ii*2)],raw_obstacle_vertices[c+(ii*2)+1]);
		}
		c += sim->obstacle_vertices_count[i]*2;
	}

  	//add a polygon body+shape
	sim->obstacle_bodies = malloc(sizeof(cpBody *) * sim->obstacle_number);
	sim->obstacle_shapes = malloc(sizeof(cpShape *) * sim->obstacle_number);

	for(int i=0; i < sim->obstacle_number; ++i) {
		sim->obstacle_bodies[i] = cpSpaceGetStaticBody(sim->space);
		sim->obstacle_shapes[i] = cpPolyShapeNew(sim->obstacle_bodies[i], sim->obstacle_vertices_count[i],
												 sim->obstacles[i], cpTransformIdentity, radius);
		cpShapeSetFriction(sim->obstacle_shapes[i], 1.0f);
		cpSpaceAddShape(sim->space, sim->obstacle_shapes[i]);
		cpShapeSetUserData(sim->obstacle_shapes[i], &obstacle_type);
	}

	sim->players=malloc(NUM_PLAYERS*(sizeof(neuralnet)+player_net_param));
	for(int i=0;i<NUM_PLAYERS;++i) {
		sim->players[i] = neuralnet_init(N_INPUT, N_OUTPUT, MAX_NEURON, MAX_LINKS);
	} 
	//sim->players=malloc(sizeof(neuralnet *)*NUM_PLAYERS);

	sim->render = 1;
	sim->time = 0.0;
	sim->input = malloc(sizeof(float)*N_INPUT);
	sim->output = malloc(sizeof(float)*N_OUTPUT);

	for(int i=0;i<N_INPUT;++i){sim->input[i]=0.0;}
	for(int i=0;i<N_INPUT;++i){sim->output[i]=32.0;}
}

void FreeSim(SimulationState *sim) {

	free(sim->output);
	free(sim->input);

	for(int i=0;i<NUM_PLAYERS;++i){
		neuralnet_free(sim->players[i]);
	}

	free(sim->players);

	for(int i=sim->obstacle_number-1; i >= 0; --i) {
		cpShapeFree(sim->obstacle_shapes[i]);
		free(sim->obstacles[i]);
	}

	free(sim->obstacle_shapes);
	free(sim->obstacle_bodies);
	free(sim->obstacles);
	for(int i = NUM_PLAYERS - 1; i >= 0; --i) {
		cpShapeFree(sim->player_shapes[i]);
		cpBodyFree(sim->player_bodies[i]);
	}

	free(sim->player_shapes);
	free(sim->player_bodies);
	cpSpaceFree(sim->space);
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
}
  
void timer( int value )
{
    glutTimerFunc( 16, timer, 0 );
    glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	if(render_sim.time < render_sim.arena_time_max) {
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			render_sim.pos = cpBodyGetPosition(render_sim.player_bodies[current_player]);
			render_sim.vel = cpBodyGetVelocity(render_sim.player_bodies[current_player]);

			render_sim.input[0] = render_sim.pos.x;
			render_sim.input[1] = render_sim.pos.y;
			render_sim.input[2] = render_sim.vel.x;
			render_sim.input[3] = render_sim.vel.y;

			float_forward_pass(render_sim.players[current_player],render_sim.input,render_sim.output);

			render_sim.vel.x = render_sim.force_multiplier*render_sim.output[0];
			render_sim.vel.y = render_sim.force_multiplier*render_sim.output[1];
			//cpBodySetVelocity(sim.player_bodies[current_player], vel);
			cpBodySetForce(render_sim.player_bodies[current_player], render_sim.vel);
		}
		cpSpaceStep(render_sim.space, render_sim.timeStep);
		render_sim.time += render_sim.timeStep;

		//sleep(sim.timeStep);
		char buffer[64];
		snprintf(buffer, sizeof buffer, "%f", render_sim.time);
		glutSetWindowTitle(buffer);


		//RENDERING

		glColor3f(0.0f, 0.0f, 1.0f);
		for(int i=0; i < render_sim.obstacle_number; ++i) {
			glBegin(GL_POLYGON);
			for(int ii=0; ii < render_sim.obstacle_vertices_count[i]; ++ii) {
				glVertex3f(render_sim.obstacles[i][ii].x, render_sim.obstacles[i][ii].y, 0.0f);
			}
			glEnd();
		}
		glColor3f(0.0f, 1.0f, 0.0f);
		glPointSize(40.0);
		glBegin(GL_POINTS);
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			render_sim.pos = cpBodyGetPosition(render_sim.player_bodies[current_player]);
			glVertex3f(render_sim.pos.x, render_sim.pos.y, 0.0f);
		}
		glPointSize(5.0);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(render_sim.target.x, render_sim.target.y, 0.0f);
		glEnd(); 
		glutSwapBuffers();

	} else if(render_sim.render == 1){
		char *name = "SIMULATION ENDED";
		glutSetWindowTitle(name);
		// Clean up our ARENA objects
		FreeSim(&render_sim);
		render_sim.render = 0;
	}  
}

void* SimulationThread(void *arg) {

	SimulationState *sim = (SimulationState *) arg;
	float *outcome = malloc(sizeof(float)*NUM_PLAYERS);

	for(int p_i = 0; p_i < NUM_PLAYERS; ++p_i) {
		cpBodySetPosition(sim->player_bodies[p_i], sim->spawn);
		outcome[p_i] = 0.0;
	}

	for(cpFloat time = 0; time < sim->arena_time_max; time += sim->timeStep){
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
			sim->vel = cpBodyGetVelocity(sim->player_bodies[current_player]);

			//outcome[current_player] += -sqrt(((pos.x - target.x)*(pos.x - target.x)) +
			//								((pos.y - target.y)*(pos.y - target.y)));

			if(DEBUG_PRINT) {
				printf("Time is %5.2f. player %d body is at (%5.2f, %5.2f). It's velocity is (%5.2f, %5.2f)\n",
				time, current_player, sim->pos.x, sim->pos.y, sim->vel.x, sim->vel.y);
			}

			sim->input[0] = sim->pos.x;
			sim->input[1] = sim->pos.y;
			sim->input[2] = sim->vel.x;
			sim->input[3] = sim->vel.y;

			float_forward_pass(sim->players[current_player],sim->input,sim->output);
			//printf("outputs: %.2f, %.2f\n", output[0], output[1]);
			sim->vel.x = sim->force_multiplier*sim->output[0];
			sim->vel.y = sim->force_multiplier*sim->output[1];
			//printf("\t\tforce x: %f, force y: %f\n", sim->vel.x, sim->vel.y);

			//cpBodySetVelocity(sim->player_bodies[current_player], vel);
			cpBodySetForce(sim->player_bodies[current_player], sim->vel);
		}
		cpSpaceStep(sim->space, sim->timeStep);
	}

	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
		outcome[current_player] = (-sqrt(((sim->pos.x - sim->target.x)*(sim->pos.x - sim->target.x)) +
										 ((sim->pos.y - sim->target.y)*(sim->pos.y - sim->target.y))));
	}

	return outcome;
}

int main2(int argc, char** argv) {
	neuralnet *p1 = neuralnet_init(1,1,10,10);
	neuralnet *p2 = neuralnet_init(1,1,10,10);
	neuralnet_replace(p2,p1);
	neuralnet_free(p1);
	neuralnet_free(p2);	
	return 0;
}

int main(int argc, char** argv) {

  	/* NEURAL NETWORK INITIALIZATION */

	srand(time(NULL));

	GA *ga = GA_init(popsize,N_INPUT,N_OUTPUT,MAX_NEURON,MAX_LINKS);
	FILE * ffit,* fsig;  
	ffit=fopen("fit.txt","w");
	fsig=fopen("sig.txt","w");

	int sizeA=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(int);
	int sizeW=ga->MAX_NEURON*ga->MAX_NEURON*sizeof(long double);  
	int sizea=ga->MAX_NEURON*2*sizeof(long double); 
	//INITIALIZING THE SIMULATION HERE! 
	SimulationState sims[MAX_THREADS];
	for(int t=0;t<MAX_THREADS;++t) {
		InitSim(&sims[t], sizeA+sizea+sizeW);
	}

	InitSim(&render_sim, sizeA+sizea+sizeW);

	pthread_t threads[MAX_THREADS];
	void *res;

	//
	assert(ga->n%NUM_PLAYERS==0);
	int schedule[ga->n];

	for (int i=0;i<Learning_time;++i) {  
		printf("round: %d\n",i);

		/*RACE FUNCTION (GA *ga,int Ngame)*/
		ga->Nimproved_fit=0;
		for (int i=0;i<ga->n;++i){
			ga->copy_fit[i]=0;
		}
		//could generate full schedule for all games rather than just for one Ngame
		for(int i=0;i<Ngame;++i){
			permute(schedule,ga->n);
			int j = 0;
			while(j < (ga->n)/NUM_PLAYERS) {
				if(j + MAX_THREADS < (ga->n)/NUM_PLAYERS) {
					if(DEBUG_THREADS) { printf("creating threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<NUM_PLAYERS;++ii){
							neuralnet_replace(sims[t].players[ii],ga->pop[schedule[NUM_PLAYERS*j+ii]]);
							//sims[t].players[ii] = ga->pop[schedule[NUM_PLAYERS*j+ii]];
						}
						if(DEBUG_THREADS) { printf("\t\tready\n");}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					if(DEBUG_THREADS) { printf("joining threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						pthread_join(threads[t], &res);
						for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
							//printf("res: %f\n", ((float *) res)[current_player]);
	    					ga->copy_fit[schedule[NUM_PLAYERS*j+current_player]]+=((float *) res)[current_player];
	    				}
						free(res);
					}
					j = j + MAX_THREADS;
				}
				else {  //case when games left is less than MAX_THREADS
					int leftover = ((ga->n)/NUM_PLAYERS) - 1 - j;
					for(int t=0;t<leftover;++t) {
						for(int ii=0;ii<NUM_PLAYERS;++ii){
							//sims[t].players[ii] = ga->pop[schedule[NUM_PLAYERS*j+ii]];
							neuralnet_replace(sims[t].players[ii],ga->pop[schedule[NUM_PLAYERS*j+ii]]);
						}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					for(int t=0;t<leftover;++t) {
						pthread_join(threads[t], &res);
						for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
	    					ga->copy_fit[schedule[NUM_PLAYERS*j+current_player]]+=((float *) res)[current_player];
	    				}
						free(res);
					}
					j = ((ga->n)/NUM_PLAYERS);
				}
			}
			mutate_sigma(ga);
			GA_mutate_weights(ga,CHANCE_MUTATION);
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

		tournament_selection(ga);
		out_fit(ffit,ga);  
		out_sig(fsig,ga);
	}


  	/*TRAINING DONE*/

	printf("Training done \n\n");
	neuralnet_write(ga->pop[max_fit(ga)]);

	//set up the render_sim arena, get best agents
	for(int i = 0; i < NUM_PLAYERS; ++i) {
		render_sim.players[i] = neuralnet_init(N_INPUT,N_OUTPUT,MAX_NEURON,MAX_LINKS);
		neuralnet_replace(render_sim.players[i],ga->pop[n_best(ga,i)]);
		cpBodySetPosition(render_sim.player_bodies[i], render_sim.spawn);
	}

	/*clean the training equipment*/
	out_fit(ffit,ga);
	out_sig(fsig,ga);
	GA_free(ga);

	for(int t=0;t<MAX_THREADS;++t) {
			FreeSim(&sims[t]);
	}



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