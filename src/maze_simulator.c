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

#define NUM_PLAYERS 16
#define PLAYERS_PER_TEAM 8
#define N_INPUT		10
#define N_OUTPUT	4
#define AGENT_SIZE 0.5
#define CHANCE_MUTATION 0.8
#define AGENT_HEALTH 30.0
#define AGENT_HIT_REWARD 100.0
#define MAX_THREADS 12


#define DEBUG_PRINT 0
#define DEBUG_THREADS 0
#define DEBUG_SCHEDULING 0
#define DEBUG_RES 0

SimulationState render_sim;
cpVect render_pos;
cpVect render_vel;

char *name;



/*
int popsize=64,
MAX_NEURON=100,
MAX_LINKS=(100*100)/4,
Ngame=5,
Learning_time=3; */


int popsize=160,
MAX_NEURON=100,
MAX_LINKS=(100*100),
Ngame=5,
Learning_time=20; 

void InitSim(SimulationState *sim, int player_net_param) {
	sim->timeStep = 1.0/20.0;
	sim->arena_time_max = 20.0;
	sim->target = cpv(27.0,22.0);//cpv(40.0, 40.0);
	sim->spawn_one = cpv(10.0, 10.0);
	sim->spawn_two = cpv(40.0, 40.0);
	sim->t_s_dist = 25;
	printf("max score: %.2f\n", sim->t_s_dist);
	sim->force_multiplier = 10.0 * AGENT_SIZE;
	sim->look_multiplier = 50.0;
	sim->scan_radius = 0.1;
	sim->look_number = 3;  //number of hitscan checks within the arc
	sim->look_spread = 0.523599 / 2.0; //size in radians of the arc which the agents look within
	int obstacle_type = 0;
	int agent_type = 1;
	int team_one_type = 1;
	int team_two_type = 2;

	sim->obstacle_number = 14;


	cpFloat radius = AGENT_SIZE;
	cpFloat mass = 1;


	sim->object_infos = malloc(sizeof(ObjectInfo) * NUM_PLAYERS * sim->obstacle_number);


	cpGroup team_one_group = team_one_type;
	cpGroup team_two_group = team_two_type;

	sim->team_one_filter = cpShapeFilterNew(team_one_group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);
	sim->team_two_filter = cpShapeFilterNew(team_two_group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);

	sim->space = cpSpaceNew();

	// The moment of inertia is like mass for rotation
	// Use the cpMomentFor*(+) functions to help you approximate it.
	cpFloat moment = cpMomentForCircle(mass, 0, radius, cpvzero);

	sim->player_bodies = malloc(sizeof(cpBody *) * NUM_PLAYERS);
	sim->player_shapes = malloc(sizeof(cpShape *) * NUM_PLAYERS);

	for(int i = 0; i < NUM_PLAYERS; ++i) {
		sim->object_infos[i].object_type = agent_type;
		sim->object_infos[i].object_id = i;
		sim->player_bodies[i] = cpSpaceAddBody(sim->space, cpBodyNew(mass, moment));
		sim->player_shapes[i] = cpSpaceAddShape(sim->space, cpCircleShapeNew(sim->player_bodies[i], radius, cpvzero));
		cpShapeSetFriction(sim->player_shapes[i], 0.7);
		
		if(i > PLAYERS_PER_TEAM) {
			sim->object_infos[i].team_id = 2;
			cpShapeSetFilter(sim->player_shapes[i], sim->team_two_filter);
			cpShapeSetUserData(sim->player_shapes[i], &(sim->object_infos[i]));
		} else {
			sim->object_infos[i].team_id = 1;
			cpShapeSetFilter(sim->player_shapes[i], sim->team_one_filter);
			cpShapeSetUserData(sim->player_shapes[i], &(sim->object_infos[i]));
		}
	}

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
		sim->object_infos[i+NUM_PLAYERS].object_type = obstacle_type;
		sim->object_infos[i+NUM_PLAYERS].object_id = i;
		sim->obstacle_bodies[i] = cpSpaceGetStaticBody(sim->space);
		sim->obstacle_shapes[i] = cpPolyShapeNew(sim->obstacle_bodies[i], sim->obstacle_vertices_count[i],
												 sim->obstacles[i], cpTransformIdentity, radius);
		cpShapeSetFriction(sim->obstacle_shapes[i], 1.0f);
		cpSpaceAddShape(sim->space, sim->obstacle_shapes[i]);
		cpShapeSetUserData(sim->obstacle_shapes[i], &(sim->object_infos[i+NUM_PLAYERS]));
	}

	sim->players=malloc(NUM_PLAYERS*(sizeof(neuralnet)+player_net_param));
	sim->player_scores=malloc(NUM_PLAYERS*sizeof(float));
	sim->player_previous_look=malloc(NUM_PLAYERS*sizeof(float));
	sim->player_previous_shoot=malloc(NUM_PLAYERS*sizeof(float));
	sim->player_health=malloc(NUM_PLAYERS*sizeof(float));

	for(int i=0;i<NUM_PLAYERS;++i) {
		sim->players[i] = neuralnet_init(N_INPUT, N_OUTPUT, MAX_NEURON, MAX_LINKS);
		sim->player_scores[i] = 0;
	} 
	//sim->players=malloc(sizeof(neuralnet *)*NUM_PLAYERS);

	sim->render = 1;
	sim->time = 0.0;
	sim->input = malloc(sizeof(float)*N_INPUT);
	sim->output = malloc(sizeof(float)*N_OUTPUT);

	for(int i=0;i<N_INPUT;++i){sim->input[i]=0.0;}
	for(int i=0;i<N_OUTPUT;++i){sim->output[i]=32.0;}

}

void FreeSim(SimulationState *sim) {

	free(sim->output);
	free(sim->input);

	for(int i=0;i<NUM_PLAYERS;++i){
		neuralnet_free(sim->players[i]);
	}

	free(sim->player_scores);
	free(sim->player_health);
	free(sim->player_previous_shoot);
	free(sim->player_previous_look);
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
	free(sim->object_infos);

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
	glEnable(GL_BLEND);


}
  
void timer( int value )
{
    glutTimerFunc( 16, timer, 0 );
    glutPostRedisplay();
}

//takes [-1.0,1.0] value and change value and returns corresponding radian value between -pi and pi
float angle_helper(float angle, float change) { 
	if(angle + change < -1.0) {
		return angle_helper(angle + change + 2.0, 0.0);
	} else if (angle + change > 1.0) {
		return angle_helper(angle + change - 2.0, 0.0);
	} else {
		return (pi * (angle + change));
	}
}

void DrawLookLine(float x1, float y1, float x2, float y2) {
	glColor4f(1.0f, 1.0f, 0.0f, 0.3f);
	glBegin(GL_LINES);   
	glVertex3f(x1,y1,0.0);     
	glVertex3f(x2,y2,0.0);
	glEnd();
}

void DrawShootLine(float x1, float y1, float x2, float y2) {
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);   
	glVertex3f(x1,y1,0.0);     
	glVertex3f(x2,y2,0.0);
	glEnd();
}

void ResetRenderSim() {
	render_sim.time = 0.0;

	for(int i = 0; i < NUM_PLAYERS; ++i) {
		cpBodySetVelocity(render_sim.player_bodies[i], cpv(0,0));
		if(i >= PLAYERS_PER_TEAM) {
			cpBodySetPosition(render_sim.player_bodies[i], render_sim.spawn_one);
		} else {
			cpBodySetPosition(render_sim.player_bodies[i], render_sim.spawn_two);
		}
		render_sim.player_previous_look[i] = 0.5;
		render_sim.player_previous_shoot[i] = 0.0;
		render_sim.player_health[i] = AGENT_HEALTH;
	}
}

void display(void)
{
	cpVect scan_target;
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 0.0f);
	//glBegin(GL_LINES);   
	if(render_sim.time < render_sim.arena_time_max) {
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			if(render_sim.player_health[current_player] > 0.0) {
				render_sim.pos = cpBodyGetPosition(render_sim.player_bodies[current_player]);
				render_sim.vel = cpBodyGetVelocity(render_sim.player_bodies[current_player]);

				render_sim.input[0] = render_sim.pos.x;
				render_sim.input[1] = render_sim.pos.y;
				render_sim.input[2] = render_sim.vel.x;
				render_sim.input[3] = render_sim.vel.y;

				if(1) {
					for(int l = 0; l < render_sim.look_number; ++l) {
						scan_target = cpv(render_sim.pos.x + (render_sim.look_multiplier*cos(angle_helper(render_sim.player_previous_look[current_player],(render_sim.look_spread / render_sim.look_number)*(l - (render_sim.look_number / 2.0))))),
								render_sim.pos.y + (render_sim.look_multiplier*sin(angle_helper(render_sim.player_previous_look[current_player],(render_sim.look_spread / render_sim.look_number)*(l - (render_sim.look_number / 2.0))))));
						render_sim.scan_hit = cpSpaceSegmentQueryFirst(render_sim.space, render_sim.pos, scan_target, 
							render_sim.scan_radius, cpShapeGetFilter(render_sim.player_shapes[current_player]), &(render_sim.scan_info));
						if(render_sim.scan_hit == NULL) {
							DrawLookLine(render_sim.pos.x, render_sim.pos.y, scan_target.x, scan_target.y);
							render_sim.input[4+l] = -1;
							render_sim.input[4+(l*2)+1] = 1;
						} else {
							render_sim.scan_hit_data = cpShapeGetUserData(render_sim.scan_hit);
							render_sim.hit_type = *((ObjectInfo *)render_sim.scan_hit_data);
							render_sim.input[4+(l*2)] = render_sim.hit_type.object_type;
							render_sim.input[4+(l*2)+1] = (render_sim.scan_info).alpha;
							DrawLookLine(render_sim.pos.x, render_sim.pos.y, render_sim.scan_info.point.x, render_sim.scan_info.point.y);
						}
					}
				}

				render_sim.input[4*(render_sim.look_number*2)] = render_sim.player_previous_shoot[current_player];

				float_forward_pass(render_sim.players[current_player],render_sim.input,render_sim.output);

				render_sim.vel.x = render_sim.force_multiplier*render_sim.output[0];
				render_sim.vel.y = render_sim.force_multiplier*render_sim.output[1];
				//cpBodySetVelocity(sim.player_bodies[current_player], vel);
				cpBodySetForce(render_sim.player_bodies[current_player], render_sim.vel);

				render_sim.player_previous_look[current_player] = render_sim.output[2];

				if(render_sim.output[3] > 0.0) {
					scan_target = cpv(render_sim.pos.x + (render_sim.look_multiplier*cos(angle_helper(render_sim.player_previous_look[current_player],0.0))),
								render_sim.pos.y + (render_sim.look_multiplier*sin(angle_helper(render_sim.player_previous_look[current_player],0.0))));
					render_sim.scan_hit = cpSpaceSegmentQueryFirst(render_sim.space, render_sim.pos, scan_target,
						 	render_sim.scan_radius, cpShapeGetFilter(render_sim.player_shapes[current_player]), &(render_sim.scan_info));
					if(render_sim.scan_hit == NULL) {
							render_sim.player_previous_shoot[current_player] = -1.0;
							DrawShootLine(render_sim.pos.x, render_sim.pos.y, scan_target.x, scan_target.y);

					} else {
							render_sim.scan_hit_data = cpShapeGetUserData(render_sim.scan_hit);
							render_sim.hit_type = *((ObjectInfo *) render_sim.scan_hit_data);
							if(render_sim.hit_type.object_type == 1) { //hit an enemy agent
								render_sim.player_previous_shoot[current_player] = 1.0;
								render_sim.player_health[render_sim.hit_type.object_id] -= 1.0;
							} else {
								render_sim.player_previous_shoot[current_player] = -1.0;
							}
							DrawShootLine(render_sim.pos.x, render_sim.pos.y, render_sim.scan_info.point.x, render_sim.scan_info.point.y);
					}
				}
			}
			else {
				cpBodySetPosition(render_sim.player_bodies[current_player], cpv(-10, -10));
			}

		}
		//glEnd();

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
		glPointSize(40.0 * AGENT_SIZE);
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
		//set up the render_sim arena, get best agents
		ResetRenderSim();
		//FreeSim(&render_sim);
		//render_sim.render = 0;
	}  
}

void* SimulationThread(void *arg) {

	SimulationState *sim = (SimulationState *) arg;
	float *outcome = malloc(sizeof(float)*NUM_PLAYERS);

	for(int p_i = 0; p_i < NUM_PLAYERS; ++p_i) {
		cpBodySetVelocity(sim->player_bodies[p_i], cpv(0,0));
		if(p_i >= PLAYERS_PER_TEAM) {
			cpBodySetPosition(sim->player_bodies[p_i], sim->spawn_one);
		} else {
			cpBodySetPosition(sim->player_bodies[p_i], sim->spawn_two);
		}
		outcome[p_i] = 0.0;
		sim->player_scores[p_i] = 0.0;
		sim->player_previous_look[p_i] = 0.0;
		sim->player_health[p_i] = AGENT_HEALTH;
	}


	for(cpFloat time = 0; time < sim->arena_time_max; time += sim->timeStep){
		for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
			if(sim->player_health[current_player] > 0.0) {
				sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
				sim->vel = cpBodyGetVelocity(sim->player_bodies[current_player]);

				if(DEBUG_PRINT) {
					printf("Time is %5.2f. player %d body is at (%5.2f, %5.2f). It's velocity is (%5.2f, %5.2f)\n",
					time, current_player, sim->pos.x, sim->pos.y, sim->vel.x, sim->vel.y);
				}

				sim->input[0] = sim->pos.x;
				sim->input[1] = sim->pos.y;
				sim->input[2] = sim->vel.x;
				sim->input[3] = sim->vel.y;

				//printf("look: %f\n", sim->player_previous_look[current_player]);

				if(1) {
					for(int l = 0; l < sim->look_number; ++l) {
						sim->scan_hit = cpSpaceSegmentQueryFirst(sim->space, sim->pos, 
							cpv(sim->pos.x + (sim->look_multiplier*cos(angle_helper(sim->player_previous_look[current_player],(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0))))),
								sim->pos.y + (sim->look_multiplier*sin(angle_helper(sim->player_previous_look[current_player],(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0)))))),
						 	sim->scan_radius, cpShapeGetFilter(sim->player_shapes[current_player]), &(sim->scan_info));
						if(DEBUG_PRINT) { 
							printf("scan from %.2f,%.2f to %.2f,%.2f\n", sim->pos.x, sim->pos.y, 
						 		cos(angle_helper(sim->player_previous_look[current_player],(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0)))),
						 		sin(angle_helper(sim->player_previous_look[current_player],(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0)))));
						}
						if(sim->scan_hit == NULL) {
							sim->input[4+(l*2)] = -1;
							sim->input[4+(l*2)+1] = 1;
						} else {
							sim->scan_hit_data = cpShapeGetUserData(sim->scan_hit);
							sim->hit_type = *((ObjectInfo *) sim->scan_hit_data);
							sim->input[4+(l*2)] = sim->hit_type.object_type;
							sim->input[4+(l*2)+1] = (sim->scan_info).alpha;
						}
					}
				}

				sim->input[4*(sim->look_number*2)] = sim->player_previous_shoot[current_player];

				float_forward_pass(sim->players[current_player],sim->input,sim->output);
				//printf("outputs: %.2f, %.2f\n", output[0], output[1]);
				sim->vel.x = sim->force_multiplier*sim->output[0];
				sim->vel.y = sim->force_multiplier*sim->output[1];
				sim->player_previous_look[current_player] = sim->output[2];

				if(sim->output[3] > 0.0) {
					sim->scan_hit = cpSpaceSegmentQueryFirst(sim->space, sim->pos, 
							cpv(sim->pos.x + (sim->look_multiplier*cos(angle_helper(sim->player_previous_look[current_player],0.0))),
								sim->pos.y + (sim->look_multiplier*sin(angle_helper(sim->player_previous_look[current_player],0.0)))),
						 	sim->scan_radius, cpShapeGetFilter(sim->player_shapes[current_player]), &(sim->scan_info));
					if(sim->scan_hit == NULL) {
							sim->player_previous_shoot[current_player] = -1.0;
					} else {
							sim->scan_hit_data = cpShapeGetUserData(sim->scan_hit);
							sim->hit_type = *((ObjectInfo *) sim->scan_hit_data);
							if(sim->hit_type.object_type == 1) { //hit an enemy agent
								sim->player_previous_shoot[current_player] = 1.0;
								sim->player_scores[current_player] += AGENT_HIT_REWARD;
								sim->player_scores[sim->hit_type.object_id] -= AGENT_HIT_REWARD;
								sim->player_health[sim->hit_type.object_id] -= 1.0;
							} else {
								sim->player_previous_shoot[current_player] = -1.0;
							}
					}
				}

				//printf("\t\tforce x: %f, force y: %f\n", sim->vel.x, sim->vel.y);

				//cpBodySetVelocity(sim->player_bodies[current_player], vel);
				cpBodySetForce(sim->player_bodies[current_player], sim->vel);
			}
			else {
				cpBodySetPosition(sim->player_bodies[current_player], cpv(-10, -10));
			}
			
		}
		cpSpaceStep(sim->space, sim->timeStep);
	}

	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
		outcome[current_player] = sim->player_scores[current_player]; 
		//+ (sim->t_s_dist - sqrt(((sim->pos.x - sim->target.x)*(sim->pos.x - sim->target.x)) + ((sim->pos.y - sim->target.y)*(sim->pos.y - sim->target.y))));
	}

	return outcome;
}

int main(int argc, char** argv) {

  	/* NEURAL NETWORK INITIALIZATION */
	float high_score = -1000.0;
	printf("popsize: %d, MAX_NEURON: %d, MAX_LINKS: %d, Ngame: %d, Learning_time: %d\n", popsize, MAX_NEURON, MAX_LINKS, Ngame, Learning_time);
	printf("CURRENT HIGH SCORE: %.2f\n", high_score);

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
		printf("round: %d out of %d, %.2f done\n",i,Learning_time,((float)i/(float)Learning_time));

		/*RACE FUNCTION (GA *ga,int Ngame)*/
		ga->Nimproved_fit=0;
		for (int i=0;i<ga->n;++i){
			ga->copy_fit[i]=0;
		}
		//could generate full schedule for all games rather than just for one Ngame
		for(int i=0;i<Ngame;++i){
			permute(schedule,ga->n);
			if(DEBUG_SCHEDULING) { 
				printf("schedule: ");
				for(int db=0; db < ga->n; ++db) {
					printf(" %d,", schedule[db]);
				}
				printf("\n");
			}
			int j = 0;
			while(j < (ga->n)/NUM_PLAYERS) {
				if(DEBUG_THREADS) { printf("j: %d, ga->n/NUM_PLAYERS:%d\n\t", j, (ga->n)/NUM_PLAYERS); }
				if((ga->n)/NUM_PLAYERS >= j + MAX_THREADS) {
					if(DEBUG_THREADS) { printf("creating threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<NUM_PLAYERS;++ii){
							neuralnet_replace(sims[t].players[ii],ga->pop[schedule[(NUM_PLAYERS*(j+t))+ii]]);
							if(DEBUG_SCHEDULING) { printf(" %d -> %d,", (NUM_PLAYERS*(j+t))+ii, schedule[(NUM_PLAYERS*(j+t))+ii]); }
							//sims[t].players[ii] = ga->pop[schedule[NUM_PLAYERS*j+ii]];
						}
						if(DEBUG_THREADS) { printf("\t\tready\n");}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					if(DEBUG_THREADS) { printf("joining threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						pthread_join(threads[t], &res);
						if(DEBUG_RES) { printf("res:"); }
						for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
							if(DEBUG_RES) { printf(" %f,\n", ((float *) res)[current_player]); }
	    					ga->copy_fit[schedule[(NUM_PLAYERS*(j+t))+current_player]]+=((float *) res)[current_player];
	    					if(((float *) res)[current_player] > high_score) {
	    						high_score = ((float *) res)[current_player];
								printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    					}
	    				}
	    				if(DEBUG_RES) { printf("\n"); }
						free(res);
					}
					j = j + MAX_THREADS;
				}
				else {  //case when games left is less than MAX_THREADS
					int leftover = ((ga->n)/NUM_PLAYERS) - 1 - j;
					for(int t=0;t<leftover;++t) {
						for(int ii=0;ii<NUM_PLAYERS;++ii){
							//sims[t].players[ii] = ga->pop[schedule[NUM_PLAYERS*j+ii]];
							neuralnet_replace(sims[t].players[ii],ga->pop[schedule[(NUM_PLAYERS*(j+t))+ii]]);
							if(DEBUG_SCHEDULING) { printf(" %d -> %d,", (NUM_PLAYERS*(j+t))+ii, schedule[(NUM_PLAYERS*(j+t))+ii]); }
						}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					for(int t=0;t<leftover;++t) {
						pthread_join(threads[t], &res);
						for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
							if(DEBUG_RES) { printf(" %f,\n", ((float *) res)[current_player]); }
	    					ga->copy_fit[schedule[(NUM_PLAYERS*(j+t))+current_player]]+=((float *) res)[current_player];
	    					if(((float *) res)[current_player] > high_score) {
	    						high_score = ((float *) res)[current_player];
								printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    					}
	    				}
						free(res);
					}
					j = ((ga->n)/NUM_PLAYERS);
				}
				if(DEBUG_SCHEDULING) { printf("\n"); }
			}
			//mutate_sigma(ga);
			//GA_mutate_weights(ga,CHANCE_MUTATION);
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

		mutate_sigma(ga);
		GA_mutate_weights(ga,CHANCE_MUTATION);

		out_fit(ffit,ga);  
		out_sig(fsig,ga);
	}


  	/*TRAINING DONE*/

	printf("Training done \n\n");
	neuralnet_write(ga->pop[max_fit(ga)]);

	//get best agents
	for(int i = 0; i < NUM_PLAYERS; ++i) {
		render_sim.players[i] = neuralnet_init(N_INPUT,N_OUTPUT,MAX_NEURON,MAX_LINKS);
		neuralnet_replace(render_sim.players[i],ga->pop[n_best(ga,i)]);
	}

	ResetRenderSim();
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