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
#include "Brain_GA.h"
#include "brain.h"
#include "maze_simulator.h"

#define pi 3.1415926

#define NUM_PLAYERS 16
#define PLAYERS_PER_TEAM 8
#define N_INPUT		10
#define N_OUTPUT	4
#define AGENT_SIZE 0.5
#define CHANCE_MUTATION 0.8
#define AGENT_HEALTH 30.0
#define AGENT_HIT_REWARD 100.0
#define MAX_THREADS 8


#define FULL_MAX_NEURONS 100
#define FULL_MAX_LINKS 2500

#define DEBUG_PRINT 0
#define DEBUG_THREADS 0
#define DEBUG_TIME 0

SimulationState render_sim;

char *name;
char buffer[64];
cpVect spawns[NUM_PLAYERS/PLAYERS_PER_TEAM];
float team_colors[2][3];

int popsize=4*NUM_PLAYERS,
  nx=10,
  ny=5,
  Size_cluster=30,
  Ncluster_links=100,
  Ngame=4,
  Learning_time=20; 




void InitSim(SimulationState *sim) {
	sim->timeStep = 1.0/20.0;
	sim->arena_time_max = 20.0;
	sim->target = cpv(27.0,22.0);//cpv(40.0, 40.0);
	sim->t_s_dist = 25;
	printf("max score: %.2f\n", sim->t_s_dist);
	sim->force_multiplier = 10.0 * AGENT_SIZE;
	sim->look_multiplier = 50.0;
	sim->scan_radius = 0.1;
	sim->look_number = 3;  //number of hitscan checks within the arc
	sim->look_spread = 0.523599 / 2.0; //size in radians of the arc which the agents look within

	spawns[0] = cpv(10.0, 10.0);
	spawns[1] = cpv(40.0, 40.0);
	team_colors[0][0] = 0.0;
	team_colors[0][1] = 1.0;
	team_colors[0][2] = 0.5;
	team_colors[1][0] = 1.0;
	team_colors[1][1] = 0.4;
	team_colors[1][2] = 0.4;

	int obstacle_type = 0;
	int agent_type = 1;
	sim->obstacle_number = 14;

	cpFloat radius = AGENT_SIZE;
	cpFloat mass = 1;

	assert((NUM_PLAYERS % PLAYERS_PER_TEAM) == 0);

	sim->object_infos = malloc(sizeof(ObjectInfo) * NUM_PLAYERS * sim->obstacle_number);
	sim->space = cpSpaceNew();

	cpFloat moment = cpMomentForCircle(mass, 0, radius, cpvzero);

	sim->teams = malloc(sizeof(Team) * (NUM_PLAYERS / PLAYERS_PER_TEAM));
	for(int i = 0; i < (NUM_PLAYERS / PLAYERS_PER_TEAM); ++i) {
		if(i % 2 == 0) {
			sim->teams[i].type = 0;
		} else {
			sim->teams[i].type = 1;
		}
		sim->teams[i].id = i + 2; //skip zero and one since they are reserved for obstacles and agents
		sim->teams[i].group = sim->teams[i].id;
		sim->teams[i].filter = cpShapeFilterNew(sim->teams[i].group, CP_ALL_CATEGORIES, CP_ALL_CATEGORIES);
		sim->teams[i].spawn = spawns[i];
	}

	sim->players = malloc(sizeof(Player) * NUM_PLAYERS);
	int t = 0; int c = 0;
	for(int i = 0; i < NUM_PLAYERS; ++i) {
		sim->players[i].team = t;
		sim->object_infos[i].object_type = agent_type;
		sim->object_infos[i].object_id = i;
		sim->object_infos[i].team_id = sim->players[i].team;

		sim->players[i].body = cpSpaceAddBody(sim->space, cpBodyNew(mass, moment));
		sim->players[i].shape = cpSpaceAddShape(sim->space, cpCircleShapeNew(sim->players[i].body, radius, cpvzero));
		cpShapeSetFriction(sim->players[i].shape, 0.7);
		cpShapeSetFilter(sim->players[i].shape, sim->teams[sim->players[i].team].filter);
		cpShapeSetUserData(sim->players[i].shape, &(sim->object_infos[i]));

		sim->players[i].id = i;
		if(sim->teams[sim->players[i].team].type == 0) { //a fully connected player
			sim->players[i].br = neuralnet_init(N_INPUT, N_OUTPUT, FULL_MAX_NEURONS, FULL_MAX_LINKS);
			sim->players[i].type = 0;
		} else {
			sim->players[i].br = brain_graph_init(N_INPUT, N_OUTPUT, nx,ny,Size_cluster,Ncluster_links);
			sim->players[i].type = 1;
		}
		sim->players[i].score = 0;
		sim->players[i].prev_look = 0;
		sim->players[i].prev_shoot = 0;
		sim->players[i].health = 0;

		c = c + 1;
		if(c == PLAYERS_PER_TEAM) {
			t++;
			c = 0;
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
	c = 0;
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

	for(int i=sim->obstacle_number-1; i >= 0; --i) {
		cpShapeFree(sim->obstacle_shapes[i]);
		free(sim->obstacles[i]);
	}

	free(sim->obstacle_shapes);
	free(sim->obstacle_bodies);
	free(sim->obstacles);

	for(int i=0;i<NUM_PLAYERS;++i){
		if(sim->players[i].type == 0) {
			neuralnet_free((neuralnet *) sim->players[i].br);
		} else {
			brain_free((brain *) sim->players[i].br);
		}
		free(sim->players[i].shape);
		free(sim->players[i].body);
	}

	free(sim->players);
	free(sim->teams);
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

void ResetSim(SimulationState *sim) {
	sim->time = 0.0;

	for(int i = 0; i < NUM_PLAYERS; ++i) {
		cpBodySetVelocity(sim->players[i].body, cpv(0,0));
		cpBodySetPosition(sim->players[i].body, sim->teams[sim->players[i].team].spawn);
		sim->players[i].prev_look = 0.5;
		sim->players[i].prev_shoot = 0.0;
		sim->players[i].score = 0.0;
		sim->players[i].health = AGENT_HEALTH;
		if(sim->players[i].type == 0) {
			neuralnet_reset_act((neuralnet *)sim->players[i].br);
		} else {
			brain_reset_act((brain * )sim->players[i].br);
		}
	}
}

void StepSim(SimulationState *sim, int graphics_on, double *rec_time) {
	clock_t times[5]; 
	if(DEBUG_TIME) { times[0] = clock(); }
	//clock_t start, post_look, post_forward, post_shoot, end;
	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		if(sim->players[current_player].health > 0.0) {
			if(DEBUG_TIME) { times[1] = clock(); }

			sim->pos = cpBodyGetPosition(sim->players[current_player].body);
			sim->vel = cpBodyGetVelocity(sim->players[current_player].body);

			sim->input[0] = sim->pos.x;
			sim->input[1] = sim->pos.y;
			sim->input[2] = sim->vel.x;
			sim->input[3] = sim->vel.y;

			//printf("look: %f\n", sim->player_previous_look[current_player]);

			for(int l = 0; l < sim->look_number; ++l) {
				sim->scan_target = cpv(sim->pos.x + (sim->look_multiplier*cos(angle_helper(sim->players[current_player].prev_look,(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0))))),
										sim->pos.y + (sim->look_multiplier*sin(angle_helper(sim->players[current_player].prev_look,(sim->look_spread / sim->look_number)*(l - (sim->look_number / 2.0))))));
				sim->scan_hit = cpSpaceSegmentQueryFirst(sim->space, sim->pos, sim->scan_target, sim->scan_radius, cpShapeGetFilter(sim->players[current_player].shape), &(sim->scan_info));
				if(sim->scan_hit == NULL) {
					if(graphics_on) { DrawLookLine(sim->pos.x, sim->pos.y, sim->scan_target.x, sim->scan_target.y); }
					sim->input[4+(l*2)] = -1;
					sim->input[4+(l*2)+1] = 1;
				} else {
					sim->scan_hit_data = cpShapeGetUserData(sim->scan_hit);
					sim->hit_type = *((ObjectInfo *) sim->scan_hit_data);
					sim->input[4+(l*2)] = sim->hit_type.object_type;
					sim->input[4+(l*2)+1] = (sim->scan_info).alpha;
					if(graphics_on) { DrawLookLine(sim->pos.x, sim->pos.y, sim->scan_info.point.x, sim->scan_info.point.y); }
				}
			}
			sim->input[4*(sim->look_number*2)] = sim->players[current_player].prev_shoot;
			if(DEBUG_TIME) { times[2] = clock(); }
			if(sim->players[current_player].type == 0) {
				float_forward_pass((neuralnet * ) sim->players[current_player].br,sim->input,sim->output);
			} else {
				brain_forward_pass((brain * )sim->players[current_player].br,sim->input,sim->output);
			}
			if(DEBUG_TIME) { times[3] = clock(); }
			for(int i = 0; i < N_OUTPUT; ++i) {
				/* if(isnan(sim->output[i])) { */
				/* 	printf("output %d: %.3f\n", i, sim->output[i]); */
				/* } */
			}			

			sim->vel.x = sim->force_multiplier*sim->output[0];
			sim->vel.y = sim->force_multiplier*sim->output[1];
			sim->players[current_player].prev_look = sim->output[2];

			if(sim->output[3] > 0.0) {
				sim->scan_target = cpv(sim->pos.x + (sim->look_multiplier*cos(angle_helper(sim->players[current_player].prev_look,0.0))),
									   sim->pos.y + (sim->look_multiplier*sin(angle_helper(sim->players[current_player].prev_look,0.0))));
				sim->scan_hit = cpSpaceSegmentQueryFirst(sim->space, sim->pos, sim->scan_target,
					 	sim->scan_radius, cpShapeGetFilter(sim->players[current_player].shape), &(sim->scan_info));
				if(sim->scan_hit == NULL) {
						if(graphics_on) { DrawShootLine(sim->pos.x, sim->pos.y, sim->scan_target.x, sim->scan_target.y); }
						sim->players[current_player].prev_shoot = -1.0;
				} else {
						sim->scan_hit_data = cpShapeGetUserData(sim->scan_hit);
						sim->hit_type = *((ObjectInfo *) sim->scan_hit_data);
						if(sim->hit_type.object_type == 1) { //hit an enemy agent
							sim->players[current_player].prev_shoot = 1.0;
							sim->players[current_player].score += AGENT_HIT_REWARD;
							sim->players[sim->hit_type.object_id].score -= AGENT_HIT_REWARD;
							sim->players[sim->hit_type.object_id].health -= 1.0;
						} else {
							sim->players[current_player].prev_shoot = -1.0;
						}
						if(graphics_on) { DrawShootLine(sim->pos.x, sim->pos.y, sim->scan_info.point.x, sim->scan_info.point.y); }
				}
			}
			if(DEBUG_TIME) { 
				times[4] = clock(); 
				rec_time[1] += ((double) (times[2] - times[1])) / CLOCKS_PER_SEC;
				rec_time[2] += ((double) (times[3] - times[2])) / CLOCKS_PER_SEC;
				rec_time[3] += ((double) (times[4] - times[3])) / CLOCKS_PER_SEC;
			}
			cpBodySetForce(sim->players[current_player].body, sim->vel);
		}
		else {
			cpBodySetPosition(sim->players[current_player].body, cpv(-100, -100));
		}
		
	}
	if(DEBUG_TIME) { times[3] = clock(); }
	cpSpaceStep(sim->space, sim->timeStep);
	if(DEBUG_TIME) { times[4] = clock(); 
		rec_time[0] += ((double) (times[4] - times[0])) / CLOCKS_PER_SEC;
		rec_time[4] += ((double) (times[4] - times[3])) / CLOCKS_PER_SEC; 
	}
}

void display(void)
{
	double sim_time_rec[5]; 
	for(int i = 0; i < 5; ++i) {
		sim_time_rec[i] = 0.0;
	}

	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 0.0f);
	//glBegin(GL_LINES);   
	if(render_sim.time < render_sim.arena_time_max) {
	
		StepSim(&render_sim, 1, sim_time_rec);
		render_sim.time += render_sim.timeStep;

		//sleep(sim.timeStep);
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
			glColor3f(team_colors[render_sim.players[current_player].team][0],
					  team_colors[render_sim.players[current_player].team][1],
					  team_colors[render_sim.players[current_player].team][2]);
			render_sim.pos = cpBodyGetPosition(render_sim.players[current_player].body);
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
		ResetSim(&render_sim);
	}  
}

void* SimulationThread(void *arg) {
	SimulationState *sim = (SimulationState *) arg;
	//NUM_PLAYERS + 0: total time of simulation
	//NUM_PLAYERS + 1: time performing sensor looks
	//NUM_PLAYERS + 2: time performing network pass
	//NUM_PLAYERS + 3: time performing shoot
	//NUM_PLAYERS + 4: time performing physics space step

	//clock_t start, post_look, post_forward, post_shoot, end;

	double sim_time_rec[5]; 
	for(int i = 0; i < 5; ++i) {
		sim_time_rec[i] = 0.0;
	}

	float *outcome = malloc(sizeof(float)*(NUM_PLAYERS + 5)); //+5 is for time keeping slots

	for(int p_i = 0; p_i < NUM_PLAYERS + 5; ++p_i) {
		outcome[p_i] = 0.0;
	}

	ResetSim(sim);

	for(cpFloat time = 0; time < sim->arena_time_max; time += sim->timeStep){
		StepSim(sim, 0, sim_time_rec);
	}

	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		//sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
		outcome[current_player] = sim->players[current_player].score; 
		//+ (sim->t_s_dist - sqrt(((sim->pos.x - sim->target.x)*(sim->pos.x - sim->target.x)) + ((sim->pos.y - sim->target.y)*(sim->pos.y - sim->target.y))));
	}


	if(DEBUG_TIME) { 
		for(int i = 0; i < 5; ++i) {
			printf("time rec %d: %lf \n", i, sim_time_rec[i]);
		}
		outcome[NUM_PLAYERS + 0] += sim_time_rec[0];
		outcome[NUM_PLAYERS + 1] += sim_time_rec[1];
		outcome[NUM_PLAYERS + 2] += sim_time_rec[2];
		outcome[NUM_PLAYERS + 3] += sim_time_rec[3];
		outcome[NUM_PLAYERS + 4] += sim_time_rec[4]; 
	}
	return outcome;
}


int main(int argc, char** argv) {
	printf("each network memory: %lf MB\ntotal memory used by the %d networks: %lf MB\nsize if using fully connected: %lf MB\n",
	 (double) (((nx * ny * (Size_cluster * Size_cluster) + (((nx - 1)*ny) + ((ny - 1)*nx)) + Ncluster_links) * sizeof(long double) * 2)/1000000.0),
	 popsize, (double) (((nx * ny * (Size_cluster * Size_cluster) + (((nx - 1)*ny) + ((ny - 1)*nx)) + Ncluster_links) * sizeof(long double) * popsize * 2)/1000000.0),
	 (double) (nx*nx*ny*ny*Size_cluster*Size_cluster*sizeof(long double)*popsize*2) / 1000000.0);
	int num_teams = (NUM_PLAYERS / PLAYERS_PER_TEAM);
	int pop_per_team = popsize / num_teams;

	clock_t start, end;
	clock_t threading_start, threading_join;
	int num_sims = Learning_time * Ngame * (popsize/NUM_PLAYERS);
	int sim_count = 0;
	float *sims_time_rec = malloc(sizeof(float) * num_sims * 5);
	double time_join_loss;
	double time_total = 0.0;
	double longest_sim;
	double round_time = 0.0;

  	/* NEURAL NETWORK INITIALIZATION */
	float high_score = -1000.0;
	/* printf("popsize: %d, MAX_NEURON: %d, MAX_LINKS: %d, Ngame: %d, Learning_time: %d\n", popsize, MAX_NEURON, MAX_LINKS, Ngame, Learning_time); */
	printf("CURRENT HIGH SCORE: %.2f\n", high_score);

	srand(time(NULL));

	assert(popsize % NUM_PLAYERS == 0);

	Species *animal_kingdom = malloc(sizeof(Species) * num_teams);
	for(int i = 0; i < num_teams; ++i) {
		if(i % 2 == 1) {
			//fully connected team
			animal_kingdom[i].type = 1;
			animal_kingdom[i].ga = Brain_GA_graph_init(pop_per_team,N_INPUT,N_OUTPUT,nx,ny,Size_cluster,Ncluster_links);
			animal_kingdom[i].schedule = malloc(sizeof(int)*((GA *) animal_kingdom[i].ga)->n);

		} else {
			animal_kingdom[i].type = 0;
			animal_kingdom[i].ga = GA_init(pop_per_team,N_INPUT,N_OUTPUT,FULL_MAX_NEURONS,FULL_MAX_LINKS);
			animal_kingdom[i].schedule = malloc(sizeof(int)*((Brain_GA *)animal_kingdom[i].ga)->n);
		}
	}

	FILE * ffit,* fsig;  
	ffit=fopen("fit.txt","w");
	fsig=fopen("sig.txt","w");


	//INITIALIZING THE SIMULATION HERE! 
	SimulationState sims[MAX_THREADS];
	for(int t=0;t<MAX_THREADS;++t) {
		InitSim(&sims[t]);
	}

	InitSim(&render_sim);

	pthread_t threads[MAX_THREADS];
	void *res;


	for (int li=0;li<Learning_time;++li) {  
		printf("round: %d out of %d, %.2f done\n",li,Learning_time,((float)li/(float)Learning_time));
		start = clock();

		//RESET FITNESS
		for(int ii=0;ii<num_teams;++ii) {
			if(animal_kingdom[ii].type == 0) {
				((GA *)(animal_kingdom[ii].ga))->Nimproved_fit=0;
				for (int i=0;i<((GA *)(animal_kingdom[ii].ga))->n;++i) {
					((GA *)(animal_kingdom[ii].ga))->copy_fit[i]=0;
				}
			} else {
				((Brain_GA *)(animal_kingdom[ii].ga))->Nimproved_fit=0;
				for (int i=0;i<((Brain_GA *)(animal_kingdom[ii].ga))->n;++i) {
					((Brain_GA *)(animal_kingdom[ii].ga))->copy_fit[i]=0;
				}
			}
		}

		for(int i=0;i<Ngame;++i){
			//GENERATE SCHEDULES, SET WHO PLAYED TO ZERO
			for(int p = 0; p < num_teams; ++p) {
				if(animal_kingdom[p].type == 0) {
					permute(animal_kingdom[p].schedule, ((GA *)(animal_kingdom[p].ga))->n);
				} else {
					permute(animal_kingdom[p].schedule, ((Brain_GA *)(animal_kingdom[p].ga))->n);
				}
			}
			int j = 0;
			while(j < popsize/NUM_PLAYERS) {
				if(DEBUG_THREADS) { printf("j: %d, popsize/NUM_PLAYERS:%d\n\t", j, popsize/NUM_PLAYERS); }
				threading_start = clock();
				longest_sim = -1.0;
				if(popsize/NUM_PLAYERS >= j + MAX_THREADS) {
					if(DEBUG_THREADS) { printf("creating threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<num_teams;++ii){
							for(int tb=0; tb<PLAYERS_PER_TEAM; ++tb) { //go through each team
								if(DEBUG_THREADS) { printf("\t\tsim %d player %d is getting species %d schedule place %d\n", t, (ii*PLAYERS_PER_TEAM)+tb, ii, (PLAYERS_PER_TEAM*(j+t))+tb); }
								if(sims[t].teams[ii].type == 0) {
									neuralnet_replace(((neuralnet *)(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br)), ((GA *)(animal_kingdom[ii].ga))->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);
								}
								else {
									brain_replace(((brain *)(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br)), ((Brain_GA *)(animal_kingdom[ii].ga))->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);
								}
							}
						}
						if(DEBUG_THREADS) { printf("\t\tready\n");}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					if(DEBUG_THREADS) { printf("joining threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
					  
						pthread_join(threads[t], &res);
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int current_team=0; current_team < num_teams; ++current_team) {
							for(int current_player = 0; current_player < PLAYERS_PER_TEAM; ++current_player) {
								if(DEBUG_THREADS) { printf("\t\tspecies %d at schedule place %d recieves results from player %d\n", current_team, (PLAYERS_PER_TEAM*(j+t))+current_player, current_player+(PLAYERS_PER_TEAM*current_team)); }
								if(sims[t].teams[current_team].type == 0) {
									((GA *) (animal_kingdom[current_team].ga))->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								} else {
									((Brain_GA *) (animal_kingdom[current_team].ga))->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								}

								if(((float *) res)[current_player] > high_score) {
	    							high_score = ((float *) res)[current_player];
									printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    						}
							}
						}
						if(DEBUG_TIME) { 
							if(((float *) res)[NUM_PLAYERS + 0] > longest_sim) {
								longest_sim = ((float *) res)[NUM_PLAYERS + 0];
							}
							for(int t_r = 0; t_r < 5; ++t_r) {
								sims_time_rec[(sim_count * 5) + t_r] = ((float *) res)[NUM_PLAYERS + t_r];
							}
						}
						printf("sim rec: %d\n", sim_count); 

						sim_count++;
						free(res);
					}
					j = j + MAX_THREADS;
				}
				else {  //case when games left is less than MAX_THREADS
					int leftover = (popsize/NUM_PLAYERS) - 1 - j;
					if(DEBUG_THREADS) { printf("leftovers: %d\n", leftover); }
					for(int t=0;t<leftover;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<num_teams;++ii){
							for(int tb=0; tb<PLAYERS_PER_TEAM; ++tb) { //go through each team
								if(DEBUG_THREADS) { printf("\t\tsim %d player %d is getting species %d schedule place %d\n", t, (ii*PLAYERS_PER_TEAM)+tb, ii, (PLAYERS_PER_TEAM*(j+t))+tb); }
								if(sims[t].teams[ii].type == 0) {
									neuralnet_replace( ((neuralnet *)(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br)), ((GA *)(animal_kingdom[ii].ga))->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);
								}
								else {
									brain_replace(((brain *)(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br)), ((Brain_GA *)(animal_kingdom[ii].ga))->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);
								}							}
						}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					for(int t=0;t<leftover;++t) {
						pthread_join(threads[t], &res);
						for(int current_team=0; current_team < num_teams; ++current_team) {
							for(int current_player = 0; current_player < PLAYERS_PER_TEAM; ++current_player) {
								if(DEBUG_THREADS) { printf("\t\tspecies %d at schedule place %d recieves results from player %d\n", current_team, (PLAYERS_PER_TEAM*(j+t))+current_player, current_player+(PLAYERS_PER_TEAM*current_team)); }
								if(sims[t].teams[current_team].type == 0) {
									((GA *) (animal_kingdom[current_team].ga))->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								} else {
									((Brain_GA *) (animal_kingdom[current_team].ga))->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								}								
								if(((float *) res)[current_player] > high_score) {
	    							high_score = ((float *) res)[current_player];
									printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    						}
							}
						}
						if(DEBUG_TIME) { 
							if(((float *) res)[NUM_PLAYERS + 0] > longest_sim) {
								longest_sim = ((float *) res)[NUM_PLAYERS + 0];
							}
							for(int t_r = 0; t_r < 5; ++t_r) {
								sims_time_rec[(sim_count * 5) + t_r] = ((float *) res)[NUM_PLAYERS + t_r];
							} 
						}
						printf("sim rec: %d\n", sim_count);
						sim_count++;
						free(res);
					}
					j = popsize/NUM_PLAYERS;
				}
				threading_join = clock();
				if(DEBUG_TIME) { time_join_loss += ((((double) threading_join - threading_start)/CLOCKS_PER_SEC) - longest_sim); }
			}
			//mutate_sigma(ga);
			//GA_mutate_weights(ga,CHANCE_MUTATION);
		}


		//SWAP FITNESS
		for(int ii=0;ii<num_teams;++ii) {
			if(animal_kingdom[ii].type == 0) {
				for (int i=0;i<((GA *) (animal_kingdom[ii].ga))->n;++i) {
					switch(((GA *) (animal_kingdom[ii].ga))->copy_fit[i] > ((GA *) (animal_kingdom[ii].ga))->fit_array[i]){
						case 1:
						++(((GA *) (animal_kingdom[ii].ga))->Nimproved_fit);
						break;
					}
				}
			} else {
				for (int i=0;i<((Brain_GA *)(animal_kingdom[ii].ga))->n;++i) {
					switch(((Brain_GA *)(animal_kingdom[ii].ga))->copy_fit[i] > ((Brain_GA *)(animal_kingdom[ii].ga))->fit_array[i]){
						case 1:
						++(((Brain_GA *)(animal_kingdom[ii].ga))->Nimproved_fit);
						break;
					}
				}
			}
		}

		for(int ii=0;ii<num_teams;++ii) {
			if(animal_kingdom[ii].type == 0) {
				memcpy(((GA *)(animal_kingdom[ii].ga))->fit_array,((GA *)(animal_kingdom[ii].ga))->copy_fit,((GA *)(animal_kingdom[ii].ga))->n*sizeof(float));
				tournament_selection(((GA *)(animal_kingdom[ii].ga)));
				mutate_sigma(((GA *)(animal_kingdom[ii].ga)));
				GA_mutate_weights(((GA *)(animal_kingdom[ii].ga)),CHANCE_MUTATION);
			} else {
				memcpy(((Brain_GA *)(animal_kingdom[ii].ga))->fit_array,((Brain_GA *)(animal_kingdom[ii].ga))->copy_fit,((Brain_GA *)(animal_kingdom[ii].ga))->n*sizeof(float));
				Brain_GA_tournament_selection(((Brain_GA *)(animal_kingdom[ii].ga)));
				Brain_GA_mutate_sigma(((Brain_GA *)(animal_kingdom[ii].ga)));
				Brain_GA_mutate_weights(((Brain_GA *)(animal_kingdom[ii].ga)),CHANCE_MUTATION);
				Brain_GA_mutate_table(((Brain_GA *)(animal_kingdom[ii].ga)),CHANCE_MUTATION);
				Brain_GA_out_fit(ffit,((Brain_GA *)(animal_kingdom[ii].ga)));  
				Brain_GA_out_sig(fsig,((Brain_GA *)(animal_kingdom[ii].ga)));
			}
		}
		end = clock();
		round_time = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("round %d took %lf secs\n", li, round_time);
		printf("estimated time left:\n\tsecs: %f\n\tmins: %f\n\thours: %f\n\tdays: %f\n", 
			(round_time*(Learning_time - li)), 
			(round_time*(Learning_time - li)/60),
			(round_time*(Learning_time - li)/(60*60)),
			(round_time*(Learning_time - li)/(60*60*24)));
		time_total += round_time;
	}


  	/*TRAINING DONE*/
	//NUM_PLAYERS + 0: total time of simulation
	//NUM_PLAYERS + 1: time performing sensor looks
	//NUM_PLAYERS + 2: time performing network pass
	//NUM_PLAYERS + 3: time performing shoot
	//NUM_PLAYERS + 4: time performing physics space step
	printf("Training done \n\n");
	printf("Total Training Time: %lf\n", time_total);
	if(DEBUG_TIME) { printf("Time lost from waiting for join: %lf\n", time_join_loss); }
	//num_sims
	double avg_sim_time = 0.0;
	double avg_sensor_time = 0.0;
	double avg_network_time = 0.0;
	double avg_shooting_time = 0.0;
	double avg_physics_time = 0.0;

	if(DEBUG_TIME) { 
		for(int i = 0; i < num_sims; ++i) {
			avg_sim_time 		+= sims_time_rec[(5*i)+0];
			avg_sensor_time 	+= sims_time_rec[(5*i)+1];
			avg_network_time 	+= sims_time_rec[(5*i)+2];
			avg_shooting_time 	+= sims_time_rec[(5*i)+3];
			avg_physics_time 	+= sims_time_rec[(5*i)+4];
		}

		avg_sim_time = avg_sim_time / (double) num_sims;
		avg_sensor_time = avg_sensor_time / (double) num_sims;
		avg_network_time = avg_network_time / (double) num_sims;
		avg_shooting_time = avg_shooting_time / (double) num_sims;
		avg_physics_time = avg_physics_time / (double) num_sims;

		printf("average simulation step time: %lf\n", avg_sim_time);
		printf("average sensor time: %lf\n", avg_sensor_time);
		printf("average network time: %lf\n", avg_network_time);
		printf("average shooting time: %lf\n", avg_shooting_time);
		printf("average physics time: %lf\n", avg_physics_time);
	}


	brain_write(((Brain_GA *)(animal_kingdom[1].ga))->pop[Brain_GA_max_fit(((Brain_GA *)(animal_kingdom[1].ga)))]);

	for(int current_team=0; current_team < num_teams; ++current_team) {
		for(int current_player = 0; current_player < PLAYERS_PER_TEAM; ++current_player) {
			if(animal_kingdom[current_team].type == 0) {
				render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br = neuralnet_init(N_INPUT,N_OUTPUT,FULL_MAX_NEURONS,FULL_MAX_LINKS);
				neuralnet_replace(((neuralnet *) (render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br)),((GA *)(animal_kingdom[current_team].ga))->pop[n_best(((GA *)(animal_kingdom[current_team].ga)), current_player)]);
			} else {
				render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br = brain_graph_init(N_INPUT,N_OUTPUT,nx,ny,Size_cluster,Ncluster_links);
				brain_replace(((brain *)(render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br)),((Brain_GA *)(animal_kingdom[current_team].ga))->pop[Brain_GA_n_best(((Brain_GA *)(animal_kingdom[current_team].ga)), current_player)]);
			}
		}
	}

	printf("agents placed into render_sim\n");
	ResetSim(&render_sim);
	/*clean the training equipment*/
	Brain_GA_out_fit(ffit,((Brain_GA *)(animal_kingdom[0].ga)));
	Brain_GA_out_sig(fsig,((Brain_GA *)(animal_kingdom[0].ga)));

	for(int current_team=0; current_team < num_teams; ++current_team) {
		if(animal_kingdom[current_team].type == 0) {
			GA_free(((GA *)(animal_kingdom[current_team].ga)));
		} else {
			Brain_GA_free(((Brain_GA *)(animal_kingdom[current_team].ga)));
		}
	}

	for(int t=0;t<MAX_THREADS;++t) {
		FreeSim(&sims[t]);
	}

	printf("sims freed");

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
