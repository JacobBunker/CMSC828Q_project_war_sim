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
#define MAX_THREADS 10


#define DEBUG_PRINT 0
#define DEBUG_THREADS 0
#define DEBUG_RES 0

SimulationState render_sim;

char *name;
char buffer[64];
cpVect spawns[NUM_PLAYERS/PLAYERS_PER_TEAM];
float team_colors[2][3];

int popsize=5*NUM_PLAYERS,
  nx=10,
  ny=10,
  Size_cluster=19,
  Ncluster_links=200,
  Ngame=5,
  Learning_time=5; 




void InitSim(SimulationState *sim, int player_net_param) {
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
		sim->players[i].br = brain_graph_init(N_INPUT, N_OUTPUT, nx,ny,Size_cluster,Ncluster_links );
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
	  brain_free(sim->players[i].br);
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
		brain_reset_act(sim->players[i].br);
	}
}

void StepSim(SimulationState *sim, int graphics_on) {
	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		if(sim->players[current_player].health > 0.0) {
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

			brain_forward_pass(sim->players[current_player].br,sim->input,sim->output);
			//printf("outputs: %.2f, %.2f\n", output[0], output[1]);
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

			cpBodySetForce(sim->players[current_player].body, sim->vel);
		}
		else {
			cpBodySetPosition(sim->players[current_player].body, cpv(-100, -100));
		}
			
	}
	cpSpaceStep(sim->space, sim->timeStep);
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 0.0f);
	//glBegin(GL_LINES);   
	if(render_sim.time < render_sim.arena_time_max) {
	
		StepSim(&render_sim, 1);
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
		// Clean up our ARENA objects
		//set up the render_sim arena, get best agents
		ResetSim(&render_sim);
		//FreeSim(&render_sim);
		//render_sim.render = 0;
	}  
}

void* SimulationThread(void *arg) {
	SimulationState *sim = (SimulationState *) arg;
	float *outcome = malloc(sizeof(float)*NUM_PLAYERS);

	for(int p_i = 0; p_i < NUM_PLAYERS; ++p_i) {
		outcome[p_i] = 0.0;
	}

	ResetSim(sim);

	for(cpFloat time = 0; time < sim->arena_time_max; time += sim->timeStep){
		StepSim(sim, 0);
	}

	for(int current_player=0;current_player<NUM_PLAYERS;++current_player) {
		//sim->pos = cpBodyGetPosition(sim->player_bodies[current_player]);
		outcome[current_player] = sim->players[current_player].score; 
		//+ (sim->t_s_dist - sqrt(((sim->pos.x - sim->target.x)*(sim->pos.x - sim->target.x)) + ((sim->pos.y - sim->target.y)*(sim->pos.y - sim->target.y))));
	}

	return outcome;
}


int main(int argc, char** argv) {

	int num_teams = (NUM_PLAYERS / PLAYERS_PER_TEAM);
	int pop_per_team = popsize / num_teams;

  	/* NEURAL NETWORK INITIALIZATION */
	float high_score = -1000.0;
	/* printf("popsize: %d, MAX_NEURON: %d, MAX_LINKS: %d, Ngame: %d, Learning_time: %d\n", popsize, MAX_NEURON, MAX_LINKS, Ngame, Learning_time); */
	printf("CURRENT HIGH SCORE: %.2f\n", high_score);

	srand(time(NULL));

	assert(popsize % NUM_PLAYERS == 0);

	Species *animal_kingdom = malloc(sizeof(Species) * num_teams);
	for(int i = 0; i < num_teams; ++i) {
	  animal_kingdom[i].ga = Brain_GA_graph_init(pop_per_team,N_INPUT,N_OUTPUT,nx,ny,Size_cluster,Ncluster_links);
		animal_kingdom[i].schedule = malloc(sizeof(int)*animal_kingdom[i].ga->n);
	}

	FILE * ffit,* fsig;  
	ffit=fopen("fit.txt","w");
	fsig=fopen("sig.txt","w");

	

	//INITIALIZING THE SIMULATION HERE! 
	SimulationState sims[MAX_THREADS];
	for(int t=0;t<MAX_THREADS;++t) {
		InitSim(&sims[t], animal_kingdom[0].ga->pop[0]->br_size);
	}

	InitSim(&render_sim, animal_kingdom[0].ga->pop[0]->br_size);

	pthread_t threads[MAX_THREADS];
	void *res;

	for (int li=0;li<Learning_time;++li) {  
		printf("round: %d out of %d, %.2f done\n",li,Learning_time,((float)li/(float)Learning_time));

		//RESET FITNESS
		for(int ii=0;ii<num_teams;++ii) {
			animal_kingdom[ii].ga->Nimproved_fit=0;
			for (int i=0;i<animal_kingdom[ii].ga->n;++i) {
				animal_kingdom[ii].ga->copy_fit[i]=0;
			}
		}

		for(int i=0;i<Ngame;++i){
			//GENERATE SCHEDULES, SET WHO PLAYED TO ZERO
			for(int p = 0; p < num_teams; ++p) {
				permute(animal_kingdom[p].schedule, animal_kingdom[p].ga->n);
			}
			int j = 0;
			while(j < popsize/NUM_PLAYERS) {
				if(DEBUG_THREADS) { printf("j: %d, popsize/NUM_PLAYERS:%d\n\t", j, popsize/NUM_PLAYERS); }
				if(popsize/NUM_PLAYERS >= j + MAX_THREADS) {
					if(DEBUG_THREADS) { printf("creating threads...\n");}
					for(int t=0;t<MAX_THREADS;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<num_teams;++ii){
							for(int tb=0; tb<PLAYERS_PER_TEAM; ++tb) { //go through each team
								if(DEBUG_THREADS) { printf("\t\tsim %d player %d is getting species %d schedule place %d\n", t, (ii*PLAYERS_PER_TEAM)+tb, ii, (PLAYERS_PER_TEAM*(j+t))+tb); }

								brain_replace(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br, animal_kingdom[ii].ga->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);

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
								animal_kingdom[current_team].ga->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								if(((float *) res)[current_player] > high_score) {
	    							high_score = ((float *) res)[current_player];
									printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    						}
							}
						}
						free(res);
					}
					j = j + MAX_THREADS;
				}
				else {  //case when games left is less than MAX_THREADS
					int leftover = (popsize/NUM_PLAYERS) - 1 - j;
					printf("leftovers: %d\n", leftover);
					for(int t=0;t<leftover;++t) {
						if(DEBUG_THREADS) { printf("\tthread%d\n",t);}
						for(int ii=0;ii<num_teams;++ii){
							for(int tb=0; tb<PLAYERS_PER_TEAM; ++tb) { //go through each team


								if(DEBUG_THREADS) { printf("\t\tsim %d player %d is getting species %d schedule place %d\n", t, (ii*PLAYERS_PER_TEAM)+tb, ii, (PLAYERS_PER_TEAM*(j+t))+tb); }

								brain_replace(sims[t].players[(ii*PLAYERS_PER_TEAM)+tb].br, animal_kingdom[ii].ga->pop[animal_kingdom[ii].schedule[(PLAYERS_PER_TEAM*(j+t))+tb]]);

							}
						}
						pthread_create(&threads[t], NULL, SimulationThread, &sims[t]);
					}
					for(int t=0;t<leftover;++t) {
						pthread_join(threads[t], &res);
						for(int current_team=0; current_team < num_teams; ++current_team) {
							for(int current_player = 0; current_player < PLAYERS_PER_TEAM; ++current_player) {
								if(DEBUG_THREADS) { printf("\t\tspecies %d at schedule place %d recieves results from player %d\n", current_team, (PLAYERS_PER_TEAM*(j+t))+current_player, current_player+(PLAYERS_PER_TEAM*current_team)); }
								animal_kingdom[current_team].ga->copy_fit[animal_kingdom[current_team].schedule[(PLAYERS_PER_TEAM*(j+t))+current_player]] += ((float *) res)[current_player+(PLAYERS_PER_TEAM*current_team)];
								if(((float *) res)[current_player] > high_score) {
	    							high_score = ((float *) res)[current_player];
									printf("CURRENT HIGH SCORE: %.2f\n", high_score);
	    						}
							}
						}
						free(res);
					}
					j = popsize/NUM_PLAYERS;
				}
			}
			//mutate_sigma(ga);
			//GA_mutate_weights(ga,CHANCE_MUTATION);
		}


		//SWAP FITNESS
		for(int ii=0;ii<num_teams;++ii) {
			for (int i=0;i<animal_kingdom[ii].ga->n;++i) {
				switch(animal_kingdom[ii].ga->copy_fit[i] > animal_kingdom[ii].ga->fit_array[i]){
					case 1:
					++(animal_kingdom[ii].ga->Nimproved_fit);
					break;
				}
			}
		}

		for(int ii=0;ii<num_teams;++ii) {
			memcpy(animal_kingdom[ii].ga->fit_array,animal_kingdom[ii].ga->copy_fit,animal_kingdom[ii].ga->n*sizeof(float));
			Brain_GA_tournament_selection(animal_kingdom[ii].ga);
			Brain_GA_mutate_sigma(animal_kingdom[ii].ga);
			Brain_GA_mutate_weights(animal_kingdom[ii].ga,CHANCE_MUTATION);
		        Brain_GA_mutate_table(animal_kingdom[ii].ga,CHANCE_MUTATION);
			
		}
		Brain_GA_out_fit(ffit,animal_kingdom[0].ga);  
		Brain_GA_out_sig(fsig,animal_kingdom[0].ga);
	}


  	/*TRAINING DONE*/

	printf("Training done \n\n");
	brain_write(animal_kingdom[0].ga->pop[Brain_GA_max_fit(animal_kingdom[0].ga)]);

	for(int current_team=0; current_team < num_teams; ++current_team) {
		for(int current_player = 0; current_player < PLAYERS_PER_TEAM; ++current_player) {
		  render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br = brain_graph_init(N_INPUT,N_OUTPUT,nx,ny,Size_cluster,Ncluster_links);
			brain_replace(render_sim.players[(current_team*PLAYERS_PER_TEAM) + current_player].br,animal_kingdom[current_team].ga->pop[Brain_GA_n_best(animal_kingdom[current_team].ga, current_player)]);
		}
	}

	ResetSim(&render_sim);
	/*clean the training equipment*/
	Brain_GA_out_fit(ffit,animal_kingdom[0].ga);
	Brain_GA_out_sig(fsig,animal_kingdom[0].ga);

	for(int current_team=0; current_team < num_teams; ++current_team) {
		Brain_GA_free(animal_kingdom[current_team].ga);
	}

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
