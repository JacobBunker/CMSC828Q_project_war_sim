#ifndef maze_simulator_H
#define maze_simulator_H


typedef struct {
	cpSpace *space;

	cpVect **obstacles;
	cpBody **obstacle_bodies;
	cpShape **obstacle_shapes;

	cpBody **player_bodies;
	cpShape **player_shapes;

	int obstacle_number;
	int *obstacle_vertices_count;
	float *raw_obstacle_vertices;

	cpFloat timeStep;
	cpFloat arena_time_max;
	cpVect target;
	cpVect spawn;

	neuralnet **players;
	float time;
	int render;
	cpVect pos;
	cpVect vel;
	float *input;
	float *output;

} SimulationState;

#endif