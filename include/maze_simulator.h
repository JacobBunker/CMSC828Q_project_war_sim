#ifndef maze_simulator_H
#define maze_simulator_H

typedef struct {
	int object_type;
	int object_id;
	int team_id;
} ObjectInfo;

typedef struct {
	cpSpace *space;

	cpVect **obstacles;
	cpBody **obstacle_bodies;
	cpShape **obstacle_shapes;

	cpBody **player_bodies;
	cpShape **player_shapes;

	ObjectInfo *object_infos;

	int obstacle_number;
	int *obstacle_vertices_count;
	float *raw_obstacle_vertices;
	float force_multiplier;
	float look_multiplier;

	cpFloat timeStep;
	cpFloat arena_time_max;
	cpVect target;
	cpVect spawn_one, spawn_two;
	float t_s_dist;

	neuralnet **players;
	float *player_scores;
	float *player_previous_look;
	float *player_previous_shoot;
	float *player_health;
	float time;
	int render;
	cpVect pos;
	cpVect vel;
	float *input;
	float *output;
	float scan_radius;

	cpShapeFilter team_one_filter;
	cpShapeFilter team_two_filter;
	cpSegmentQueryInfo scan_info;
	cpDataPointer scan_hit_data;
  	cpVect scan_target;
  	cpShape *scan_hit;
  	ObjectInfo hit_type;

  	float look_spread;
  	int look_number;

} SimulationState;

#endif