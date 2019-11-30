#ifndef maze_simulator_H
#define maze_simulator_H

typedef struct {
  // change GA *ga;
        Brain_GA *ga;
	int team;
	int size;
	int *schedule;
	int *who_played;
} Species;

typedef struct {
	int object_type;
	int object_id;
	int team_id;
} ObjectInfo;

typedef struct {
	int id;
	cpGroup group;
	cpShapeFilter filter;
	cpVect spawn;
} Team;

typedef struct {
	int id;
	int team;
	brain *br;
	float score;
	float prev_look;
	float prev_shoot;
	float health;

	cpBody *body;
	cpShape *shape;
} Player;

typedef struct {
	cpSpace *space;

	cpVect **obstacles;
	cpBody **obstacle_bodies;
	cpShape **obstacle_shapes;

	ObjectInfo *object_infos;

	Team *teams;
	Player *players;

	int obstacle_number;
	int *obstacle_vertices_count;
	float *raw_obstacle_vertices;
	float force_multiplier;
	float look_multiplier;

	cpFloat timeStep;
	cpFloat arena_time_max;
	cpVect target;
	float t_s_dist;

	float time;
	int render;
	cpVect pos;
	cpVect vel;
	float *input;
	float *output;
	float scan_radius;

	cpSegmentQueryInfo scan_info;
	cpDataPointer scan_hit_data;
  	cpVect scan_target;
  	cpShape *scan_hit;
  	ObjectInfo hit_type;

  	float look_spread;
  	int look_number;

} SimulationState;

#endif

