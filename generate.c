#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define WALL_WIDTH 2.0
#define SEED 41
#define PI 3.14159265

typedef struct
{
	float x;
	float y;
} Vector2;

typedef struct
{
	int point_number;
	Vector2 **points;
} Obstacle;

typedef struct
{
	int spawn_count;
	int obstacle_count;
	float spawn_radius;
	Vector2 **spawns;
	Obstacle **obstacles;
} Map;

void GenerateWalls(float size, Map *map)
{
	int i = 0;
	// Define Left Wall
	map->obstacles[i] = malloc(sizeof(Obstacle));
	map->obstacles[i]->point_number = 4;
	map->obstacles[i]->points = malloc(4 * sizeof(char *));
	map->obstacles[i]->points[0] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[0]->x = 0.0;
	map->obstacles[i]->points[0]->y = 0.0;
	map->obstacles[i]->points[1] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[1]->x = WALL_WIDTH;
	map->obstacles[i]->points[1]->y = 0.0;
	map->obstacles[i]->points[2] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[2]->x = WALL_WIDTH;
	map->obstacles[i]->points[2]->y = size;
	map->obstacles[i]->points[3] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[3]->x = 0.0;
	map->obstacles[i]->points[3]->y = size;
	i++;

	//Define Bottom Wall
	map->obstacles[i] = malloc(sizeof(Obstacle));
	map->obstacles[i]->point_number = 4;
	map->obstacles[i]->points = malloc(4 * sizeof(char *));
	map->obstacles[i]->points[0] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[0]->x = 0.0;
	map->obstacles[i]->points[0]->y = 0.0;
	map->obstacles[i]->points[1] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[1]->x = size;
	map->obstacles[i]->points[1]->y = 0.0;
	map->obstacles[i]->points[2] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[2]->x = size;
	map->obstacles[i]->points[2]->y = WALL_WIDTH;
	map->obstacles[i]->points[3] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[3]->x = 0.0;
	map->obstacles[i]->points[3]->y = WALL_WIDTH;
	i++;

	//Define Right Wall
	map->obstacles[i] = malloc(sizeof(Obstacle));
	map->obstacles[i]->point_number = 4;
	map->obstacles[i]->points = malloc(4 * sizeof(char *));
	map->obstacles[i]->points[0] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[0]->x = size - WALL_WIDTH;
	map->obstacles[i]->points[0]->y = 0.0;
	map->obstacles[i]->points[1] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[1]->x = size;
	map->obstacles[i]->points[1]->y = 0.0;
	map->obstacles[i]->points[2] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[2]->x = size;
	map->obstacles[i]->points[2]->y = size;
	map->obstacles[i]->points[3] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[3]->x = size - WALL_WIDTH;
	map->obstacles[i]->points[3]->y = size;
	i++;

	//Define Top Wall
	map->obstacles[i] = malloc(sizeof(Obstacle));
	map->obstacles[i]->point_number = 4;
	map->obstacles[i]->points = malloc(4 * sizeof(char *));
	map->obstacles[i]->points[0] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[0]->x = 0.0;
	map->obstacles[i]->points[0]->y = size - WALL_WIDTH;
	map->obstacles[i]->points[1] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[1]->x = size;
	map->obstacles[i]->points[1]->y = size - WALL_WIDTH;
	map->obstacles[i]->points[2] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[2]->x = size;
	map->obstacles[i]->points[2]->y = size;
	map->obstacles[i]->points[3] = malloc(sizeof(Vector2));
	map->obstacles[i]->points[3]->x = 0.0;
	map->obstacles[i]->points[3]->y = size;
}

int touching_circles(float x1, float y1, float x2,
					 float y2, float r1, float r2)
{
	float distSq = (x1 - x2) * (x1 - x2) +
				 (y1 - y2) * (y1 - y2);
	float radSumSq = (r1 + r2) * (r1 + r2);
	return distSq < radSumSq;
}

int touching_wall_lower(float c, float r, float size){
	return (c-r < WALL_WIDTH);
}
int touching_wall_upper(float c, float r, float size){
	return (c+r > size-WALL_WIDTH);
}

void GenerateSpawns(int spawns, float spawn_radius, float size, Map *map){
	int i = 0;
	while (i < spawns)
	{
		float x = WALL_WIDTH + spawn_radius + (float)rand() / (float)(RAND_MAX / (size-2*(WALL_WIDTH+spawn_radius)));
		float y = WALL_WIDTH + spawn_radius + (float)rand() / (float)(RAND_MAX / (size-2*(WALL_WIDTH+spawn_radius)));
		int j = 0;
		int touching = 0;
		while(j<i && !touching){
			touching = touching || touching_circles(x,y, map->spawns[j]->x, map->spawns[j]->y, spawn_radius, spawn_radius);
			j+= 1;
		}
		if (!touching){
			map->spawns[i] = malloc(sizeof(Vector2));
			map->spawns[i]->x = x;
			map->spawns[i]->y = y;
			i += 1;
		}
	}
}

void GenerateObstacle(float size, float spawn_radius, int spawns, int i, Map *map){
	int found = 0;
	while(!found){
		float r = WALL_WIDTH/2 + (float)rand() / (float)(RAND_MAX / (WALL_WIDTH/2));
		float x = r + WALL_WIDTH + (float)rand() / (float)(RAND_MAX / (size-2*(WALL_WIDTH+r)));
		float y = r + WALL_WIDTH + (float)rand() / (float)(RAND_MAX / (size-2*(WALL_WIDTH+r)));
		int j = 0;
		int touching = 0;
		while(j<spawns && !touching){
			touching = touching || touching_circles(x,y, map->spawns[j]->x, map->spawns[j]->y, r, spawn_radius);
			j+= 1;
		}
		if (!touching){
			map->obstacles[i] = malloc(sizeof(Obstacle));
			int points = 3 + (int) rand()/ (int)(RAND_MAX/4);
			map->obstacles[i]->points = malloc(points * sizeof(char *));
			float alpha = 0.0;
			for(int k = 0; k < points; k++){
				map->obstacles[i]->points[k] = malloc(sizeof(Vector2));
				map->obstacles[i]->points[k]->x = x + r*cos(alpha);
				map->obstacles[i]->points[k]->y = y + r*sin(alpha);
				alpha = alpha + ((float)rand() / (float)(RAND_MAX / (2.0*PI-alpha- PI/12)));
			}
			map->obstacles[i]->point_number = points;
			found = 1;
		}
	}
}

void GenerateObstacles(float size, float spawn_radius, int density, int spawns, Map *map){
	for(int i = 4; i < density+4; i++){
		GenerateObstacle(size, spawn_radius, spawns, i, map);
		
	}
}



void GenerateMap(float size, int spawns, float spawn_radius, int density, Map *map)
{
	map->spawn_count = spawns;
	map->obstacle_count = density + 4;
	map->spawn_radius = spawn_radius;
	map->spawns = malloc(spawns * sizeof(char*));
	map->obstacles = malloc((4 + density) * sizeof(char *));
	srand(SEED);
	GenerateWalls(size, map);
	GenerateSpawns(spawns, spawn_radius, size, map);
	GenerateObstacles(size, spawn_radius, density, spawns, map);
	

}