#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>


#define NUM_PARTICLES 1

//see https://www.toptal.com/game/video-game-physics-part-i-an-introduction-to-rigid-body-dynamics
//for resources on this system

typedef struct {
	float x;
	float y;
} Vector2;

typedef struct {
	Vector2 position;
	Vector2 velocity;
	float mass;
} Particle;

Particle particles[NUM_PARTICLES];

void PrintParticles() {
	for (int i = 0; i < NUM_PARTICLES; ++i) {
		Particle *particle = &particles[i];
		printf("particle[%i] (%.2f, %.2f)\n", i, particle->position.x, particle->position.y);
	}
}


//create particles with random positions, zero velocities, and 1 kg mass
void InitializeParticles() {
	for (int i = 0; i < NUM_PARTICLES; ++i) {
		particles[i].position = (Vector2){arc4random_uniform(50), arc4random_uniform(50)};
		particles[i].velocity = (Vector2){0, 0};
		particles[i].mass = 1;
	}
}

//just apply gravity (mass times 9.81 m/s^2) to each particle
Vector2 ComputeForce(Particle *particle) {
	return (Vector2){0, particle->mass * -9.81};
}

void RunSimulation() {
	float totalSimulationTime = 10;
	float currentTime = 0.0;
	float dt = 1; //physics will be calculated every second

	InitializeParticles();
	PrintParticles();

	while (currentTime < totalSimulationTime) {
		//because we have no realtime constraints we will be going full blast
		//without sleeping or trying to match the framerate

		for (int i = 0; i < NUM_PARTICLES; ++i) {
			Particle *particle = &particles[i];
			Vector2 force = ComputeForce(particle);
			Vector2 acceleration = (Vector2){force.x / particle->mass, force.y / particle->mass};

			//important that velocity is calculated before position because otherwise
			//the error factor will tend towards increasing speeds rather than decreasing,
			//leading to explosive behaviors. learn more about what's happening here: 
			//https://web.archive.org/web/20190405191806/https://gafferongames.com/post/integration_basics/
			particle->velocity.x += acceleration.x * dt;
			particle->velocity.y += acceleration.y * dt;
			particle->position.x += particle->velocity.x * dt;
			particle->position.y += particle->velocity.y * dt;
		}

		PrintParticles();
		currentTime += dt;
	}
}

int main( int argc, const char* argv[] )
{
	printf("\nWAR GAME SIMULATION\n\n");
	RunSimulation();
}
