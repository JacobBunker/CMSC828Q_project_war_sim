#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <bsd/stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>


#define NUM_RIGID_BODIES 2
#define TOLERANCE 0.00001
// the tolerance should be something positive close to zero (ex. 0.00001)

//see https://www.toptal.com/game/video-game-physics-part-i-an-introduction-to-rigid-body-dynamics
//for resources on this system

typedef struct {
	float x;
	float y;
} Vector2;

struct Node {
	Vector2 data;
	struct Node* next;
};


typedef struct {
	Vector2 normal;
	int index;
	double distance;
} Edge;

typedef struct {
	float width;
	float height;
	float mass;
	float momentOfInertia;
	Vector2 vertices[4];
} BoxShape;

// Calculates the inertia of a box shape and stores it in the momentOfInertia variable.
void CalculateBoxInertia(BoxShape *boxShape) {
    float m = boxShape->mass;
    float w = boxShape->width;
    float h = boxShape->height;
    boxShape->momentOfInertia = m * (w * w + h * h) / 12;
}

// Two dimensional rigid body
typedef struct {
    Vector2 position;
    Vector2 linearVelocity;
    float angle;
    float angularVelocity;
    Vector2 force;
    float torque;
    BoxShape shape;
    int vertex_number;
} RigidBody;

RigidBody rigidBodies[NUM_RIGID_BODIES];

void PrintRigidBodies() {
    for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        printf("body[%i] p = (%.2f, %.2f), a = %.2f\n", i, rigidBody->position.x, rigidBody->position.y, rigidBody->angle);
		FILE * fp;
		fp = fopen ("/home/jack/828/war_sim/positions.txt","a");
 		fprintf (fp, "%i %.2f %.2f %.2f\n", i, rigidBody->position.x, rigidBody->position.y, rigidBody->angle);  
  	 	fclose (fp);
    }
}

// Initializes rigid bodies with random positions and angles and zero linear and angular velocities.
// They're all initialized with a box shape of random dimensions.
void InitializeRigidBodies() {
    for (int i = 0; i < NUM_RIGID_BODIES - 1; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        rigidBody->position = (Vector2){arc4random_uniform(50), arc4random_uniform(50)};
        rigidBody->angle = arc4random_uniform(360)/360.f * M_PI * 2;
        rigidBody->linearVelocity = (Vector2){0, 0};
        rigidBody->angularVelocity = 0;
        
        BoxShape shape;
        shape.mass = 10;
        shape.width = 1 + arc4random_uniform(2);
        shape.height = 1 + arc4random_uniform(2);

        //top right
        shape.vertices[0].x = shape.width/2.0;
        shape.vertices[0].y = shape.height/2.0;
        //top left
        shape.vertices[1].x = -shape.width/2.0;
        shape.vertices[1].y = shape.height/2.0;
        //bottom right
        shape.vertices[2].x = shape.width/2.0;
        shape.vertices[2].y = -shape.height/2.0;
        //bottom left
        shape.vertices[3].x = -shape.width/2.0;
        shape.vertices[3].y = -shape.height/2.0;
        rigidBody->vertex_number = 4;
        CalculateBoxInertia(&shape);
        rigidBody->shape = shape;
		FILE * fp;
		fp = fopen ("/home/jack/828/war_sim/positions.txt","a");
 		fprintf (fp, "b %.2f %.2f\n", shape.vertices[0].x, shape.vertices[0].y);  
 		fprintf (fp, "b %.2f %.2f\n", shape.vertices[1].x, shape.vertices[1].y);  
 		fprintf (fp, "b %.2f %.2f\n", shape.vertices[2].x, shape.vertices[2].y);  
 		fprintf (fp, "b %.2f %.2f\n", shape.vertices[3].x, shape.vertices[3].y);  
  	 	fclose (fp);
    }

    //create the ground
    RigidBody *rigidBody = &rigidBodies[NUM_RIGID_BODIES-1];
    rigidBody->position = (Vector2){0.0, -25.0};
    rigidBody->angle = 0.0;
    rigidBody->linearVelocity = (Vector2){0, 0};
    rigidBody->angularVelocity = 0;
    
    BoxShape shape;
    shape.mass = 1000000;
    shape.width = 1000;
    shape.height = 50;

    //top right
    shape.vertices[0].x = shape.width/2.0;
    shape.vertices[0].y = shape.height/2.0;
    //top left
    shape.vertices[1].x = -shape.width/2.0;
    shape.vertices[1].y = shape.height/2.0;
    //bottom right
    shape.vertices[2].x = shape.width/2.0;
    shape.vertices[2].y = -shape.height/2.0;
    //bottom left
    shape.vertices[3].x = -shape.width/2.0;
    shape.vertices[3].y = -shape.height/2.0;
    rigidBody->vertex_number = 4;
    CalculateBoxInertia(&shape);
    rigidBody->shape = shape;
}

// Applies a force at a point in the body, inducing some torque.
void ComputeForceAndTorque(RigidBody *rigidBody) {
    Vector2 f = (Vector2){0, -100};
    rigidBody->force = f;
    // r is the 'arm vector' that goes from the center of mass to the point of force application
    Vector2 r = (Vector2){rigidBody->shape.width / 2, rigidBody->shape.height / 2};
    rigidBody->torque = r.x * f.y - r.y * f.x;
}

void FreeList(struct Node* head)
{
   	struct Node* tmp;

   	while (head != NULL)
   	{
       	tmp = head;
       	head = head->next;
       	free(tmp);
    }
}


Vector2 VectorSub(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	return c;
}

Vector2 VectorNegate(Vector2 a) {
	Vector2 b;
	b.x = -a.x;
	b.y = -a.y;
	return b;
}

Vector2 VectorTripleProduct(Vector2 a, Vector2 b, Vector2 c) {
	Vector2 d;
	d.x = a.x * b.x * c.x;
	d.y = a.y * b.y * c.y;
	return d;
}

double VectorNorm(Vector2 a) {
	return sqrt((a.x*a.x) + (a.y*a.y));
}

Vector2 VectorDivideConstant(Vector2 a, double b) {
	Vector2 c;
	c.x = a.x / b;
	c.y = a.y / b;
	return c;
}

Vector2 VectorNormalize(Vector2 a){
	return VectorDivideConstant(a, VectorNorm(a));
}


Vector2 AffineTransform(Vector2 vertex, Vector2 position, float angle) {
	double cosA = cos((double) (angle*M_PI/180.0));
	double sinA = sin((double) (angle*M_PI/180.0));
	Vector2 output;
	output.x = ((vertex.x * cosA) - (vertex.y * sinA)) + position.x;
	output.y = ((vertex.x * sinA) + (vertex.y * cosA)) + position.y;
	return output;
}

Vector2 GetFarthestPointInDirection(Vector2 *vertices, int count, Vector2 d) {
	float highest = -FLT_MAX;
	Vector2 support = (Vector2) {0, 0};

	for (int i = 0; i < count; i++) {
		Vector2 v = vertices[i];
		float dot = v.x * d.x + v.y * d.y;

		if (dot > highest) {
			highest = dot;
			support = v;
		}
	}

	return support;
}

Vector2 Support(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 d) {
	//d is a vector direction
	//get points on edge of the shapes in opposite directions
	Vector2 zero;
	zero.x = 0.0; zero.y = 0.0;
	Vector2 p1 = GetFarthestPointInDirection((rigidBodyOne->shape).vertices, rigidBodyOne->vertex_number, AffineTransform(d, zero, rigidBodyOne->angle));
	Vector2 p2 = GetFarthestPointInDirection((rigidBodyTwo->shape).vertices, rigidBodyTwo->vertex_number, AffineTransform(VectorNegate(d), zero, rigidBodyTwo->angle));
	//perform Minkowski Difference
	Vector2 p3 = VectorSub(AffineTransform(p1, rigidBodyOne->position, rigidBodyOne->angle), AffineTransform(p2, rigidBodyTwo->position, rigidBodyTwo->angle));
	//p3 is now a point in the Minkowski space ont he edge of the Minkowski Difference
	return p3;
}


//methods to handle simplex
int GetEmpty(int list[3]) {
	if(list[0] < 0) {
		return 0;
	}
	if(list[1] < 0) {
		return 1;
	}
	if(list[2] < 0) {
		return 2;
	}
	return -1;
}

int GetLast(int list[3]) {
	if(list[0] > list[1]) {
		if(list[0] > list[2]) {
			return 0;
		}
		else {
			return 2;
		}
	} else {
		if(list[1] > list[2]) {
			return 1;
		}
		else {
			return 2;
		}
	}
}

//get the next to last addition
int GetB(int list[3], int a) {
	if(a == 0) {
		if(list[1] > list[2]) {
			return 1;
		}
		else {
			return 2;
		}
	}
	if(a == 1) {
		if(list[0] > list[2]) {
			return 0;
		}
		else {
			return 2;
		}
	}
	if(a == 2) {
		if(list[0] > list[1]) {
			return 0;
		}
		else {
			return 1;
		}
	}
	return -1;
}

//get the least recent addition
int GetC(int list[3], int a) {
	if(a == 0) {
		if(list[1] < list[2]) {
			return 1;
		}
		else {
			return 2;
		}
	}
	if(a == 1) {
		if(list[0] < list[2]) {
			return 0;
		}
		else {
			return 2;
		}
	}
	if(a == 2) {
		if(list[0] < list[1]) {
			return 0;
		}
		else {
			return 1;
		}
	}
	return -1;
}

int GetSize(int list[3]) {
	int i = 0;
	if(list[0] >= 0) {
		i++;
	}
	if(list[1] >= 0) {
		i++;
	}
	if(list[2] >= 0) {
		i++;
	}
	return i;
}

double DotProduct(Vector2 a, Vector2 b) {
	return (a.x*b.x) + (a.y*b.y);
}


int GJK(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 *simplex) {
	Vector2 a,b,c,ao,ab,ac,abPerp,acPerp;
	int t, last, b_i, c_i;

	Vector2 d;
	d.x = 0.0;
	d.y = 1.0; //arbitrary init currently
	//Vector2 simplex[3];
	int simplex_record[3];
	int counter = 1;
	simplex_record[0] = 0;
	simplex_record[1] = -1;
	simplex_record[2] = -1;
	simplex[0] = Support(rigidBodyOne, rigidBodyTwo, d);
	d = VectorNegate(d);
	while(1) {
		t = GetEmpty(simplex_record);
		simplex[t] = Support(rigidBodyOne, rigidBodyTwo, d);
		simplex_record[t] = counter;
		counter++;
		last = GetLast(simplex_record);
		a = simplex[last];
		if((DotProduct(a, d) <= 0) && (GetSize(simplex_record) == 3)) {
			//if the point added last was not past the origin in the direction of d
			//then the Minkowski Sum cannot possibly contain the origin since
			//the last point added is on the edge of the Minkowski Difference
			return 0;
		} else {
			//otherwise we need to determine if the origin is in the current simplex
			ao = VectorNegate(a);
			if(GetSize(simplex_record) == 3) {
				//triangle case
				//get b and c
				b_i = GetB(simplex_record, last);
				c_i = GetC(simplex_record, last);
				b = simplex[b_i];
				c = simplex[c_i];
				//compute the edges
				ab = VectorSub(b, a);
				ac = VectorSub(c, a);
				//compute the normals
				abPerp = VectorTripleProduct(ac, ab, ab);
				acPerp = VectorTripleProduct(ab, ac, ac);
				//is the origin in R4
				if(DotProduct(abPerp, ao) > 0) {
					//remove point C
					simplex_record[c_i] = -1;
					//set the new direction to abPerp
					d = abPerp;
				} else {
					//is the origin in R3
					if(DotProduct(acPerp, ao) > 0) {
						//remove point B
						simplex_record[b_i] = -1;
						//set the new direction to acPerp
						d = acPerp;
					} else {
						//otherwise we know it's in R5 so we can return true
						return 1;
					}
				}
			} else {
				//line segment case
				b_i = GetB(simplex_record, last);
				b = simplex[b_i];
				//compute AB
				ab = VectorSub(b, a);
				//get the perp to AB in the direction of the origin
				abPerp = VectorTripleProduct(ab, ao, ab);
				//set the direction to abPerp
				d = abPerp;
			}
		}
	}

	return 0;
}


Edge FindClosestEdge(struct Node *simplex_head) {
	Edge closest;
	closest.distance = DBL_MAX;
	Vector2 a, b, e, oa, n;
	double d;

	struct Node* current = simplex_head;
	struct Node* next;

	int last_switch = 0;
	int i = 1;

	while(last_switch == 0) {
		next = current->next;
		if(next == NULL) { //do wrap around
			next = simplex_head;
			last_switch = 1;
			i = 0;
		}
		//get current point and the next one
		a = current->data;
		b = next->data;
		//create the edge vector
		e = VectorSub(b, a);
		oa = a; //a - origin
		//get the vector from the edge towards the origin
		n = VectorTripleProduct(e, oa, e);
		//normalize the vector
		n = VectorNormalize(n);
		//calculate distance from origin to the edge
		d = DotProduct(n, a);
		//check for min distance
		if (d < closest.distance) {
			closest.distance = d;
			closest.normal = n;
			closest.index = i;
		}
		current = next;
		i++;
	}

	return closest;
}



Edge EPA(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo) {
	struct Node* head = NULL;
	struct Node* second = NULL;
	struct Node* third = NULL;

	Edge e;
	Vector2 p;
	double d;
	int i;

	Vector2 simplex[3];
	int gjk_result = GJK(rigidBodyOne, rigidBodyTwo, simplex);
	printf("GJK complete with result %d\n", gjk_result);
	if(gjk_result == 0) {
		//no collision, so report negative distance
		e.distance = -1.0;
		return e;
	}

	head = (struct Node *)malloc(sizeof(struct Node));
	second = (struct Node *)malloc(sizeof(struct Node));
	third = (struct Node *)malloc(sizeof(struct Node));
	head->data = simplex[0];
	head->next = second;
	second->data = simplex[1];
	second->next = third;
	third->data = simplex[2];
	third->next = NULL;

	while(1) {
		//obtain the edge closest to the origin on the Minkowski Difference
		e = FindClosestEdge(head);
		p = Support(rigidBodyOne, rigidBodyTwo, e.normal);
		//check distance from origin to the edge against
		//distance p along e.normal
		d = DotProduct(p, e.normal);
		if(d - e.distance < TOLERANCE) {
			//if the difference is less than the tolerance then we can assume
			//that we cannot expand the simplex any further, and have our solution
			e.distance = d; //re-using struct for convenience
			FreeList(head);
			return e;
		} else {
			printf("over tolerance with value: %.6f\n", (d - e.distance));
			//have not reached the edge of the Minkowski Difference
			// continue expanding by adding the new point to the simlpex
			//in between the points that made the closest edge
			struct Node* current = NULL;
			struct Node* previous = NULL;
			struct Node* temp = (struct Node *)malloc(sizeof(struct Node));

			current = head;
			i = 0;
			while(i < e.index) {
				previous = current;
				current = current->next;
				i++;
			}
			if(i == 0) { //betwen the last point and the first
				temp->data = head->data;
				temp->next = head->next;
				head->data = p;
				head->next = temp;
			} else { //between the current and previous
				temp->data = p;
				temp->next = current;
				previous->next = temp;
			}
		}
	}
}

void RunRigidBodySimulation() {
    float totalSimulationTime = 10; // The simulation will run for 10 seconds.
    float currentTime = 0; // This accumulates the time that has passed.
    float dt = 1; // Each step will take one second.
    
    InitializeRigidBodies();
    PrintRigidBodies();
    
    while (currentTime < totalSimulationTime) {
            
        for (int i = 0; i < NUM_RIGID_BODIES - 1; ++i) {
            RigidBody *rigidBody = &rigidBodies[i];
            ComputeForceAndTorque(rigidBody);
            Vector2 linearAcceleration = (Vector2){rigidBody->force.x / rigidBody->shape.mass, rigidBody->force.y / rigidBody->shape.mass};
            rigidBody->linearVelocity.x += linearAcceleration.x * dt;
            rigidBody->linearVelocity.y += linearAcceleration.y * dt;

            //brute force detect collisions (very inefficient- should change this!)
            for (int ii = 0; ii < NUM_RIGID_BODIES; ++ii) {
            	if(i != ii) { //can't collide with self
            		Edge e = EPA(rigidBody, &rigidBodies[ii]);
            		printf("EPA result between %d and %d: %0.4f\n", i, ii, e.distance);
            	}
			}

            rigidBody->position.x += rigidBody->linearVelocity.x * dt;
            rigidBody->position.y += rigidBody->linearVelocity.y * dt;
            float angularAcceleration = rigidBody->torque / rigidBody->shape.momentOfInertia;
            rigidBody->angularVelocity += angularAcceleration * dt;
            rigidBody->angle += rigidBody->angularVelocity * dt;
        }
        
        PrintRigidBodies();
        currentTime += dt;
    }
}

int main( int argc, const char* argv[] )
{
	printf("\nWAR GAME SIMULATION\n\n");
	RunRigidBodySimulation();
}
