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

//collision stuff 
#define MAX_SI_ITERATIONS 30
#define SI_TOLERANCE 0.001
#define SI_BIAS_FACTOR 0.0001


// the tolerance should be something positive close to zero (ex. 0.00001)

//see https://www.toptal.com/game/video-game-physics-part-i-an-introduction-to-rigid-body-dynamics
//for resources on this system

typedef struct {
	float x;
	float y;
} Vector2;

typedef struct {
	float a;
	float b;
	float c;
	float d;
} Vector4;

typedef struct {
	Vector2 a;
	Vector2 b;
	Vector2 c;
	Vector2 d;
} Vector4x2;

typedef struct {
	Vector4 a;
	Vector4 b;
	Vector4 c;
	Vector4 d;
} Vector4x4;

typedef struct {
	Vector2 r_a;
	Vector2 n;
} Constraint;

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
Constraint constraints[NUM_RIGID_BODIES];

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
        shape.vertices[3].x = shape.width/2.0;
        shape.vertices[3].y = -shape.height/2.0;
        //bottom left
        shape.vertices[2].x = -shape.width/2.0;
        shape.vertices[2].y = -shape.height/2.0;
        rigidBody->vertex_number = 4;
        CalculateBoxInertia(&shape);
        rigidBody->shape = shape;
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
    shape.vertices[3].x = shape.width/2.0;
    shape.vertices[3].y = -shape.height/2.0;
    //bottom left
    shape.vertices[2].x = -shape.width/2.0;
    shape.vertices[2].y = -shape.height/2.0;
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

int PrintList(struct Node* head) {
	int i = 0;
	struct Node* current;
	current = head;
	while(current != NULL) {
		printf("Node %d: %.2f %.2f\n", i, current->data.x, current->data.y);
		current = current->next;
		i++;
	}
	return i;
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

Vector2 VectorAdd(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	return c;
}

Vector2 VectorSub(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	//printf("vectorsub a: %.2f, %.2f\n \tb: %.2f, %.2f\n\tc: %.2f, %.2f\n", a.x, a.y, b.x, b.y, c.x, c.y);
	return c;
}

Vector2 VectorNegate(Vector2 a) {
	Vector2 b;
	b.x = -a.x;
	b.y = -a.y;
	return b;
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

Vector2 VectorMultiplyConstant(Vector2 a, double b) {
	Vector2 c;
	c.x = a.x * b;
	c.y = a.y * b;
	return c;
}

Vector2 Vec2Vec2Multiply(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x * b.x;
	c.y = a.y * b.y;
	return c;
}

Vector2 VectorNormalize(Vector2 a){
	return VectorDivideConstant(a, VectorNorm(a));
}


//Note that the following triple product expansion is used:
//(A x B) x C = B(C.dot(A)) – A(C.dot(B)) to evaluate the triple product.
Vector2 VectorTripleProduct(Vector2 a, Vector2 b, Vector2 c) {
	return VectorSub(VectorMultiplyConstant(b, DotProduct(c, a)), VectorMultiplyConstant(a, DotProduct(c, b)));;
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

Vector2 Support(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 d, Vector2 *simplex_shape_points) {
	//d is a vector direction
	//get points on edge of the shapes in opposite directions
	Vector2 zero;
	zero.x = 0.0; zero.y = 0.0;

	printf("Support Function input d: %.2f, %.2f\n", d.x, d.y);
	Vector2 p1 = GetFarthestPointInDirection((rigidBodyOne->shape).vertices, rigidBodyOne->vertex_number, AffineTransform(d, zero, rigidBodyOne->angle));
	printf("Support p1: %.2f, %.2f\n", p1.x, p1.y);
	Vector2 p2 = GetFarthestPointInDirection((rigidBodyTwo->shape).vertices, rigidBodyTwo->vertex_number, AffineTransform(VectorNegate(d), zero, rigidBodyTwo->angle));
	printf("Support p2: %.2f, %.2f\n", p2.x, p2.y);
	//perform Minkowski Difference
	Vector2 p3 = VectorSub(AffineTransform(p1, rigidBodyOne->position, rigidBodyOne->angle), AffineTransform(p2, rigidBodyTwo->position, rigidBodyTwo->angle));
	//p3 is now a point in the Minkowski space ont he edge of the Minkowski Difference
	printf("Support p3: %.2f, %.2f\n", p3.x, p3.y);
	simplex_shape_points[0] = p1;
	simplex_shape_points[1] = p2;
	return p3;
}



int GJK(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 *simplex_output, Vector2 *closest_points) {
	Vector2 a,b,c,ao,ab,ac,abPerp,acPerp;
	double lambda_one, lambda_two;
	int t, last, b_i, c_i;

	Vector2 d;
	d.x = 1.0;
	d.y = -1.0; //arbitrary init currently
	Vector2 simplex[3];
	Vector2 simplex_shape_points[3][2];
	int simplex_record[3];
	int counter = 1;
	simplex_record[0] = 0;
	simplex_record[1] = -1;
	simplex_record[2] = -1;
	simplex[0] = Support(rigidBodyOne, rigidBodyTwo, d, simplex_shape_points[0]);
	d = VectorNegate(d);
	while(1) {
		t = GetEmpty(simplex_record);
		simplex[t] = Support(rigidBodyOne, rigidBodyTwo, d, simplex_shape_points[t]);
		simplex_record[t] = counter;
		counter++;
		last = GetLast(simplex_record);
		a = simplex[last];
		printf("last: %.2f, %.2f\n", a.x, a.y);
		printf("proj: %f\n", DotProduct(a, d));
		if((DotProduct(a, d) <= 0)) {
			//if the point added last was not past the origin in the direction of d
			//then the Minkowski Sum cannot possibly contain the origin since
			//the last point added is on the edge of the Minkowski Difference
			return 0;
		} else {
			//otherwise we need to determine if the origin is in the current simplex
			ao = VectorNegate(a);
			printf("ao: %.2f, %.2f\n", ao.x, ao.y);

			if(GetSize(simplex_record) == 3) {
				//triangle case
				//get b and c
				b_i = GetB(simplex_record, last);
				c_i = GetC(simplex_record, last);
				b = simplex[b_i];
				c = simplex[c_i];
				printf("a: %.2f, %.2f\n", a.x, a.y);
				printf("b: %.2f, %.2f\n", b.x, b.y);
				printf("c: %.2f, %.2f\n", c.x, c.y);
				//compute the edges
				ab = VectorSub(b, a);
				printf("ab: %.2f, %.2f\n", ab.x, ab.y);
				ac = VectorSub(c, a);
				printf("ac: %.2f, %.2f\n", ac.x, ac.y);
				//compute the normals
				abPerp = VectorTripleProduct(ac, ab, ab);
				acPerp = VectorTripleProduct(ab, ac, ac);
				printf("abPerp: %.2f, %.2f\n", abPerp.x, abPerp.y);
				printf("acPerp: %.2f, %.2f\n", acPerp.x, acPerp.y);

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

						//calculate closest points
						printf("calculating closest points:\n\tpoints corresponding to minkowski point A (%.2f, %.2f):\n\t\tS1:%.2f,%.2f\n\t\tS2:%.2f,%.2f\n\tpoints corresponding to minkowski point B (%.2f, %.2f):\n\t\tS1:%.2f,%.2f\n\t\tS2:%.2f,%.2f\n", 
								a.x, a.y, simplex_shape_points[last][0].x, simplex_shape_points[last][0].y, simplex_shape_points[last][1].x, simplex_shape_points[last][1].y,
								b.x, b.y, simplex_shape_points[b_i][0].x, simplex_shape_points[b_i][0].y, simplex_shape_points[b_i][1].x, simplex_shape_points[b_i][1].y);

						//ab is equal to L, so we will use that
						//λ2 = -L · A / (L · L)
						//λ1 = 1 – λ2
						lambda_two = -(DotProduct(ab, a) / DotProduct(ab, ab));
						lambda_one = 1 - lambda_two;

						if(lambda_one < 0) {
							closest_points[0] = simplex_shape_points[b_i][0];
							closest_points[1] = simplex_shape_points[b_i][1];
						} else if (lambda_two < 0) {
							closest_points[0] = simplex_shape_points[last][0];
							closest_points[1] = simplex_shape_points[last][1];
						}
						else {
							closest_points[0] = VectorAdd(VectorMultiplyConstant(simplex_shape_points[last][0], lambda_one), VectorMultiplyConstant(simplex_shape_points[b_i][0], lambda_two));
							closest_points[1] = VectorAdd(VectorMultiplyConstant(simplex_shape_points[last][1], lambda_one), VectorMultiplyConstant(simplex_shape_points[b_i][1], lambda_two));
						}

						simplex_output[0] = c;
						simplex_output[1] = b;
						simplex_output[2] = a;
						return 1;
					}
				}
			} else {
				//line segment case
				b_i = GetB(simplex_record, last);
				b = simplex[b_i];
				printf("b: %.2f, %.2f\n", b.x, b.y);
				printf("a: %.2f, %.2f\n", a.x, a.y);
				//compute AB
				ab = VectorSub(a, b); //last minus b
				printf("ab: %.2f, %.2f\n", ab.x, ab.y);
				//get the perp to AB in the direction of the origin
				abPerp = VectorTripleProduct(ab, ao, ab);
				printf("abPerp: %.2f, %.2f\n", abPerp.x, abPerp.y);
				Vector2 abPerpNorm = VectorNormalize(abPerp);
				printf("abPerpNorm: %.5f, %.5f\n", abPerpNorm.x, abPerpNorm.y);
				//set the direction to abPerp
				d = abPerp;
			}
		}
	}

	return 0;
}


Edge FindClosestEdge(struct Node *simplex_head, int winding) {
	Edge closest;
	closest.distance = DBL_MAX;
	Vector2 a, b, e, n;
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
		printf("a: %d: %.2f %.2f\n", i-1, a.x, a.y);
		b = next->data;
		printf("b: %d: %.2f %.2f\n", i, b.x, b.y);
		//create the edge vector
		e = VectorSub(b, a);
		printf("e: %.2f %.2f\n", e.x, e.y);

		//oa = a; //a - origin
		//printf("oa: %.2f %.2f\n", a.x, a.y);

		//per-product to get normal of the edge
		if(winding == 0) {
			//clockwise polygon
			n.x = e.y;
			n.y = -e.x;
		} else {
			//counter-clockwise polygon
			n.x = -e.y;
			n.y = e.x;
		}
		//get the vector from the edge towards the origin
		//n = VectorTripleProduct(e, oa, e);

		printf("n1: %.2f, %.2f\n", n.x, n.y);

		//normalize the vector
		n = VectorNormalize(n);
		printf("n2: %.2f, %.2f\n", n.x, n.y);
		//calculate distance from origin to the edge
		d = DotProduct(n, a);
		printf("d: %.4f\n", d);
		//check for min distance
		if (d < closest.distance) {
			closest.distance = d;
			closest.normal = n;
			closest.index = i;
		}
		current = next;
		i++;
	}

	printf("closest:\n\tnorm: %.2f, %.2f\n\tdistance: %.4f\n\tindex: %d\n", closest.normal.x, closest.normal.y, closest.distance, closest.index);

	return closest;
}



Edge EPA(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 *closest_points) {
	struct Node* head = NULL;
	struct Node* second = NULL;
	struct Node* third = NULL;

	Edge e;
	Vector2 p;
	double d;
	int i;
	int winding;

	Vector2 simplex_shape_points_temp[2];
	Vector2 simplex[3];
	int gjk_result = GJK(rigidBodyOne, rigidBodyTwo, simplex, closest_points);
	printf("GJK complete with result:\n\tcollision: %d\n\tclosest point A: %.2f, %.2f\n\tclosest point B: %.2f, %.2f\n", 
				gjk_result, closest_points[0].x, closest_points[0].y, closest_points[1].x, closest_points[1].y);
	printf("simplex:\n \t1: %.2f, %.2f\n\t2: %.2f, %.2f\n\t3: %.2f, %.2f\n", simplex[0].x, simplex[0].y, simplex[1].x, simplex[1].y, simplex[2].x, simplex[2].y);
	
	if(gjk_result == 0) {
		//no collision, so report negative distance
		e.distance = -1.0;
		return e;
	}

	//detect winding with cross product
	//((xi - xi-1),(yi - yi-1)) x ((xi+1 - xi),(yi+1 - yi)) 
	// = (xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
	double cross_product = ((simplex[0].x - simplex[2].x) * (simplex[1].y - simplex[0].y)) -
							 ((simplex[0].y - simplex[2].y) * (simplex[1].x - simplex[0].x));
	printf("cross product of simplex: %f\n", cross_product);
	if(cross_product > 0.0) {
		//counter-clockwise polygon
		winding = 1;
	} else {
		//clockwise polygon
		winding = 0;
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

	printf("first list:\n");
	PrintList(head);

	while(1) {
		//obtain the edge closest to the origin on the Minkowski Difference
		e = FindClosestEdge(head, winding);
		p = Support(rigidBodyOne, rigidBodyTwo, e.normal, simplex_shape_points_temp);
		//check distance from origin to the edge against
		//distance p along e.normal
		d = DotProduct(p, e.normal);
		printf("current distance: %f\n", d);
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
			printf("new list:\n");
			PrintList(head);
		}
	}
}




void PrintRigidBodies(float time, FILE *fp) {
	Vector2 a,b,c,d;
	printf("CURRENT TIME: %0.2f\n", time);
    for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        //printf("body[%i] p = (%.2f, %.2f), a = %.2f\n", i, rigidBody->position.x, rigidBody->position.y, rigidBody->angle);
        a = AffineTransform(rigidBody->shape.vertices[0], rigidBody->position, rigidBody->angle);
        b = AffineTransform(rigidBody->shape.vertices[1], rigidBody->position, rigidBody->angle);
        c = AffineTransform(rigidBody->shape.vertices[3], rigidBody->position, rigidBody->angle);
        d = AffineTransform(rigidBody->shape.vertices[2], rigidBody->position, rigidBody->angle);
        fprintf(fp, "%.2f Object %d %f %f %f %f %f %f %f %f\n", time, i, a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y);
    }
    //printf("\n");
}


void RunRigidBodySimulation() {

	FILE *fp;
	fp = fopen("/home/jbunker/c_workspace/war_games/sim_results.txt", "w");

    float totalSimulationTime = 20; // The simulation will run for 10 seconds.
    float currentTime = 0; // This accumulates the time that has passed.
    float dt = 0.25; // Each step will take one second.

    int constraint_count;
    int constraint_iterator;
	double max_impulse;
    
    InitializeRigidBodies();
    PrintRigidBodies(currentTime, fp);
    
    while (currentTime < totalSimulationTime) {

        for (int i = 0; i < NUM_RIGID_BODIES - 1; ++i) {
            RigidBody *rigidBody = &rigidBodies[i];
            ComputeForceAndTorque(rigidBody);
            Vector2 linearAcceleration = (Vector2){rigidBody->force.x / rigidBody->shape.mass, rigidBody->force.y / rigidBody->shape.mass};
            rigidBody->linearVelocity.x += linearAcceleration.x * dt;
            rigidBody->linearVelocity.y += linearAcceleration.y * dt;

            constraint_count = 0;

            //brute force detect collisions (very inefficient- should change this!)
            for (int ii = 0; ii < NUM_RIGID_BODIES; ++ii) {
            	if(i != ii) { //can't collide with self
            		Vector2 contact_points[2];
            		Edge e = EPA(rigidBody, &rigidBodies[ii], contact_points);
            		printf("EPA result between %d and %d:\n", i, ii);
            		if(e.distance < 0.0) {
            			printf("\t no collision found during GJK\n");
            		}
            		else {
            			printf("\tdistance: %f\n\tnorm: %.2f, %.2f\n\tcontact_point A: %.2f, %.2f\n\tcontact_point B: %.2f, %.2f\n", e.distance, e.normal.x, e.normal.y,
									 contact_points[0].x, contact_points[0].y, contact_points[1].x, contact_points[1].y);
            			constraints[constraint_count].r_a.x = contact_points[0].x;
            			constraints[constraint_count].r_a.y = contact_points[0].y;
            			constraints[constraint_count].n.x = e.normal.x;
            			constraints[constraint_count].n.y = e.normal.y;
            			constraint_count++;
            		}
            	}
			}

			Vector2 m_inv_c1;
			Vector2 m_inv_c2;

			Vector4x2 J;
			Vector4 Jv2;
			Vector4x2 MinvJt;
			Vector4x4 JMinvJt;
			Vector4 sigma;
			Vector2 Pc;
			Vector2 impulse;

			double magnitude;

			//iterative constraint solving
			int ii = 0;
			while(ii < MAX_SI_ITERATIONS) {
				max_impulse = -1.0;
				constraint_iterator = 0;
				while(constraint_iterator < constraint_count) {

					//in this section I convert a big equation into c instructions
					//if you're wondering why I didn't just use a library or make my own methods for the math, it's because:
					// A: this is the easiest way to ensure correctness while handling transposes and such
					// B: it's pretty efficient

					//create J a 4x1 matrix [-n, -(p-x1) x n, n, (p-x2) x n]
					//where n is the collision normal and p is the collision point

					J.a = VectorNegate(constraints[constraint_iterator].n);
					J.b = Vec2Vec2Multiply(VectorNegate(constraints[constraint_count].r_a), constraints[constraint_iterator].n);
					J.c = VectorNegate(constraints[constraint_iterator].n);
					J.d = Vec2Vec2Multiply(constraints[constraint_count].r_a, constraints[constraint_iterator].n);

					//need Jv2, which is the 4x2 matrix J multiplied by the 2x1 matrix v2 resulting in a 4x1 matrix
					Jv2.a = (J.a.x * rigidBody->linearVelocity.x) + (J.a.y * rigidBody->linearVelocity.y);
					Jv2.b = (J.b.x * rigidBody->linearVelocity.x) + (J.b.y * rigidBody->linearVelocity.y);
					Jv2.c = (J.c.x * rigidBody->linearVelocity.x) + (J.c.y * rigidBody->linearVelocity.y);
					Jv2.d = (J.d.x * rigidBody->linearVelocity.x) + (J.d.y * rigidBody->linearVelocity.y);

					//add b
					Jv2.a = Jv2.a + SI_BIAS_FACTOR;
					Jv2.b = Jv2.b + SI_BIAS_FACTOR;
					Jv2.c = Jv2.c + SI_BIAS_FACTOR;
					Jv2.d = Jv2.d + SI_BIAS_FACTOR;

					//now we need the inverse mass
					m_inv_c1.x = 1.0 / rigidBody->shape.mass;
					m_inv_c1.y = 0;
					m_inv_c2.x = 0;
					m_inv_c2.y = 1.0 / rigidBody->shape.momentOfInertia;

					//get M_inverse * J_t giving us a 2x4 matrix (in memory as its transpose)
					MinvJt.a.x = (m_inv_c1.y * J.a.y) + (m_inv_c2.y * J.a.x);
					MinvJt.a.y = (m_inv_c1.x * J.a.y) + (m_inv_c2.x * J.a.x);
					MinvJt.b.x = (m_inv_c1.y * J.b.y) + (m_inv_c2.y * J.b.x);
					MinvJt.b.y = (m_inv_c1.x * J.b.y) + (m_inv_c2.x * J.b.x);
					MinvJt.c.x = (m_inv_c1.y * J.c.y) + (m_inv_c2.y * J.c.x);
					MinvJt.c.y = (m_inv_c1.x * J.c.y) + (m_inv_c2.x * J.c.x);
					MinvJt.d.x = (m_inv_c1.y * J.d.y) + (m_inv_c2.y * J.d.x);
					MinvJt.d.y = (m_inv_c1.x * J.d.y) + (m_inv_c2.x * J.d.x);

					//get J * MinvJt giving us a 4x4 matrix

					JMinvJt.a.a = (J.a.x * MinvJt.a.y) + (J.a.y * MinvJt.a.x);
					JMinvJt.a.b = (J.b.x * MinvJt.a.y) + (J.b.y * MinvJt.a.x);
					JMinvJt.a.c = (J.c.x * MinvJt.a.y) + (J.c.y * MinvJt.a.x);
					JMinvJt.a.d = (J.d.x * MinvJt.a.y) + (J.d.y * MinvJt.a.x);

					JMinvJt.b.a = (J.a.x * MinvJt.b.y) + (J.a.y * MinvJt.b.x);
					JMinvJt.b.b = (J.b.x * MinvJt.b.y) + (J.b.y * MinvJt.b.x);
					JMinvJt.b.c = (J.c.x * MinvJt.b.y) + (J.c.y * MinvJt.b.x);
					JMinvJt.b.d = (J.d.x * MinvJt.b.y) + (J.d.y * MinvJt.b.x);

					JMinvJt.c.a = (J.a.x * MinvJt.c.y) + (J.a.y * MinvJt.c.x);
					JMinvJt.c.b = (J.b.x * MinvJt.c.y) + (J.b.y * MinvJt.c.x);
					JMinvJt.c.c = (J.c.x * MinvJt.c.y) + (J.c.y * MinvJt.c.x);
					JMinvJt.c.d = (J.d.x * MinvJt.c.y) + (J.d.y * MinvJt.c.x);

					JMinvJt.d.a = (J.a.x * MinvJt.d.y) + (J.a.y * MinvJt.d.x);
					JMinvJt.d.b = (J.b.x * MinvJt.d.y) + (J.b.y * MinvJt.d.x);
					JMinvJt.d.c = (J.c.x * MinvJt.d.y) + (J.c.y * MinvJt.d.x);
					JMinvJt.d.d = (J.d.x * MinvJt.d.y) + (J.d.y * MinvJt.d.x);

					//divide by one (WARNING: MIGHT NEED TO BE TAKING INVERSE INSTEAD)

					JMinvJt.a.a = 1 / JMinvJt.a.a;
					JMinvJt.a.b = 1 / JMinvJt.a.b;
					JMinvJt.a.c = 1 / JMinvJt.a.c;
					JMinvJt.a.d = 1 / JMinvJt.a.d;

					JMinvJt.b.a = 1 / JMinvJt.b.a;
					JMinvJt.b.b = 1 / JMinvJt.b.b;
					JMinvJt.b.c = 1 / JMinvJt.b.c;
					JMinvJt.b.d = 1 / JMinvJt.b.d;

					JMinvJt.c.a = 1 / JMinvJt.c.a;
					JMinvJt.c.b = 1 / JMinvJt.c.b;
					JMinvJt.c.c = 1 / JMinvJt.c.c;
					JMinvJt.c.d = 1 / JMinvJt.c.d;

					JMinvJt.d.a = 1 / JMinvJt.d.a;
					JMinvJt.d.b = 1 / JMinvJt.d.b;
					JMinvJt.d.c = 1 / JMinvJt.d.c;
					JMinvJt.d.d = 1 / JMinvJt.d.d;

					//now get JMinvJt * Jv2, resulting in a 1x4 matrix sigma

					sigma.a = (JMinvJt.a.a * Jv2.a) + (JMinvJt.b.a * Jv2.b) + (JMinvJt.c.a * Jv2.c) + (JMinvJt.d.a * Jv2.d);
					sigma.b = (JMinvJt.a.b * Jv2.a) + (JMinvJt.b.b * Jv2.b) + (JMinvJt.c.b * Jv2.c) + (JMinvJt.d.b * Jv2.d);
					sigma.c = (JMinvJt.a.c * Jv2.a) + (JMinvJt.b.c * Jv2.b) + (JMinvJt.c.c * Jv2.c) + (JMinvJt.d.c * Jv2.d);
					sigma.d = (JMinvJt.a.d * Jv2.a) + (JMinvJt.b.d * Jv2.b) + (JMinvJt.c.d * Jv2.c) + (JMinvJt.d.d * Jv2.d);

					//negate it

					sigma.a = -sigma.a;
					sigma.b = -sigma.b;
					sigma.c = -sigma.c;
					sigma.d = -sigma.d;

					//now multiply Jt by sigma, resulting in a 1x2 vector Pc

					Pc.x = (J.a.y * sigma.a) + (J.b.y * sigma.b) + (J.c.y * sigma.c) + (J.d.y * sigma.d);
					Pc.y = (J.a.x * sigma.a) + (J.b.x * sigma.b) + (J.c.x * sigma.c) + (J.d.x * sigma.d);

					//multiply Minv by sigma to get the 1x2 vector impulse

					impulse.x = (m_inv_c1.x * Pc.x) + (m_inv_c2.x * Pc.y);
					impulse.y = (m_inv_c1.y * Pc.x) + (m_inv_c2.y * Pc.y); 

					//add the impulse to the velocity
					rigidBody->linearVelocity.x += impulse.x;
					rigidBody->linearVelocity.y += impulse.y;

					//get magnitude of impulse
					magnitude = sqrt((impulse.x * impulse.x) + (impulse.y * impulse.y));
					if(magnitude > max_impulse) {
						max_impulse = magnitude;
					}

					constraint_iterator++;
				}

				if(max_impulse <= SI_TOLERANCE) { 
					//corrective impulses small enough to quit
					ii = MAX_SI_ITERATIONS;
				}
				ii++;
			}



            rigidBody->position.x += rigidBody->linearVelocity.x * dt;
            rigidBody->position.y += rigidBody->linearVelocity.y * dt;
            float angularAcceleration = rigidBody->torque / rigidBody->shape.momentOfInertia;
            rigidBody->angularVelocity += angularAcceleration * dt;
            rigidBody->angle += rigidBody->angularVelocity * dt;
        }
        
        currentTime += dt;
        PrintRigidBodies(currentTime, fp);
    }

    fclose(fp);
}

int main( int argc, const char* argv[] )
{
	printf("\nWAR GAME SIMULATION\n\n");
	RunRigidBodySimulation();
}
