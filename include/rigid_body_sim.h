#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#define NUM_RIGID_BODIES 5
#define TOLERANCE 0.00001
#include "neuralnet.h"

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
void CalculateBoxInertia(BoxShape *boxShape);
float * Rigidbody_get_input(RigidBody * rb);
void InitializeRigidBodies();
void ComputeForceAndTorque(RigidBody *rigidBody) ;
void FreeList(struct Node* head) ;
int PrintList(struct Node* head) ;
int GetEmpty(int list[3]) ;

int GetLast(int list[3]) ;
int GetB(int list[3], int a) ;
int GetC(int list[3], int a) ;
int GetSize(int list[3]) ;
double DotProduct(Vector2 a, Vector2 b) ;
Vector2 VectorSub(Vector2 a, Vector2 b) ;
Vector2 VectorNegate(Vector2 a) ;

double VectorNorm(Vector2 a) ;
Vector2 VectorDivideConstant(Vector2 a, double b) ;
Vector2 VectorMultiplyConstant(Vector2 a, double b);
Vector2 VectorNormalize(Vector2 a);


Vector2 VectorTripleProduct(Vector2 a, Vector2 b, Vector2 c) ;
Vector2 AffineTransform(Vector2 vertex, Vector2 position, float angle);
Vector2 GetFarthestPointInDirection(Vector2 *vertices, int count, Vector2 d); 
Vector2 Support(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 d) ;

int GJK(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo, Vector2 *simplex);

Edge FindClosestEdge(struct Node *simplex_head);
Edge EPA(RigidBody *rigidBodyOne, RigidBody *rigidBodyTwo);

void PrintRigidBodies(float time, FILE *fp);
void RigidBody_get_input(RigidBody * rb,float * input);
void apply_nn_mvt(neuralnet *player,RigidBody *rb,float *x,float *y);
int RunRigidBodySimulation(neuralnet ** players,int print);





 
#endif
