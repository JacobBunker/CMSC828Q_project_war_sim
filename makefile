CC=gcc
CFLAGS=-g -Wall -O3 

all: main test sim

sim: sim.o rigid_body_sim.o neuralnet.o arrays.o GA.o tictactoe.o
	$(CC) $(CFLAGS) -o  sim  sim.o rigid_body_sim.o arrays.o GA.o neuralnet.o tictactoe.o -lm -lbsd


main: main.o neuralnet.o arrays.o tictactoe.o GA.o rigid_body_sim.o
	$(CC) $(CFLAGS) -o main main.o neuralnet.o arrays.o GA.o tictactoe.o rigid_body_sim.o -lm  -lbsd

test: test.o neuralnet.o arrays.o tictactoe.o GA.o rigid_body_sim.o
	$(CC) $(CFLAGS) -o test test.o neuralnet.o arrays.o GA.o tictactoe.o rigid_body_sim.o  -lm  -lbsd

sim.o:  sim.c rigid_body_sim.h rigid_body_sim.c neuralnet.c neuralnet.h arrays.h arrays.c GA.c GA.h
	$(CC) $(CFLAGS) -c sim.c  -lm -lbsd

rigid_body_sim.o: rigid_body_sim.h rigid_body_sim.c neuralnet.c neuralnet.h arrays.h arrays.c
	$(CC) $(CFLAGS) -c -lbsd rigid_body_sim.c  -lm -lbsd  

main.o: main.c neuralnet.h neuralnet.c arrays.h arrays.c
	$(CC) $(CFLAGS) -c main.c
test.o: test.c neuralnet.h neuralnet.c arrays.h arrays.c
	$(CC) $(CFLAGS) -c test.c

tictactoe.o:tictactoe.h tictactoe.c neuralnet.h neuralnet.c arrays.h arrays.c
	$(CC) $(CFLAGS) -c tictactoe.c

neuralnet.o: neuralnet.h neuralnet.c arrays.h arrays.c
	$(CC) $(CFLAGS) -c neuralnet.c

GA.o:  arrays.h arrays.c rigid_body_sim.h rigid_body_sim.c GA.h GA.c
	$(CC) $(CFLAGS) -c GA.c  

arrays.o: arrays.h arrays.c
	$(CC) $(CFLAGS) -c arrays.c 
clean:
	rm main *o *~
