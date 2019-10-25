all: main.c 
	gcc -g -Wall -o /home/jack/828/war_sim/main.out main.c -lm -lbsd

clean: 
	$(RM) /home/jack/828/war_sim/main.out
