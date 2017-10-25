#include "alleleTraj.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ranlib.h"
#include <unistd.h>
#include <fcntl.h>

unsigned long devrand(void);
#define MAXTRAJ 10000000

void createTrajectory(int N,double alpha,double dt, double *currentTrajectory);

int main (int argc, const char * argv[]) {
	double freq, t, dt, N, alpha, s;
	int i,j;
	double sum = 0;
	long seed1, seed2;
	double *currentTraj;
	
	seed1 = (long) (devrand() % 2147483399);
	seed2 = (long) (devrand() % 2147483399);
	setall(seed1, seed2 );
	currentTraj = malloc(sizeof(double) * MAXTRAJ);
	N = 20000;
	dt = 1.0 / (400.0 * N);
	s = 0.0002;
	alpha = 2*N*s;
	createTrajectory(N,alpha,dt,currentTraj);
	
	i=0;
	while(currentTraj[i] > (1.0/(2*N))){
		printf("%lf\n",currentTraj[i++]);
	}
	printf("%d\n",i);
	
}

void createTrajectory(int N,double alpha,double dt, double *currentTrajectory){
	float freq = ((2*N) - 1.0) / (2*N) ; 
	int i=0;
	while(freq > (1.0/(2*N))){
		freq = 1.0 - genicSelectionStochasticForwards(dt, (1.0 - freq), alpha);
		currentTrajectory[i++]= freq;
	}
	
}