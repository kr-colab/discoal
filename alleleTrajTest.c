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
	N = 10000;
	dt = 1.0 / (400.0 * N);
	alpha = 1000;
    int reps = 1;
    for(j=0; j<reps;j++){
        createTrajectory(N,alpha,dt,currentTraj);
	    i=0;
        while(currentTraj[i] > (1.0/(2*N))){
            printf("%g %g\n",currentTraj[i],i*dt);
            i++;
        }
	    sum += i * dt;
    }
    printf("N: %g\n",N);
    printf("alpha: %g\n",alpha);
    printf("mean time: %g\n",sum/((float) reps));
    printf("mean time x2: %g\n",2 * N * sum/((float) reps));
}

void createTrajectory(int N,double alpha,double dt, double *currentTrajectory){
	float freq = ((2*N) - 1.0) / (2*N) ; 
	int i;
    // zero out vector
    for(i=0;i<MAXTRAJ;i++){
        currentTrajectory[i] = 0.0;
    }
    i = 0;
	while(freq > (1.0/(2*N))){
		freq = 1.0 - genicSelectionStochasticForwards(dt, (1.0 - freq), alpha);
		currentTrajectory[i++]= freq;
	}
	
}
