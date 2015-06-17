// discoal.c 
//
// a discrete sequence coalescent simulation
// A. Kern 11/16/12
//




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include "ranlib.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "alleleTraj.h"



int locusNumber; 
int untilMode = 0;
char *fileName;
double uTime;



void getParameters(int argc,const char **argv);
void usage();

int main(int argc, const char * argv[]){
	int i,j,ii,activeSweepFlag, totalSimCount;
	double currentSize,nextTime, currentFreq;
	long seed1, seed2;
	
	
	getParameters(argc,argv);
	double N = EFFECTIVE_POPN_SIZE; // effective population size
	
	seed1 = (long) abs(devrand() % 2147483399);
	seed2 = (long) abs(devrand() % 2147483399);
	setall(seed1, seed2 );
	
	//Hudson style header
	for(i=0;i<argc;i++)printf("%s ",argv[i]);
	printf("\n%ld %ld\n", seed1, seed2);
	
	currentTime=0;
	nextTime=999;
	currentSize=1.0;
	i = 0;
        totalSimCount = 0;
	while(i < sampleNumber){
		initialize();
		j=0;
		activeSweepFlag = 0;
		for(j=0;j<eventNumber && alleleNumber > 1;j++){
			if(j == eventNumber - 1){
				nextTime = MAXTIME;
				}
			else{
				nextTime = events[j+1].time;
			}
		//	printf("type: %c size: %f currentTime: %f nextTime: %f alleleNumber: %d\n",events[j].type,events[j].popnSize,
		//		currentTime,nextTime,alleleNumber);	
			switch(events[j].type){
				case 'n':
				currentTime = events[j].time;
				currentSize = events[j].popnSize;
				if(activeSweepFlag == 0){
					currentTime = neutralPhase(&breakPoints[0], currentTime, nextTime, currentSize);
				}
				else{
					currentTime = sweepPhaseEventsWithRecurrentMut(&breakPoints[0], currentTime, nextTime, sweepSite, \
						 currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
				}
			//				printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],\
										sweepPopnSizes[0]);
				break;
				case 's':
				assert(activeSweepFlag == 0);
				activeSweepFlag = 1;
				currentTime = events[j].time;
				currentSize = events[j].popnSize;
				currentFreq = 1.0 - (1.0 / (2.0 * N * currentSize));
			//	printf("event%d currentTime: %f nextTime: %f popnSize: %f\n",j,currentTime,nextTime,currentSize);
				currentTime = sweepPhaseEventsWithRecurrentMut(&breakPoints[0], currentTime, nextTime, sweepSite, \
					 currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
			//	printf("currentFreqAfter: %f alleleNumber:%d\n",currentFreq,alleleNumber);
			//	printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],\
							sweepPopnSizes[0]);
						
				break;
				case 'p':
				currentTime = events[j].time;
				currentSize = events[j].popnSize;
				mergePopns(1, 0);
				if(activeSweepFlag == 0){
				//	printf("here\n");
			//		printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],\
						sweepPopnSizes[0]);
					currentTime = neutralPhase(&breakPoints[0], currentTime, nextTime, currentSize);
				//		printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],\
							sweepPopnSizes[0]);
				}
				else{
					currentTime = sweepPhaseEventsWithRecurrentMut(&breakPoints[0], currentTime, nextTime, sweepSite, \
						 currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
				}
				break;
			}
			
		}
		//finish up the coalescing action!
		if(alleleNumber > 1){
			currentTime = neutralPhase(&breakPoints[0], currentTime, MAXTIME, currentSize);
		}
		//assign root
	//	root = nodes[0];
		//add Mutations
		if(untilMode==0)
			dropMutations();
		else
			dropMutationsUntilTime(uTime);	

		if(condRecMode == 0){
			//Hudson style output
			makeGametesMS(argc,argv);
                        i++;
		}
		else{
			if(condRecMet == 1){
                                i++;
				makeGametesMS(argc,argv);
				condRecMet = 0;
			}
	
		}
		freeTree(nodes[0]);
                totalSimCount += 1;
	}
        if(condRecMode == 1)
        {
            fprintf(stderr, "Needed run %d simulations to get %d with a recombination event within the specified bounds.\n", totalSimCount, i);
        }
	return(0);
}



void getParameters(int argc,const char **argv){
	int args;
	int i;
	
	if( argc < 3){
		usage();
	}
	//locusNumber = atoi(argv[1]);
	sampleSize = atoi(argv[1]);
	sampleNumber = atoi(argv[2]);
	nSites = atoi(argv[3]);
	args = 4;

	npops = 1;
	popnSizes[0]=sampleSize;
	popnSizes[1]=0;
	sampleSizes[0]=sampleSize;
	sampleSizes[1]=0;
	rho = 0.0;
	my_gamma = 0.0;
	gcMean = 0;
	theta = 0.0;
	alpha = 0.0;
	lambda = 0.0;
	tau = 0;
	ancestralSizeRatio = 1.0;
	f0=0.0;
        uA=0.0;
	sweepSite = 0.5;
	tDiv=666;
	priorTheta=priorRho=priorAlpha=priorTau=priorX=priorF0=priorUA=0;

	eventFlag = 1;
	effectiveSampleSize = sampleSize;
	finiteOutputFlag = 0;
	outputStyle = 'h';
	mask = 0;
	
	
	eventNumber = 0;
	events[eventNumber].time = 0.0;
	events[eventNumber].popnSize = 1.0;
	events[eventNumber].type = 'n';
	eventNumber++;

	condRecMode= 0;
	while(args < argc){
		switch(argv[args][1]){
			case 'S' :
			runMode = 'S';
			fileName = argv[++args];
			break;
			case 's' :
			segSites =  atoi(argv[++args]);
			break;
			case 't' :
			theta = atof(argv[++args]);
			break;
			case 'r' :
			rho = atof (argv[++args]);
			break;
			case 'g' :
			my_gamma = atof (argv[++args]);
			gcMean = atoi (argv[++args]);
			break;
			case 'a' :
			alpha = atof(argv[++args]);
			break;
			case 'x' :
			sweepSite = atof(argv[++args]);
			break;
			case 'd' :
			tDiv =  atof(argv[++args]);
			events[eventNumber].time = tDiv;
			events[eventNumber].type = 'p'; //pop split
			events[eventNumber].popnSize = 1.0; //anc size
			eventNumber++;
			break;
			case 'p' :
			npops = atoi(argv[++args]);
			for(i=0;i<npops;i++)sampleSizes[i]=atoi(argv[++args]);
			break;
			case 'e' :
			events[eventNumber].time = atof(argv[++args]);
			events[eventNumber].popnSize = atof(argv[++args]);
			events[eventNumber].type = 'n'; //neutral change
			eventNumber++;
			break;
			case 'w':
			switch(argv[args][2]){
				case 'd':
				sweepMode = 'd';
				break;
				case 's':
				sweepMode = 's';
				break;
				case 'n':
				sweepMode = 'N';
				break;
			}
			tau = atof(argv[++args]);
			events[eventNumber].time = tau;
			events[eventNumber].popnSize = events[eventNumber  - 1].popnSize;
			events[eventNumber].type = 's'; //sweep event
			eventNumber++;
			break;
			case 'f':
			f0 = atof(argv[++args]);
			break;
                        case 'u':
                        uA = atof(argv[++args]);
                        break;
			case 'P' :
			  switch(argv[args][2]){
				case 't':
				  priorTheta = 1;
				  pThetaLow=atof(argv[++args]);
				  pThetaUp=atof(argv[++args]);
				break;
				case 'r':
                                  if (strlen(argv[args]) == 4 && argv[args][3] == 'e'){
					priorRho = 2;
					pRhoMean=atof(argv[++args]);
                                  }
                                  else{
					priorRho = 1;
					pRhoLow=atof(argv[++args]);
					pRhoUp=atof(argv[++args]);
                                  }
				break;
				case 'a':
				  priorAlpha = 1;
				  pAlphaLow=atof(argv[++args]);
				  pAlphaUp=atof(argv[++args]);
				break;
				case 'u':
                                	if (strlen(argv[args]) == 4 && argv[args][3] == 'A'){
						priorUA = 1;
						pUALow=atof(argv[++args]);
						pUAUp=atof(argv[++args]);
					}
                                	else{
						priorTau = 1;
						pTauLow=atof(argv[++args]);
						pTauUp=atof(argv[++args]);
					}
				break;
				case 'x':
				  priorX = 1;
				  pXLow=atof(argv[++args]);
				  pXUp=atof(argv[++args]);
				break;
				case 'f':
				  priorF0 = 1;
				  pF0Low=atof(argv[++args]);
				  pF0Up=atof(argv[++args]);
				break;
				case 'e':
				switch(argv[args][3]){
					case '1':
					priorE1 = 1;
					pE1TLow=atof(argv[++args]);
					pE1THigh=atof(argv[++args]);
					pE1SLow=atof(argv[++args]);
					pE1SHigh=atof(argv[++args]);
					events[eventNumber].type = 'n';
					eventNumber++;
					break;
					case '2':
					priorE2 = 1;
					pE2TLow=atof(argv[++args]);
					pE2THigh=atof(argv[++args]);
					pE2SLow=atof(argv[++args]);
					pE2SHigh=atof(argv[++args]);
					events[eventNumber].type = 'n';
					eventNumber++;
					break;
				}
				break;
			}
                        break;
			case 'U' :
			untilMode= 1;
			uTime=atof(argv[++args]);
			break; 
			case 'C' :
			condRecMode= 1;
			condRecMet = 0;
			lSpot=atoi(argv[++args]);
			rSpot=atoi(argv[++args]);
			break;	  
		}
		args++;
	}
	sortEventArray(events,eventNumber);

	for(i=0;i<eventNumber;i++){
		if(events[i].type == 's'){
			events[i].popnSize = events[i-1].popnSize;
		}
	//	printf("event %d time: %f popnSize: %f\n",i,events[i].time,events[i].popnSize);
	}
//	printf("%s %f %f\n",&sweepMode,tau,f0);
}


void usage(){
	fprintf(stderr,"usage: discoal sampleSize sampleNumber nSites -ws tau -e time1 size1 -e time2 size2\n");
	fprintf(stderr,"parameters: \n");

	fprintf(stderr,"\t -t theta \n");
	fprintf(stderr,"\t -r rho (=zero if not specified)\n");
	fprintf(stderr,"\t -g conversionRate tractLengthMean (gene conversion)\n");
	fprintf(stderr,"\t -p npops sampleSize1 sampleSize2\n");
	fprintf(stderr,"\t -d divTime\n");	
	fprintf(stderr,"\t -ws tau (sweep happend tau generations ago- stochastic sweep)\n");  
	fprintf(stderr,"\t -wd tau (sweep happend tau generations ago- deterministic sweep)\n"); 
	fprintf(stderr,"\t -wn tau (sweep happend tau generations ago- neutral sweep)\n");
	fprintf(stderr,"\t -f first frequency at which selection acts on allele (F0; sweep models only)\n");
	fprintf(stderr,"\t -uA rate at which adaptive mutation recurs during the sweep phase (sweep models only)\n");
	fprintf(stderr,"\t -a alpha (=2Ns)\n");
	fprintf(stderr,"\t -x sweepSite (0-1)\n");	
	fprintf(stderr,"\t -Pt low high (prior on theta)\n");
	fprintf(stderr,"\t -Pr low high (prior on rho)\n");
	fprintf(stderr,"\t -Pre mean (prior on rho -- exponentially distributed)\n");
	fprintf(stderr,"\t -Pa low high (prior on alpha)\n");
	fprintf(stderr,"\t -Pu low high (prior on tau; sweep models only)\n");
	fprintf(stderr,"\t -PuA low high (prior on uA; sweep models only)\n");
	fprintf(stderr,"\t -Px low high (prior on sweepSite; sweep models only)\n");
	fprintf(stderr,"\t -Pf low high (prior on F0; sweep models only)\n");
	fprintf(stderr,"\t -Pe1 lowTime highTime lowSize highSize (priors on first demographic move time and size)\n");
	fprintf(stderr,"\t -Pe2 lowTime highTime lowSize highSize (priors on second demographic move time and size)\n");
	fprintf(stderr,"\t -U time (only record mutations back to specified time)\n");
	
	
	fprintf(stderr,"\t Note: all time units are 2N_0 generations!!!\n");
	exit(1);
}
