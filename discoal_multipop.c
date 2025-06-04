// discoal_multipop.c 
//
// a discrete sequence coalescent simulation
// A. Kern 10/3/14
//




#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include "ranlib.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "alleleTraj.h"



int locusNumber; 
int leftRhoFlag=0;
int untilMode = 0;
const char *fileName;
double uTime;
double *currentSize;
long seed1, seed2;
//float *currentTrajectory;

void getParameters(int argc,const char **argv);
void usage();

int main(int argc, const char * argv[]){
	int i,j,k, totalSimCount;
	float tempSite;
	int lastBreak;
	double nextTime, currentFreq, probAccept;
	
	
	
	getParameters(argc,argv);
	double N = EFFECTIVE_POPN_SIZE; // effective population size
	setall(seed1, seed2 );

	//Hudson style header
	for(i=0;i<argc;i++)printf("%s ",argv[i]);
	printf("\n%ld %ld\n", seed1, seed2);
	
	i = 0;
        totalSimCount = 0;
	trajectoryCapacity = TRAJSTEPSTART;
	currentTrajectory = malloc(sizeof(float) * trajectoryCapacity);
	assert(currentTrajectory);
	
	// Initialize global trajectory generator (lazy approach)
	activeTrajectoryGen = NULL;

	while(i < sampleNumber){
		currentTime=0;
		nextTime=999;
		currentSize[0]=1.0;
		currentFreq = 1.0 - (1.0 / (2.0 * N * currentSize[0])); //just to initialize the value
//		printf("popnsize[0]:%d",popnSizes[0]);
		maxTrajSteps = trajectoryCapacity;
		
		
		
		initialize();

		j=0;
		activeSweepFlag = 0;
		for(j=0;j<eventNumber && alleleNumber > 1;j++){
			currentEventNumber=j; //need this annoying global for trajectory generation
			if(j == eventNumber - 1){
				nextTime = MAXTIME;
				}
			else{
				nextTime = events[j+1].time;
			}
			// printf("type: %c popID: %d size: %f currentTime: %f nextTime: %f alleleNumber: %d\n",events[j].type,events[j].popID,events[j].popnSize,
			//	currentTime,nextTime,alleleNumber);	
			switch(events[j].type){
				case 'n':
				currentTime = events[j].time;
				currentSize[events[j].popID] = events[j].popnSize;
				//for(i=0;i<npops;i++)
				//	for(j=0;j<npops;j++) printf("%f\n",migMat[i][j]);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
			//	printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],
			//							sweepPopnSizes[0]);
				break;
				case 's':
				assert(activeSweepFlag == 0);
				currentTime = events[j].time;
				if (partialSweepMode == 1){
					currentFreq = MIN(partialSweepFinalFreq,1.0 - (1.0 / (2.0 * N * currentSize[0])));
				}
				else{
					currentFreq = 1.0 - (1.0 / (2.0 * N * currentSize[0]));
				}
			//	printf("event%d currentTime: %f nextTime: %f popnSize: %f\n",j,currentTime,nextTime,currentSize);

				//generate a proposed trajectory using lazy approach
				TrajectoryGenerator *candidateGen = initializeTrajectoryGenerator(currentEventNumber, currentSize, sweepMode, currentFreq, alpha, f0, currentTime);
				probAccept = calculateTrajectoryAcceptance(candidateGen);
				while(ranf()>probAccept){
					cleanupTrajectoryGenerator(candidateGen);
					candidateGen = initializeTrajectoryGenerator(currentEventNumber, currentSize, sweepMode, currentFreq, alpha, f0, currentTime);
					probAccept = calculateTrajectoryAcceptance(candidateGen);
					//printf("probAccept: %lf\n",probAccept);
				}
				cleanupTrajectoryGenerator(candidateGen);
				
				currentTime = sweepPhaseEventsLazyTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
				//printf("currentFreqAfter: %f alleleNumber:%d currentTime:%f\n",currentFreq,alleleNumber,currentTime);
				//printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],
				//			sweepPopnSizes[0]);
				if (currentTime < nextTime)
                                        currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
						
				break;
				case 'p': //merging populations
				currentTime = events[j].time;
				//printf("here at P flag time=%f\n",currentTime);
				mergePopns(events[j].popID, events[j].popID2);

				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                        		currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				break;
				case 'a':
				currentTime = events[j].time;
				admixPopns(events[j].popID, events[j].popID2, events[j].popID3, events[j].admixProp);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				break;
				case 'A':
				currentTime = events[j].time;
				//printAllActiveNodes();
				addAncientSample(events[j].lineageNumber, events[j].popID, events[j].time, activeSweepFlag, currentFreq);
				//printAllActiveNodes();
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				break;
			}
			
		}
		//finish up the coalescing action!
		if(alleleNumber > 1){
			currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, MAXTIME, currentSize);
		}
		//assign root
	//	root = nodes[0];
		//add Mutations
		if(untilMode==0)
			dropMutations();
		else
			dropMutationsUntilTime(uTime);	

		if(condRecMode == 0){
			if(treeOutputMode == 1){
				//output newick trees
                qsort(breakPoints, breakNumber, sizeof(breakPoints[0]), compare_floats);
				lastBreak = 0;
				printf("\n//\n");
				for(k=0;k<breakNumber;k++){
					tempSite = ((float) breakPoints[k] / nSites) - (0.5/nSites) ; //padding
					if(breakPoints[k] - lastBreak > 0){
						printf("[%d]",breakPoints[k] - lastBreak);
                        //printf("%g\n",allNodes[findRootAtSite(breakPoints[k])]->time);
                        printTreeAtSite(tempSite); 
						lastBreak = breakPoints[k];
					}
				}
				printf("[%d]",nSites- lastBreak);
				//printf("%g\n",allNodes[findRootAtSite(1.0-(1.0/nSites))]->time);
				printTreeAtSite(1.0 - (1.0/nSites)); 

			}
			else{
				//Hudson style output
				//errorCheckMutations();
				makeGametesMS(argc,argv);
			}
			//printf("rep: %d\n",i);
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
		cleanupBreakPoints();
                totalSimCount += 1;
	}
        if(condRecMode == 1)
        {
            fprintf(stderr, "Needed run %d simulations to get %d with a recombination event within the specified bounds.\n", totalSimCount, i);
        }
	free(currentTrajectory);
	free(currentSize);
	
	// Clean up global trajectory generator if allocated
	if(activeTrajectoryGen != NULL) {
		cleanupTrajectoryGenerator(activeTrajectoryGen);
		activeTrajectoryGen = NULL;
	}
	
	return(0);
}



void getParameters(int argc,const char **argv){
	int args;
	int i,j;
	double migR;
	int selCheck,nChangeCheck;
	
	if( argc < 3){
		usage();
	}

	sampleSize = atoi(argv[1]);
	if(sampleSize > 254){
		#ifndef BIG
		printf("Error: sampleSize > 254. recompile discoal and set the -DBIG flag\n");
		exit(666);
		#endif
	}
	sampleNumber = atoi(argv[2]);
	nSites = atoi(argv[3]);
	if(nSites>MAXSITES){
		printf("Error: number of sites set higher than current compilation limit. Please reduce the number of sites or change the MAXSITES define and recompile\n");
		exit(666);
	}
	args = 4;

	npops = 1;
	popnSizes[0]=sampleSize;
	popnSizes[1]=0;
	sampleSizes[0]=sampleSize;
	sampleSizes[1]=0;
	leftRho = 0.0;
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

	seed1 = (long) (devrand() % 2147483399);
	seed2 = (long) (devrand() % 2147483399);
	

	EFFECTIVE_POPN_SIZE = 1000000;
	sweepSite = 0.5;
	tDiv=666;
        gammaCoRatioMode = 0;
	priorTheta=priorRho=priorAlpha=priorTau=priorX=priorF0=priorUA=priorC=0;

	eventFlag = 1;
	effectiveSampleSize = sampleSize;
	finiteOutputFlag = 0;
	outputStyle = 'h';
	mask = 0;
	migFlag = 0;
	deltaTMod = 40;
	recurSweepMode = 0;
	treeOutputMode= 0;
	partialSweepMode = 0;
	softSweepMode = 0;
	ancSampleFlag = 0;
	ancSampleSize = 0;
	hidePartialSNP = 0;
	
	//set up first bogus event
	eventNumber = 0;
	events[eventNumber].time = 0.0;
	events[eventNumber].popID = 0;
	events[eventNumber].popnSize = 1.0;
	events[eventNumber].type = 'n';
	eventNumber++;
	currentSize = malloc(sizeof(double) * MAXPOPS);

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
			case 'i' :
			deltaTMod = atof(argv[++args]);
			break;
			case 'r' :
			rho = atof (argv[++args]);
			break;
			case 'g' :
                        if (argv[args][2] == 'r')
                        {
			  gammaCoRatio = atof (argv[++args]);
			  gcMean = atoi (argv[++args]);
                          gammaCoRatioMode = 1;
                        }
                        else
                        {
			  my_gamma = atof (argv[++args]);
			  gcMean = atoi (argv[++args]);
                        }
			break;
			case 'a' :
			alpha = atof(argv[++args]);
			break;
			case 'x' :
			sweepSite = atof(argv[++args]);
			break;
			case 'M' :
			if(npops==1){
				fprintf(stderr,"Error: attempting to set migration but only one population! Be sure that 'm' flags are specified after 'p' flag\n");
				exit(1);
			}
			migR = atof(argv[++args]);
			for(i=0;i<npops;i++){
				for(j=0;j<npops;j++){
					if(i!=j){
						migMatConst[i][j]=migR;
					}
					else{
						migMatConst[i][j]= 0.0;
					}
				}
			}
			migFlag = 1;
			break;
			case 'm' :
			if(npops==1){
				fprintf(stderr,"Error: attempting to set migration but only one population! Be sure that 'm' flags are specified after 'p' flag\n");
				exit(1);
			}
			i = atoi(argv[++args]);
			j = atoi(argv[++args]);
			migR = atof(argv[++args]);
			migMatConst[i][j]=migR;
			migFlag = 1;
			break;
			case 'p' :
			npops = atoi(argv[++args]);
			if(npops > MAXPOPS){
				fprintf(stderr,"Error: too many populations defined. Current maximum number = %d. Change MAXPOPS define in discoal.h and recompile... if you dare\n",MAXPOPS);
				exit(1);
			}
			for(i=0;i<npops;i++){
				sampleSizes[i]=atoi(argv[++args]);
				currentSize[i] = 1.0;
			}
			
			break;
			case 'e' :
				switch(argv[args][2]){
					case 'n':
					events[eventNumber].time = atof(argv[++args]) * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popnSize = atof(argv[++args]);
					events[eventNumber].type = 'n'; //size change
					eventNumber++;
					break;
					case 'd' :
					tDiv =  atof(argv[++args]);
					events[eventNumber].time = tDiv * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].type = 'p'; //pop split
					eventNumber++;
					break;
					case 'a' :
					events[eventNumber].time = atof(argv[++args]) * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].popID3 = atoi(argv[++args]);
					events[eventNumber].admixProp = atof(argv[++args]);
					events[eventNumber].type = 'a'; //admix split
					eventNumber++;
					break;
				}
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
				tau = atof(argv[++args]) * 2.0;
				events[eventNumber].time = tau;
				events[eventNumber].type = 's'; //sweep event
				eventNumber++;
				break;
			case 'l':
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
				sweepSite = -1.0;
				tau = atof(argv[++args]) * 2.0;
				leftRho = atof(argv[++args]) * 2.0;
				leftRhoFlag=1;
				events[eventNumber].time = tau;
				events[eventNumber].type = 's'; //sweep event
				eventNumber++;
				break;
			case 'f':
			f0 = atof(argv[++args]);
			softSweepMode = 1;
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
				case 'c':
				  priorC = 1;
			          partialSweepMode = 1;
				  pCLow=atof(argv[++args]);
				  pCUp=atof(argv[++args]);
				break;
                                case 'r':
                                  if (strlen(argv[args]) == 4 && argv[args][3] == 'e'){
                                        priorRho = 2;
                                        pRhoMean=atof(argv[++args]);
                                        pRhoUp=atof(argv[++args]);
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
					  pTauLow=atof(argv[++args]) * 2.0;
					  pTauUp=atof(argv[++args]) * 2.0;
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
					pE1TLow=atof(argv[++args])*2.0;
					pE1THigh=atof(argv[++args])*2.0;
					pE1SLow=atof(argv[++args]);
					pE1SHigh=atof(argv[++args]);
					events[eventNumber].type = 'n';
					eventNumber++;
					break;
					case '2':
					priorE2 = 1;
					pE2TLow=atof(argv[++args])*2.0;
					pE2THigh=atof(argv[++args])*2.0;
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
			uTime=atof(argv[++args])*2.0;
			break; 
			case 'd' :
			seed1=atoi(argv[++args]);
			seed2=atoi(argv[++args]);
			break;
			case 'N' :
			EFFECTIVE_POPN_SIZE=atoi(argv[++args]);
			break;
			case 'T' :
			treeOutputMode= 1;
			break;
			case 'C' :
			condRecMode= 1;
			condRecMet = 0;
			lSpot=atoi(argv[++args]);
			rSpot=atoi(argv[++args]);
			break;
			case 'R' :
			recurSweepMode = 1;
			sweepMode = 's';
			recurSweepRate = atof(argv[++args]);
			break;
			case 'L' :
			recurSweepMode = 1;
			sweepMode = 's';
			sweepSite = -1.0;
			recurSweepRate = atof(argv[++args]);
            if (recurSweepRate <= 0)
            {
                fprintf(stderr,"recurSweepRate must be > 0\n");
                exit(0);
            }
			break;
			case 'c' :
			partialSweepMode = 1;
			sweepMode = 's';
			partialSweepFinalFreq = atof(argv[++args]);
            if (partialSweepFinalFreq <= 0.0 || partialSweepFinalFreq >= 1.0)
            {
                fprintf(stderr,"partialSweepFinalFreq must be > 0 and < 1.0\n");
                exit(1);
            }
			break;
			case 'h' :
			hidePartialSNP = 1;
			break;
			case 'A' :
			events[eventNumber].lineageNumber = atoi(argv[++args]);
			events[eventNumber].popID = atoi(argv[++args]);	 
			events[eventNumber].time = atof(argv[++args]) * 2.0; 
			ancSampleSize += events[eventNumber].lineageNumber;
			events[eventNumber].type = 'A'; //ancient sample
			ancSampleFlag = 1;
			eventNumber++;
			assert(events[eventNumber].lineageNumber < sampleSize);
			break;
			 
		}
		args++;
	}
	sortEventArray(events,eventNumber);

	//make sure events are kosher
	selCheck = 0;
	nChangeCheck=0;
	for(i=0;i<eventNumber;i++){
		//printf("event %d: type is %c\n", i, events[i].type);
	 		if(events[i].type == 's'){
				selCheck = 1;
	 		}
			if(events[i].type == 'n'){
				nChangeCheck += 1;
	 		}
	}
	if(selCheck == 1){
		if(recurSweepMode == 1){
			printf("Error with event specification: a single sweep event has been found but recurrentSweep mode has been specified\n");
			exit(666);
		}
		if(nChangeCheck > 1 && sweepMode=='d'){
			printf("Error with event specification: you chose 1 or more population size changes with a deterministic sweep. Please us -ws flag instead\n");
			exit(666);
		}
		if(softSweepMode == 1 && partialSweepMode == 1){
			if(f0 >= partialSweepFinalFreq){
				printf("Error with event specification: you specified a partial soft sweep but final frequency of partial sweep <= f_0\n");
				exit(666);
			}
		}
	}
	if(leftRhoFlag && sweepSite >= 0.0){
		printf("Error with event specification: you chose leftRho mode but the sweep site is within the locus\n");
		exit(666);
	}
	if(softSweepMode == 1 && recurSweepMode == 1){
		printf("Error with event specification: currently recurrent soft sweeps are not implemented. this will be a future addition\n");
		exit(666);
	}
	
}
		


void usage(){
	fprintf(stderr,"usage: discoal sampleSize numReplicates nSites -ws tau\n");
	fprintf(stderr,"parameters: \n");

	fprintf(stderr,"\t -t theta \n");
	fprintf(stderr,"\t -r rho (=zero if not specified)\n");
	fprintf(stderr,"\t -g conversionRate tractLengthMean (gene conversion)\n");
	fprintf(stderr,"\t -gr conversionToCrossoverRatio tractLengthMean (gene conversion where initiation rate = rho*conversionToCrossoverRatio)\n");
	fprintf(stderr,"\t -p npops sampleSize1 sampleSize2 etc.\n");
	fprintf(stderr,"\t -en time popnID size (changes size of popID)\n");	
	fprintf(stderr,"\t -ed time popnID1 popnID2 (joins popnID1 into popnID2)\n");
	fprintf(stderr,"\t -ea time daughterPopnID founderPopnID1 founderPopnID2 admixProp (admixture-- back in time daughterPopnID into two founders)\n");
	
	fprintf(stderr,"\t -ws tau (sweep happend tau generations ago- stochastic sweep)\n");  
	fprintf(stderr,"\t -wd tau (sweep happend tau generations ago- deterministic sweep)\n"); 
	fprintf(stderr,"\t -wn tau (sweep happend tau generations ago- neutral sweep)\n");
	fprintf(stderr,"\t -ls tau leftRho (stochastic sweep some genetic distance to the left of the simulated window--specified by leftRho=4Nr)\n");
	fprintf(stderr,"\t\t similarly, ld and ln simulate deterministic and neutral sweeps to the left of the window, respectively\n");
	fprintf(stderr,"\t -f first frequency at which selection acts on allele (F0; sweep models only)\n");
	fprintf(stderr,"\t -uA rate at which adaptive mutation recurs during the sweep phase (sweep models only)\n");
	fprintf(stderr,"\t -N sweepEffectivePopnSize (sweep models only)\n");	
	fprintf(stderr,"\t -a alpha (=2Ns)\n");
	fprintf(stderr,"\t -x sweepSite (0-1)\n");
	fprintf(stderr,"\t -c partialSweepFinalFrequency (partial sweeps)\n");	
	fprintf(stderr,"\t -i dt (sweep time increment scalar; default 400 -> 1/400N)\n");
	
	fprintf(stderr,"\t -M migRate (sets all rates to migRate)\n");
	fprintf(stderr,"\t -m popnID1 popnID2 migRate (sets migRate from popnID1 to popnID2)\n");
	fprintf(stderr,"\t -A sampleSize popnID time (ancient sample from popnID at specified time)\n");
	
	fprintf(stderr,"\t -Pt low high (prior on theta)\n");
	fprintf(stderr,"\t -Pr low high (prior on rho)\n");
        fprintf(stderr,"\t -Pre mean upperBound (prior on rho -- exponentially distributed but truncated at an upper bound)\n");
	fprintf(stderr,"\t -Pa low high (prior on alpha)\n");
	fprintf(stderr,"\t -Pu low high (prior on tau; sweep models only; still must use \"-ws tau\" and \"tau\" will be ignored)\n");
	fprintf(stderr,"\t -PuA low high (prior on uA; sweep models only)\n");
	fprintf(stderr,"\t -Px low high (prior on sweepSite; sweep models only)\n");
	fprintf(stderr,"\t -Pf low high (prior on F0; sweep models only)\n");
	fprintf(stderr,"\t -Pc low high (prior on partialSweepFinalFreq; sweep models only)\n");
	fprintf(stderr,"\t -Pe1 lowTime highTime lowSize highSize (priors on first demographic move time and size)\n");
	fprintf(stderr,"\t -Pe2 lowTime highTime lowSize highSize (priors on second demographic move time and size)\n");
	//fprintf(stderr,"\t -U time (only record mutations back to specified time)\n");
	fprintf(stderr,"\t -R rhhRate (recurrent hitch hiking mode at the locus; rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -L rhhRate (recurrent hitch hiking mode to the side of locus; leftRho is ~Unif(0,4Ns); rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -h (hide selected SNP in partial sweep mode)\n");
	fprintf(stderr,"\t -T (tree output mode)\n");
	fprintf(stderr,"\t -d seed1 seed2 (set random number generator seeds)\n");
	
	exit(1);
}
