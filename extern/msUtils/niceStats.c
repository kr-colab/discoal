/******* niceStats.c ********
for calculating sample stats from MS output 
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000



void usage();



int maxsites = 100000 ;

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ,nwins, siteIdx;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , nadv, iss, h, exponent1, exponent2;
	double pi , th,  z, H, mfda, tajD,w, wins[50],max, min, temp_site, thetaA, thetaHPi, tajDX, achazD, h1, h2, h12, thetaW;
        double winsH1[50], winsH2[50], winsH12[50];
	char dum[50], astr[100] ;
        double *harmonicSums;
        double ehh, rehh;
        int *haplotype_counts;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);

	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	//print header
	nwins=5;
	//printf("pi\tss\tthetaH\ttajD\tfayWuH\tHapCount");
        //output used for soft shoulder analysis:
	//printf("pi\tss\tthetaH\ttajD\tfayWuH\tHapCount\tH1\tomegaCenter\tZnS");
        //output used for spatial svm:
	printf("pi\tss\tthetaH\ttajD\tfayWuH\tmaxFDA\tHapCount\tH1\tH12\tH2/H1\tOmega\tZnS");
        /*for (exponent1 = -20; exponent1 <= 20; exponent1++){
            for (exponent2 = -20; exponent2 < exponent1; exponent2++) {
                printf("\tachazsD_%d_%d", exponent1,exponent2);
            }
        }*/
        /*for (exponent1 = -50; exponent1 <= 50; exponent1++){
            printf("\tthetaA_%d", exponent1);
        }*/
	//for(i=0;i<nwins;i++)printf("\tH12_win%d\tH2/H1_win%d",i,i);
	printf("\n");
	while( howmany-count++ ) {

/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
		}while ( line[0] != '/' );

		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);

			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}
		/* analyse sample ( do stuff with segsites and list) */

	        haplotype_counts = (int *)malloc( nsam*sizeof( int ) ) ;

		iss = segSites(segsites,nsam,list);
		pi = nucdiv(nsam, segsites, list) ;
		mfda = maxFDA(nsam, segsites, list);
		th = thetah(nsam, segsites, list) ;
		h = nHaplotypes(segsites,nsam,list);
                getHaplotypeFreqSpec(segsites, nsam, list, haplotype_counts);
                h1 = petrovH1(haplotype_counts,nsam);
                h2 = petrovH2(haplotype_counts,nsam);
                h12 = petrovH12(haplotype_counts,nsam);
		H = th-pi;
		tajD = tajd(nsam,segsites,pi);
                //ehh = meanEHH(segsites, posit, 0.25, nsam, list);
                //rehh = meanREHH(segsites, posit, 0.25, nsam, list);
                //thetaHPi = achazThetaHPi(nsam, segsites, list);
                //thetaW = achazThetaExponentWeights(nsam,segsites,list,-1);
                //tajDX = achazTajimasDExtreme(nsam, segsites, list);
		z = ZnS( segsites,  nsam,  list);
		/*OmegaCenter*/
		//Get snp that is nearest to our fixation but on the right (omegaCenter treats this focal SNP as being in the right window)
		/*
                min = (fabs(0.5 - posit[0]));
                siteIdx = segsites-1;
		for (i = 0; i < segsites; ++i){
			temp_site = posit[i];
				if (temp_site >= 0.5 && fabs(0.5 - temp_site) <= min){
					siteIdx = i;
					min = fabs(0.5-temp_site);
				}
		}
		w = omegaCenter(siteIdx, segsites, nsam, list);
		*/
                w = omegaMax(segsites, nsam,list);
		
		//printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%f\t%f", pi, iss, th ,tajD,H, h, z, w);
		//output used for soft shoulders manuscript:
                //printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%f\t%f\t%f", pi, iss, th ,tajD,H, h, h1, w, z);
                //output used for spatial SVM
		printf("%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf", pi, iss, th, tajD, H, mfda, h, h1, h12, h2/h1, w, z);
		//printf("%lf\t%d\t%lf\t%lf\t%lf\t%d", pi, iss, th ,tajD,H, h);
                
                /*for (exponent1 = -50; exponent1 <= 50; exponent1++){
                    thetaA = achazThetaExponentWeights(nsam,segsites,list,exponent1);
                    printf("\t%f", thetaA);
                }*/

                /*harmonicSums = compute_HarmonicSums(nsam);
                for (exponent1 = -20; exponent1 <= 20; exponent1++){
                    for (exponent2 = -20; exponent2 < exponent1; exponent2++) {
                        achazD = achazNeutTestExponentWeights(nsam,segsites,list,exponent1,exponent2,harmonicSums);
                        printf("\t%f", achazD);
                    }
                }*/


		//window stats
		//tajdWindow( nwins, posit, wins, nsam, segsites,list);
                //petrovHStatsWindow(segsites, nwins, posit, winsH12, winsH1, winsH2, nsam, list);
                /*
                //get max for normalizing windows
		max=0.0;
		for( i=0; i<nwins ; i++){
			if (wins[i] > max){
				max = wins[i];
			}
		} 
                //print windows
                
		for( i=0; i<nwins ; i++){
			if (max == 0.0){ // Prevent "-nans" in output. 
 				printf("\t%f", 0.0);
			}
			else{
				printf("\t%f\t%f",winsH12[i],winsH2[i]/winsH1[i]);
				//printf("\t%f",wins[i]/max);
			} 
	        }*/
        //free(harmonicSums);
        free(haplotype_counts);
	printf("\n");
 }
	return(0);
}
