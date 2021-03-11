/*  Index of phases:
 *    P   - 0 -> ray1.info
 *    PP  - 1 -> ray2.info
 *    PcP - 2 -> ray3.info
 *    S   - 3 -> ray4.info
 *    SS  - 4 -> ray5.info
 *    ScS - 5 -> ray6.info
 *  you can add phases yourself, just give their indicators according to comments in rayinformation.c to make sure you get correct information
 *
 *  IN:
 *    sdepth - depth of source
 *    rdepth - depth of receivers
 *    Vr[0]  - Vs of near-receiver area
 *    Vr[1]  - Vp of near-receiver area
 *    Vm[3]  - Vs of PP/SS rebound point structure
 *    Vm[4]  - Vp of PP/SS rebound point structure
 *    Vm[5]  - rho of PP/SS rebound point structure
 *    modelname - name of input standary model
 *  OUT:
 *    first line is number of rays followed by lines consist of
 *    gcarc(deg), spherical parameter(/sec), geometrical spreading factor(1./R in /km), traveltime(sec), t star(t* in sec), reflection coefficient,take-off angle 
 *    
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rayinformation.c"

#define NUM 6

int main(){
	double sdepth;		// depth of source 
	double ddepth[NUM];	// depth of discontinuities for each ray
	double rdepth;		// depth of receivers
	int indr[NUM][4];	// indicators of ray (see above and rayinformation.c)
	char modelname[20];	// name of input standary 1-D model

	double Vr[2];		// structure of near-receiver layer
	double Vm[6];		// structure of PP/SS rebound point
	double temp[6];		// structure of CMB (you can modify code below to change it to another discontinuity)
	int i,layers;
	double depth,Vp,Vs;
	double rho,Qp,Qs;
	double depth_cmb;
	double dh[6];           //differential of hh when generating ray information 
	char outfile[20];
	int error;

	FILE *fp;

	scanf ("%lf",&sdepth);			// depth of source
	scanf ("%lf",&rdepth);			// depth of receivers
	scanf ("%lf%lf",&Vr[0],&Vr[1]);			// Vs and Vp near the receiver
	scanf ("%lf%lf%lf",&Vm[3],&Vm[4],&Vm[5]);	// Vs, Vp and Rho of PP/SS rebound point
	Vm[0] = 0.0;	//Vs of the layer above surface
	Vm[1] = 0.0;	//Vp of the layer above surface
	Vm[2] = 0.0;	//rho of the layer above surface
	scanf ("%s",modelname);
	dh[0] = 0.0005;
	dh[1] = 0.0005;
	dh[2] = 0.0001;
        dh[3] = 0.0005;
	dh[4] = 0.0005;
        dh[5] = 0.0001;

	if ( strncmp("none",modelname,4) == 0 ){
		sprintf (modelname,"%s","prem");
	}

	fp = fopen (modelname,"r");
	if ( fp == NULL ){
		fprintf (stderr,"Error in opening %s file\n",modelname);
		exit (-1);
	}
	if ( fscanf (fp,"%d",&layers) != 1 ){
		fprintf (stderr,"Error in reading number of layers in stdmodel file\n");
		exit (-1);
	}

	for ( i=0;i<layers;i++ ){
		if ( fscanf (fp,"%lf%lf%lf%lf%lf%lf",&depth,&Vp,&Vs,&rho,&Qp,&Qs) != 6 ){
			fprintf (stderr,"Error in reading information from stdmodel\n");
			exit (-1);
		}
		if (Vs < 0.1){		//CMB 
			depth_cmb = depth;
			temp[3] = Vs;
			temp[4] = Vp;
			temp[5] = rho;
			break;
		}
		temp[0] = Vs;	//structure information of CMB
		temp[1] = Vp;
		temp[2] = rho;
	}

	/* P indicator */
	indr[0][0] = 1;
	indr[0][1] = 1;
	indr[0][2] = 1;
	indr[0][3] = 0;
	ddepth[0]  = rdepth;


	/* PP indicator */
	indr[1][0] = 1;
	indr[1][1] = 1;
	indr[1][2] = 1;
	indr[1][3] = 1;
	ddepth[1]  = 0.;

	/* PcP indicator */
	indr[2][0] = 1;
	indr[2][1] = 1;
	indr[2][2] = 0;
	indr[2][3] = 0;
	ddepth[2]  = depth_cmb;
	
	/* S indicator */
	indr[3][0] = 0;
	indr[3][1] = 0;
	indr[3][2] = 1;
	indr[3][3] = 0;
	ddepth[3]  = rdepth;

	/* SS indicator */
	indr[4][0] = 0;
	indr[4][1] = 0;
	indr[4][2] = 1;
	indr[4][3] = 1;
	ddepth[4]  = 0.;

	/* ScS indicator */
	indr[5][0] = 0;
	indr[5][1] = 0;
	indr[5][2] = 0;
	indr[5][3] = 0;
	ddepth[5]  = depth_cmb;

	for ( i=0;i<NUM;i++ ){
		sprintf (outfile,"ray%1d.info",i+1);
		if (i != 2 && i != 5){
			error = rayinformation(modelname,sdepth,ddepth[i],rdepth,indr[i],Vm,Vr,outfile,dh[i]);
			if (error == -1){
				fprintf (stderr,"Error in running rayinformation\n");
				continue;
			}
		}else{
			error = rayinformation(modelname,sdepth,ddepth[i],rdepth,indr[i],temp,Vr,outfile,dh[i]);
			if (error == -1){
				fprintf (stderr,"Error in running rayinformation\n");
				continue;
			}	
		}
	}

	return 0;
}

