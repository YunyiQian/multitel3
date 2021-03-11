/*  Calculating Ray Information
 *  Return:
 *    0 if succeed
 *    -1 if failed
 *  Notice:
 *    this program only calculating reflected/transmited phases which combined with two parts
 *    such as PcP (P+P), sP(Sv+P) and so on.
 *  Example:
 *    sP -- P1 = 0	//indicating first part is S wave
 *          P2 = 1	//indicating second part is P wave
 *          Turn1 = 0	//indicating first part is no-turn wave
 *          Turn2 = 1	//indicating second part is turn wave
 *
 *    ScP -- P1 = 0	//first part is S wave
 *           P2 = 1	//second part is P wave
 *           Turn1 = 0	//no-turn
 *           Turn2 = 0	//no-turn
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.c"
#include "Reflection_Coefficient.c"

#define MAX 100000
#define CMB 2891.0
#define Rearth 6371.0
#define dS  1.0e4
#define deg2rad 0.0174533
#define rad2deg 57.29578
#define gcarc_min 0
#define gcarc_max 100

int rayinformation(
	char *name,		//name of input model, "none" means using default model(prem)
	double sdepth,		//depth of the source in km
	double ddepth,		//depth of the discontinuity in km
	double rdepth,		//depth of the receiver in km
	int indicator[4],	//indicators of ray type
	double Vd[6],		//information of discontinuity structure
	double Vr[2],		//information of receiver side structure, Vr[0] is Vs, Vr[1] is Vp
	char *outfile,		//name of output file
	double dh            //ifferential of hh when generating ray information
	){


	int P1 = indicator[0];	//indicator of first part of ray is P wave(1) or S wave(0)
	int P2 = indicator[1];	//indicator of second part of ray, same as above
	int Turn1 = indicator[2];	//indicator of first part of ray is going to turn(1) or not(0)
	int Turn2 = indicator[3];	//indicator of second part of ray, same as above

	int i,j;
        int kk;
	int ks,ke;	//number of source layer and discontinuity layer in model 1
	int ls,le;	//number of discontinuity layer and receiver layer in model 2
	int error1;
	int error2;
	int num;
	double hh;	//takeoff angle of phase
	double hstart,hend;	//start and end of hh
//	double dh = 0.005;	//differential of hh when generating ray information
	double parameter;	//ray parameter in spherical model
	double xc,ti;		//great circle distance and traveltime
	double Q;		//t star in attenuation

	double V1[2],V2[2];	//V1[0/1] is Vs/Vp of incident side and V2[0/1] is Vs/Vp of another side
	double Rho1,Rho2;	//Rho1 is rho of incident side and Rho2 is another side
	double dhh,dgcarc;	//differential of takeoff angle and great circle distance when calculating geometrical spreading factor
	double R_recv,R_disc;	//Riduals of receiver and discontinuity
	double factor;
	double ii,rr,ll;	//incident, refleced angle at discontinuity and incident angle at receiver
	double temp;
	int temp1;
        int temp2;
        int temp3;
        int temp4;
	int type;
	int sum;


	double HH[MAX],PA[MAX];		//takeoff angle and ray parameter
	double DEG[MAX],XC[MAX],XC1[MAX],XC2[MAX];	//great circle distance at discontinuity and at receiver
	double QQ[MAX],TI[MAX];		//t star in attenuation and traveltime
	double DIFFU[MAX];		//geometrical spreading factor
	double REF_CO[MAX];		//reflection coefficient factor
	MODEL *m1,*m2;
	FILE *fp;

	if ( Turn1 == 0 && sdepth > ddepth ){
		hstart = PI/2. + dh;
		hend   = PI - dh;
	}else{
//		hstart = dh;
//		hend   = PI/2. - dh;
		hstart = 0;
                hend   = PI/2. - dh;
	}
	
	for (i=0,sum=0;i<4;i++){
		sum += indicator[i];
	}

	m1 = (MODEL *)malloc(sizeof(MODEL));
	model_input(name,m1);
	model_init(m1,sdepth,&ks,ddepth,&ke,P1);

	m2 = (MODEL *)malloc(sizeof(MODEL));
	model_input(name,m2);
	model_init(m2,ddepth,&ls,rdepth,&le,P2);

	fp = fopen(outfile,"w");
	if (fp == NULL){
		fprintf(stderr,"Error in opening outfile\n");
		exit(-1);
	}

	R_recv = Rearth - rdepth;
	R_disc = Rearth - ddepth;

	for (hh = hstart,num = 0;hh < hend;hh += dh){
		parameter = m1->p[ks]*sin(hh);
		if (sum != 1){
			error1 = raytrace2(m1,parameter,ks,ke,Turn1,&xc,&Q,&ti);
		}else{
			error1 = raytrace3(m1,parameter,ks,ke,Turn1,&xc,&Q,&ti);
		}

		if (error1 == -1){
			continue;
		}

		HH[num] = hh;
		PA[num] = parameter;
		XC1[num]=xc;
		XC[num] = xc;
		DEG[num] = XC[num];
		TI[num] = ti;
		QQ[num] = Q;

		if (sum != 1){
			error2 = raytrace2(m2,parameter,ls,le,Turn2,&xc,&Q,&ti);
		}else{
			error2 = raytrace3(m2,parameter,ls,le,Turn2,&xc,&Q,&ti);
		}

		if (error2 == -1){
			continue;
		}
		XC2[num] =xc;
		XC[num] += xc;
		TI[num] += ti;
		QQ[num] += Q;

		if ( XC[num] < gcarc_min || XC[num] > gcarc_max || XC[num] != XC[num] ){
//                if ( XC[num] < gcarc_min || XC[num] > gcarc_max){
			continue;
		}

		num++;
	}

//	for (i=1;i<num-1;i++){
//		printf("dist1=%f\tdist2=%f\n",XC1[i],XC2[i]);
//		}	

	for (i=1;i<num-1;i++){
		/* determine which side the ray is approach to discontinuity */
		if (Turn1 == 0 && sdepth < ddepth){
			V1[0] = Vd[0];		//V1[0] is Vs of incident side
			V1[1] = Vd[1];		//V1[1] is Vp of incident side
			Rho1  = Vd[2];		//rho of incident side
			
			V2[0] = Vd[3];		//V2[0] is Vs of another side
			V2[1] = Vd[4];		//V2[1] is Vp of another side
			Rho2  = Vd[5];		//rho of another side
		}else{
			V1[0] = Vd[3];		//the same as describe above
			V1[1] = Vd[4];
			Rho1  = Vd[5];

			V2[0] = Vd[0];
			V2[1] = Vd[1];
			Rho2  = Vd[2];
		}

                dgcarc = (XC[i+1] - XC[i-1])/2    ;
                dhh  = (HH[i+1]-HH[i])        ;



		ii = asin( PA[i]*V1[P1]/R_disc );
		rr = asin( PA[i]*V1[P2]/R_disc );
		ll = asin( PA[i]*Vr[P2]/R_recv );

                temp = rad2deg*dhh*sin(HH[i])/( sin(deg2rad*XC[i])*cos(ll)*dgcarc );
		DIFFU[i] = sqrt( fabs(temp) )/R_recv;
// printf("gcarc=%f\ttheta=%f\ti0=%f\tsin(theta)=%f\tcos(i0)=%f\tsin(gcarc)=%f\tsin(theta)/sin(gcarc)/cos(i0)=%f\tdtheta=%f\tdgcarc=%f\tdtheta/dgcarc=%f\tG=%f\n",XC[i],HH[i]*180/PI,ll*180/PI,sin(HH[i]),cos(ll),sin(deg2rad*XC[i]),sin(HH[i])/(sin(deg2rad*XC[i])*cos(ll)),dhh*180/PI,dgcarc,dhh*180/PI/dgcarc,DIFFU[i]*1e5);
		/* determine type of reflection, no Sv to Sv reflection */
		if ( P1 == 1 && P2 == 1 ){
			type = 1;
		}else if ( P1 == 1 && P2 == 0 ){
			type = 2;
		}else if ( P1 == 0 && P2 == 1 ){
			type = 3;
		}else if ( P1 == 0 && P2 == 0 ){
			type = 5;
		}

		REF_CO[i] = Reflection_Coefficient(V1[1],V1[0],V2[1],V2[0],Rho1,Rho2,ii,type);
//		if ( REF_CO[i] != REF_CO[i] ) printf ("%lf\t%d\t%d\t%d\t%d\n",V1[0],V2[0],P2,Turn1,Turn2);
	}


                                               
	
//	DIFFU[1]=(17*DIFFU[i-4]+14*DIFFU[i-3]+11*DIFFU[i-2]+8*DIFFU[i-1]+5*DIFFU[i]+2*DIFFU[i+1]-DIFFU[i+2]-4*DIFFU[i+3]-7*DIFFU[i+4])/45;
//	DIFFU[2]=(56*DIFFU[i-4]+47*DIFFU[i-3]+38*DIFFU[i-2]+29*DIFFU[i-1]+20*DIFFU[i]+11*DIFFU[i+1]+2*DIFFU[i+2]-7*DIFFU[i+3]-16*DIFFU[i+4])/180;
//	DIFFU[3]=(22*DIFFU[i-4]+19*DIFFU[i-3]+16*DIFFU[i-2]+13*DIFFU[i-1]+10*DIFFU[i]+7*DIFFU[i+1]+4*DIFFU[i+2]+DIFFU[i+3]-2*DIFFU[i+4])/90;
//	DIFFU[4]=(32*DIFFU[i-4]+29*DIFFU[i-3]+26*DIFFU[i-2]+23*DIFFU[i-1]+20*DIFFU[i]+17*DIFFU[i+1]+14*DIFFU[i+2]+11*DIFFU[i+3]+8*DIFFU[i+4])/180;
	for (i=25;i<num-16;i++){
        	DIFFU[i]=(DIFFU[i-4]+DIFFU[i-3]+DIFFU[i-2]+DIFFU[i-1]+DIFFU[i]+DIFFU[i+1]+DIFFU[i+2]+DIFFU[i+3]+DIFFU[i+4])/9;
        }//9 for all

        for (j=1;j<30;j++){
        	for (i=1;i<num-1;i++){
			if(XC[i]>=26 && XC[i]<=40)
                        	DIFFU[i]=(DIFFU[i-4]+DIFFU[i-3]+DIFFU[i-2]+DIFFU[i-1]+DIFFU[i]+DIFFU[i+1]+DIFFU[i+2]+DIFFU[i+3]+DIFFU[i+4])/9;
        	}
        }
        for (j=1;j<30;j++){
                for (i=1;i<num-1;i++){
                        if(XC[i]>=35 && XC[i]<=50)
                                DIFFU[i]=(DIFFU[i-4]+DIFFU[i-3]+DIFFU[i-2]+DIFFU[i-1]+DIFFU[i]+DIFFU[i+1]+DIFFU[i+2]+DIFFU[i+3]+DIFFU[i+4])/9;
                }
        }

        for (j=1;j<30;j++){
                for (i=1;i<num-1;i++){
                        if(XC[i]>=85 && XC[i]<=95)
                                DIFFU[i]=(DIFFU[i-4]+DIFFU[i-3]+DIFFU[i-2]+DIFFU[i-1]+DIFFU[i]+DIFFU[i+1]+DIFFU[i+2]+DIFFU[i+3]+DIFFU[i+4])/9;
                }
        }
        for (j=1;j<30;j++){
                for (i=1;i<num-1;i++){
                        if(XC[i]>=80 && XC[i]<=90)
                                DIFFU[i]=(DIFFU[i-4]+DIFFU[i-3]+DIFFU[i-2]+DIFFU[i-1]+DIFFU[i]+DIFFU[i+1]+DIFFU[i+2]+DIFFU[i+3]+DIFFU[i+4])/9;
                }
        }

//	DIFFU[num-5]=(8*DIFFU[i-4]+11*DIFFU[i-3]+14*DIFFU[i-2]+17*DIFFU[i-1]+20*DIFFU[i]+23*DIFFU[i+1]+26*DIFFU[i+2]+29*DIFFU[i+3]+32*DIFFU[i+4])/180;
//        DIFFU[num-4]=(-2*DIFFU[i-4]+DIFFU[i-3]+4*DIFFU[i-2]+7*DIFFU[i-1]+10*DIFFU[i]+13*DIFFU[i+1]+16*DIFFU[i+2]+19*DIFFU[i+3]+22*DIFFU[i+4])/90;
//        DIFFU[num-3]=(-16*DIFFU[i-4]-7*DIFFU[i-3]+2*DIFFU[i-2]+11*DIFFU[i-1]+20*DIFFU[i]+29*DIFFU[i+1]+38*DIFFU[i+2]+47*DIFFU[i+3]+56*DIFFU[i+4])/180;
//        DIFFU[num-2]=(-7*DIFFU[i-4]-4*DIFFU[i-3]-DIFFU[i-2]+2*DIFFU[i-1]+5*DIFFU[i]+8*DIFFU[i+1]+11*DIFFU[i+2]+14*DIFFU[i+3]+17*DIFFU[i+4])/45;

	fprintf(fp,"%d\n",num-2);
	for (i=1;i<num-1;i++){
		fprintf(fp,"%f\t%f\t%e\t%f\t%f\t%f\t%f\n",XC[i],PA[i],DIFFU[i]*1.e5,TI[i],QQ[i],REF_CO[i],HH[i]*rad2deg);
//		fprintf(fp,"%f\t%f\t%f\n",HH[i],TI[i],XC[i]);
	}
	fclose(fp);
	fp = NULL;
	free (m1);
	free (m2);

	return 0;
}

