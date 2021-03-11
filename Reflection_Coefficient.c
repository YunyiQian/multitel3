//***********************************//
//Usage:calculate the coefficient of reflection
//      according to SEISMIC RAY THEORY(Page 481)
//
//  type:
//	1 - return coefficient of P to P reflection
//      2 - return coefficient of P to SV reflection
//	3 - return coefficient of SV to P reflection
//	4 - return coefficient of SV to SV reflection
//	5 - return coefficient of SH to SH reflection
//
//***********************************//

#include <stdlib.h>
#include <math.h>

double Reflection_Coefficient(
	double c1_p,	//Vp of media 1(incident side)
	double c1_s,	//Vs of media 1
	double c2_p,	//Vp of media 2
	double c2_s,	//Vs of media 2
	double rho1,	//density of media 1
	double rho2,	//density of media 2
	double i1,	//incidence of plane wave
	int type	//indication of reflection mode
	){

	double p;
        double Coefficient;
	double D;
	double P1,P2,P3,P4;
	double q,X,Y,Z;

	if (type > 5 || type < 1){
		fprintf(stderr,"Error type is given\n");
		exit(-1);
	}

	if (type == 1 || type == 2){
		p = sin(i1)/c1_p;
	}else{
		p = sin(i1)/c1_s;
	}

	P1 = sqrt(1 - c1_p*c1_p*p*p);
	P2 = sqrt(1 - c1_s*c1_s*p*p);

	P3 = sqrt(1 - c2_p*c2_p*p*p);
	P4 = sqrt(1 - c2_s*c2_s*p*p);

	q = 2.*(rho2*c2_s*c2_s - rho1*c1_s*c1_s);

	X = rho2 - q*p*p;
	Y = rho1 + q*p*p;
	Z = rho2 - rho1 - q*p*p;

	D =  q*q*p*p*P1*P2*P3*P4 
	   + rho1*rho2*(c1_s*c2_p*P1*P4 + c1_p*c2_s*P2*P3)
	   + c1_p*c1_s*P3*P4*Y*Y
	   + c2_p*c2_s*P1*P2*X*X
	   + c1_p*c2_p*c1_s*c2_s*p*p*Z*Z;

	switch (type){
	case 1:
		Coefficient =  (q*q*p*p*P1*P2*P3*P4 
			     + rho1*rho2*(c1_s*c2_p*P1*P4 - c1_p*c2_s*P2*P3)
			     - c1_p*c1_s*P3*P4*Y*Y
			     + c2_p*c2_s*P1*P2*X*X
			     - c1_p*c2_p*c1_s*c2_s*p*p*Z*Z)/D;
		break;
	case 2:
		Coefficient = 2.*c1_p*p*P1*(q*P3*P4*Y + c2_p*c2_s*X*Z)/D;
		break;
	case 3:
		Coefficient = -2.*c1_s*p*P2*(q*P3*P4*Y + c2_p*c2_s*X*Z)/D;
		break;
	case 4:
		Coefficient =  (q*q*p*p*P1*P2*P3*P4
			     + rho1*rho2*(c1_p*c2_s*P2*P3 - c1_s*c2_p*P1*P4)
			     - c1_p*c1_s*P3*P4*Y*Y
			     + c2_p*c2_s*P1*P2*X*X
			     - c1_p*c2_p*c1_s*c2_s*p*p*Z*Z)/D;
		break;
	case 5:
		Coefficient = (rho1*c1_s*P2 - rho2*c2_s*P4)/(rho1*c1_s*P2 + rho2*c2_s*P4);
		break;
	}

	return Coefficient;
}
