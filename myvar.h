#ifndef _MYVAR_H
#define _MYVAR_H
#include <stdio.h>
#include <math.h>
#include "sachead.h"
#define Ra 6371.0
#define Rearth 6371.0
#define PI M_PI
#define TWOPI   6.2831853071795847693
#define FOURPI 12.5663706143591729539
#define R2D    57.2957795130823208768
#define D2R  .017453292519943295769237

FILE *Fopen(char *filename,char *mode);
#ifdef __cplusplus
extern "C" {
int jul2cal(int jday,int yy,int* mm,int* dd);
int cal2jul(int mm,int dd,int yy);
int leap(int);
int delaz(double lat1,double lon1,double lat2,double lon2,double &gcarc,double &az);
void error();
}
#else
int jul2cal(int jday,int yy,int* mm,int* dd);
int cal2jul(int mm,int dd,int yy);
int delaz(double lat1,double lon1,double lat2,double lon2,
	  double *gcarc,double *az);
int leap(int);
void error();
#endif
#endif









