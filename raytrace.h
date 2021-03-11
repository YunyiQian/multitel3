#include "boolean.h"
#define MAXLAY 1000
typedef struct MODEL_1D{	/* Model is continuous model, not layered model*/
  int nlay;			/* number of layers */

  double vp[MAXLAY],vs[MAXLAY];	/* parameter in layer */
  double rou[MAXLAY],depth[MAXLAY];
  double qp[MAXLAY],qs[MAXLAY];

  double radius[MAXLAY];	/* derived parameters */
  double B[MAXLAY],p[MAXLAY];
  double Q[MAXLAY];
  int ncmb,nicb;
} MODEL;

int model_input(char *name,MODEL *m);
void model_init(MODEL *m,double sdepth,int *ilay1,double rdepth,int *ilay2,int P);
int inslay(double depth,MODEL *m);
void model_setcmb(MODEL *m);
void mhrv_interp(MODEL *m,int P); // Mohoroviviv's interpolate 
int raytrace(MODEL *m,double p,int ks,int ke,int TURN,
		double xc[],double yc[],double ti[]);
int raytrace2(MODEL *m,double p,int ks,int ke,int TURN, 
		double *xc,double *yc,double *tc);
int raytrace3(MODEL *m,double p,int ks,int ke,int TURN, 
		double *xc,double *tc,double *Q);

