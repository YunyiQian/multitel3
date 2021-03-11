#include <stdio.h>
#include "myvar.h"
#include "raytrace.h"
/* Return:
 *	-1 if failed
 *	0 if succeed
 */
#define BIG 10000  //bigger than max S wave parameter
#define SMALL 0.1
#define HZ 1.0
/* Usage:
 *    #include "raytrace.h"
 *    model_input(filename,model)
 *  IN:
 *    filename -- name of model file
 *  OUT:
 *    model->nlay,model->depth[],vp[],vs[],rou[],qp[],qs[]
 *  RETURN:
 *    number of layers in model file if succeed
 *    -1 if failed
 */

int model_input(char *name,MODEL *m){
	FILE *fp;
	int error;
	int i,im1;

	if ( name==NULL ){
		fprintf(stderr,"No model name\n");
		return -1;
	}else{
		fp = fopen(name,"r");
		if ( fp == NULL ){
			fprintf(stderr,"Error in opening model file\n");
			return -1;
		}
	}

	fscanf(fp,"%d",&m->nlay);
	if ( m->nlay > MAXLAY ){
		fprintf(stderr,"Too many layers in model file\n");
		m->nlay = MAXLAY;
	}

	for (i=0;i<m->nlay;i++){
		fscanf(fp,"%lf%lf%lf%lf%lf%lf",&m->depth[i],&m->vp[i],&m->vs[i],&m->rou[i],&m->qp[i],&m->qs[i]);
	}
	
	fclose(fp);

	for (i=1,im1=0;i<m->nlay;i++,im1++){
		if (m->depth[i] < m->depth[im1]){
			fprintf(stderr,"Depth in model file is not correct\n");
			return -1;
		}
	}

	return m->nlay;
}

/* Usage:
 *    #include "raytrace.h"
 *    model_init(m,sdepth,ilay1,rdepth,ilay2,P)
 *  IN:
 *    m -- pointer of MODEL m
 *    sdepth -- source depth of layer we will insert
 *    ilay1 -- number of source depth
 *    rdepth -- receiver depth of layer we will insert
 *    ilay2 -- number of receiver depth
 *    P -- P wave(1) or S wave(0) raytrace
 *  OUT:
 *    ilay1 & ilay2
 */

void model_init(MODEL *m,double sdepth,int *ilay1,double rdepth,int *ilay2,int P){
	int i;
	
	if (sdepth < rdepth){
		*ilay1 = inslay(sdepth,m);
		*ilay2 = inslay(rdepth,m);
	}else{
		*ilay2 = inslay(rdepth,m);
		*ilay1 = inslay(sdepth,m);
	}
	if ( *ilay1 == -1 || *ilay2 == -1 ){
		fprintf(stderr,"Error in creating layer at given depth\n");
		exit (-1);
	}

	mhrv_interp(m,P);
	model_setcmb(m);
	return;
}

/* Usage:
 *    #include "raytrace.h"
 *    inslay(depth,m)
 *  IN:
 *    depth -- depth of layer that will be inserted in model m
 *    m -- model
 *  RETURN:
 *    inslay -- the number of layer where we inserting in model if succeed
 *    -1 if failed
 *
 */

int inslay(double depth,MODEL *m){
	int i,inslay;
	int k;
	double slope;

	if (depth < 0. || depth > m->depth[m->nlay-1]){
		fprintf(stderr,"Input depth is out of range\n");
		return -1;
	}

	for (i=0;i<m->nlay;i++){
		if (m->depth[i] > depth){
			if ( fabs(depth-m->depth[i-1]) < 0.005 ){
				inslay = i-1;
				return inslay;
			}else if( fabs(depth-m->depth[i]) < 0.005 ){
				inslay = i;
				return inslay;
			}else{
				inslay = i;
				m->nlay ++;
				if (m->nlay > MAXLAY){
					fprintf(stderr,"Too many layers in inslay()\n");
					return -1;
				}
				for (k=(m->nlay-1);k>i;k--){
					m->depth[k] = m->depth[k-1];
					m->vp[k]    = m->vp[k-1];
					m->vs[k]    = m->vs[k-1];
					m->rou[k]   = m->rou[k-1];
					m->qp[k]    = m->qp[k-1];
					m->qs[k]    = m->qs[k-1];
				}

				slope = (depth-m->depth[i-1])/(m->depth[i+1] - m->depth[i-1]);
				m->depth[i] = depth;
				m->vp[i]  = slope*m->vp[i+1]  + (1.-slope)*m->vp[i-1];
				m->vs[i]  = slope*m->vs[i+1]  + (1.-slope)*m->vs[i-1];
				m->rou[i] = slope*m->rou[i+1] + (1.-slope)*m->rou[i-1];
				m->qp[i]  = slope*m->qp[i+1]  + (1.-slope)*m->qp[i-1];
				m->qs[i]  = slope*m->qs[i+1]  + (1.-slope)*m->qs[i-1];

				return inslay;
			}	// end of if
		}// end of if
	}// end of i loop 

}

/* Usage:
 *    #include "raytrace.h"
 *    mhrv_interp(m,P)
 *  IN:
 *    m -- model
 *    P -- indicator of P wave(1) or S wave(other)
 *  MODIFIED:
 *    m.radius[l] -- radius of lth layer
 *    m.p[l] -- turn ray parameter (spherical) of lth layer
 *    m.B[l] -- interpolation character between lth and (l+1)th layers, Velo=Ar**B
 *      if at lth layer, vs == 0, then vp is used instead
 */

void mhrv_interp(MODEL *m,int P){
	int i;
	int k=0;
	double tmp,tmp1;
	
	m->radius[0] = Rearth - m->depth[0];
	if ( P == 1 || m->vs[0] < SMALL ){
		m->p[0] = m->radius[0]/m->vp[0];
		m->Q[0] = m->qp[0];
	}else{
		m->p[0] = m->radius[0]/m->vs[0];
		m->Q[0] = m->qs[0];
	}
/*
	if ( P == 1 ){
		m->p[0] = m->radius[0]/m->vp[0];
		m->Q[0] = m->qp[0];
	}else if ( m->vs[0] < SMALL ){
		m->p[0] = m->radius[0]/SMALL;
		m->Q[0] = SMALL;
	}else{
		m->p[0] = m->radius[0]/m->vs[0];
		m->Q[0] = m->qs[0];
	}
*/
	for (i=1;i<m->nlay;i++){
		m->radius[i] = Rearth - m->depth[i];
		if (m->radius[i] < 0.001){	// the center of earth
			m->p[i]   = 0.00001;
			m->Q[i]   = 100000;
			m->B[i-1] = 1.;
			m->B[i]   = 1.;
			return;
		}else{
			tmp = log(m->radius[i-1]/m->radius[i]);
		}

		if ( P == 1 ){
			m->p[i] = m->radius[i]/m->vp[i];
			m->Q[i] = m->qp[i];
			tmp1 = log(m->vp[i-1]/m->vp[i]);
		}else if( m->vs[i] < SMALL ){
			m->Q[i] = m->qp[i];
			if ( k == 1 ){
				m->p[i] = m->radius[i]/m->vp[i];
				tmp1 = log(m->vp[i-1]/m->vp[i]);
			}else{
//				m->p[i] = m->radius[i]/m->vp[i];
				m->p[i] = BIG;
				tmp1 = log(m->vs[i-1]/m->vp[i]);
				k = 1;
			}
		}else{
			m->p[i] = m->radius[i]/m->vs[i];
			m->Q[i] = m->qs[i];
			tmp1 = log(m->vs[i-1]/m->vs[i]);
		}

/*		if ( P == 1 ){
			m->p[i] = m->radius[i]/m->vp[i];
			m->Q[i] = m->qp[i];
			tmp1 = log(m->vp[i-1]/m->vp[i]);
		}else{
			if ( m->vs[i] < SMALL ){
				m->vs[i] = SMALL;
			}
			m->p[i] = m->radius[i]/m->vs[i];
			m->Q[i] = m->qs[i];
			tmp1 = log(m->vs[i-1]/m->vs[i]);
		}
*/
		if (fabs(tmp) < 1.0e-4){	//discontinuity
			m->B[i-1] = 0.;
		}else{
			tmp = tmp1/tmp;
			m->B[i-1] = 1./(1.-tmp);
		}
	}
	m->B[i] = 1.;

	return;
}

/* Usage:
 *    #include "raytrace.h"
 *    model_setcmb(m)
 *  IN:
 *    m -- model
 *  MODIFIED:
 *    m->ncmb -- the number of core-mantle boundary
 *    m->nicb -- the number of inner-outercore boundary
 */

void model_setcmb(MODEL *m){
	int i;

	m->ncmb = -1;
	m->nicb = -1;
	for (i=1;i<m->nlay;i++){
		if (m->vs[i-1] < SMALL && m->vs[i] > SMALL){
			m->nicb = i;
			return;
		}
		if (m->vs[i-1] > SMALL && m->vs[i] < SMALL){
			m->ncmb = i;
		}
	}
	return;
}

/* Usage:
 *    #include "raytrace.h"
 *    raytrace(m,p,ks,ke,TURN,xc[],yc[],ti[])
 *  PURPOSE:
 *    1-D raytrace between layer ks, ke, for spherical rayparameter p
 *    path positions are returned in xc[],yc[], time in ti[]
 *  IN:
 *    m -- model
 *    p -- spherical ray parameter in sec (double)
 *    ks, ke -- start and end layer (int)
 *    TURN -- indicator of whether ray is going to turn
 *  OUT:
 *    xc[l] -- great circle distance between ks and lth step's location in degree
 *    yc[l] -- depth of corresponding step
 *    ti[l] -- travel time of corresponding step
 *  RETURN:
 *    number of points in the ray path
 *    -1 if error occuring
 *  REVISION:
 *    9/23/93 Xiaoming Ding, first version
 *    3/10/96 by Xiaoming Ding,  get rid of global variable,and put in one package, raytrace.h and raytrace.c
 *
 *  NOTICE:
 *    In this program, first all points are computed. That's some waste, but not much.
 */

int raytrace(MODEL *m,double p,int ks,int ke,int TURN,double xc[],double yc[],double ti[]){

	int i,nturn,npt;
	int l0,il;
	int l;
	int k0;
	double w;
	double pdiv0,pdiv1;
	double psq0,psq1;
	double xsum;
	double tsum;
	double xtmp;
	double ttmp;
	double p2;

	if (ks < 0 || ke < 0 || p < 0.){
		fprintf(stderr,"Input is illegal\n");
		return -1;
	}

	for (i=0;i<m->nlay;i++){
		if (m->p[i] <= p){
			nturn = i;
			break;
		}
	}

	if (ks > nturn || nturn == 0 || ke > nturn){
		fprintf(stderr,"The wave cannot arrive the ke layer\n");
		return -1;
	}

	l0 = ks < ke ? ks : ke;
	p2 = p*p;

	xsum  = 0.;
	tsum  = 0.;
	xc[0] = 0.;
	ti[0] = 0.;
	yc[0] = m->depth[l0];

	pdiv0 = acos(p/m->p[l0]);
	psq0  = sqrt(m->p[l0]*m->p[l0]-p2);

	for (l=l0+1,i=1;l < nturn;l++,i++){
		w = m->B[l-1];
		pdiv1 = acos(p/m->p[l]);
		psq1  = sqrt(m->p[l]*m->p[l] - p2);
		xsum += w*(pdiv0 - pdiv1);
		tsum += w*(psq0 - psq1);
		yc[i] = m->depth[l];
		xc[i] = xsum*R2D;
		ti[i] = tsum;

		pdiv0 = pdiv1;
		psq0  = psq1;
	}

	w = m->B[l-1];
	xsum += w*pdiv0;
	tsum += w*psq0;
	xsum  = xc[i] = xsum*R2D;
	yc[i] = Rearth - m->radius[l-1]*pow(p/m->p[l-1],w);
	ti[i] = tsum;
	
	for (l=i-1;l>=0;l--,i++){
		xc[i] = 2.*xsum - xc[l];
		ti[i] = 2.*tsum - ti[l];
		yc[i] = yc[l];
	}

	if (TURN == 1){
		npt = 2*nturn + 1 - ks - ke;
	}else{
		if ( ks <= ke ){
			npt = ke - ks + 1;
		}else{
			npt = ks - ke + 1;
		}
	}

	if ( ks <= ke ){
		return npt;
	}else{
		if (TURN == 1){
			k0 = ks - ke;
		}else{
			k0 = 2*nturn - ks - ke;
		}

		xtmp = xc[k0];
		ttmp = ti[k0];
		for (i=0,il=k0;i < npt;i++,il++){
			xc[i] = xc[il] - xtmp;
			ti[i] = ti[il] - ttmp;
			yc[i] = yc[il];
		}

		return npt;
	}
}

/* Usage:
 *    #include "raytrace.h"
 *    raytrace2(m,p,ks,ke,TURN,&xc,&tt,&ti)
 *  IN:
 *    m -- model
 *    p -- spherical ray parameter in sec (double)
 *    ks,ke -- start and end layer (int)
 *    TURN -- indicatation of whether the ray is going to turn
 *    xc -- great circle distance between ks and ke layer
 *    tt -- t star in attenuation calculation
 *    ti -- travel time
 */

int raytrace2(MODEL *m,double p,int ks,int ke,int TURN,double *xc,double *tt,double *ti){

	int i,k;
	int l0,l1;
	int nturn;
	double w;
	double pdiv0,pdiv1;
	double psq0,psq1;
	double xsum,xtmp;
	double tsum,ttmp;
	double dt;
	double tstar,temp;
	double p2;

	if (ks < 0 || p < 0 || ke < 0){
		fprintf(stderr,"Input is illegal\n");
		return -1;
	}

	for (i=0,nturn=0;i < m->nlay;i++){
		if (m->p[i] <= p){
			nturn = i;
			break;
		}
	}

	if (p < 1.4 || nturn == 0){
		nturn = m->nlay;
	}

	if (ks >= nturn || nturn == 0 || ke >= nturn){
//if (nturn > 45)
//printf("ks=%d\tke=%d\tnturn=%d\n",ks,ke,nturn);
//		fprintf(stderr,"The wave cannot arrives the ke layer\n");
		return -1;
	}

	if ( ks > ke ){
		l0 = ke;
		l1 = ks;
	}else{
		l0 = ks;
		l1 = ke;
	}

	p2 = p*p;
	xsum  = 0.;
	tsum  = 0.;
	tstar = 0.;
	pdiv0 = acos(p/m->p[l0]);
	psq0  = sqrt(m->p[l0]*m->p[l0] - p2);

	for (i=l0+1;i<=l1;i++){
		w = m->B[i-1];
		pdiv1  = acos(p/m->p[i]);
		psq1   = sqrt(m->p[i]*m->p[i] - p2);
		xsum  += w*(pdiv0 - pdiv1);
		dt     = w*(psq0 - psq1);
		tsum  += dt;
		if ( m->Q[i-1] < SMALL ){
//			printf ("Bad Q value in raytrace2\n");
			return -1;
		}else{
			tstar += dt/m->Q[i-1];
		}
		pdiv0  = pdiv1;
		psq0   = psq1;
	}

	if (TURN == 1){
		xtmp = xsum;
		ttmp = tsum;
		temp = tstar;
		xsum = 0.;
		tsum = 0.;
		tstar = 0.;

		for (i=l1+1;i<nturn;i++){
			w = m->B[i-1];
			pdiv1  = acos(p/m->p[i]);
			psq1   = sqrt(m->p[i]*m->p[i] - p2);
			xsum  += w*(pdiv0 - pdiv1);
			dt     = w*(psq0 - psq1);
			tsum  += dt;
			if ( m->Q[i-1] < SMALL ){
//				printf ("Bad Q value in raytrace2\n");
				return -1;
			}else{
				tstar += dt/m->Q[i-1];
			}
			pdiv0  = pdiv1;
			psq0   = psq1;
		}
		w = m->B[nturn-1];
		xsum  += w*pdiv0;
		dt     = w*psq0;
		tsum  += dt;
		tstar += dt/m->Q[nturn-1];
		
		*xc = R2D*(2.*xsum + xtmp);
		*ti = 2.*tsum + ttmp;
//		*bot = Rearth - m->radius[nturn-1]*pow(p/m->p[nturn-1],w);
		*tt = 2.*tstar + temp;
	}else{
		*xc = xsum*R2D;
		*ti = tsum;
//		*bot = m->depth[l1];
		*tt = tstar;
	}
	return 0;
}


/* Usage:
 *    #include "raytrace.h"
 *    raytrace3(m,p,ks,ke,TURN,&xc,&tt,&ti)
 * Note:
 *    only used to calculate ray traveling in mantle
 *
 *  IN:
 *    m -- model
 *    p -- spherical ray parameter in sec (double)
 *    ks,ke -- start and end layer (int)
 *    TURN -- indicatation of whether the ray is going to turn
 *    xc -- great circle distance between ks and ke layer
 *    tt -- t star in attenuation calculation
 *    ti -- travel time
 *    bot -- bottle of ray
 */

int raytrace3(MODEL *m,double p,int ks,int ke,int TURN,double *xc,double *tt,double *ti){

	int i,k;
	int l0,l1;
	int nturn;
	double bot;
	double w;
	double pdiv0,pdiv1;
	double psq0,psq1;
	double xsum,xtmp;
	double tsum,ttmp;
	double dt;
	double tstar,temp;
	double p2;

	if (ks < 0 || p < 0 || ke < 0){
		fprintf(stderr,"Input is illegal\n");
		return -1;
	}

	for (i=0,nturn=0;i < m->nlay;i++){
		if (m->p[i] <= p){
			nturn = i;
			break;
		}
	}

	if (p < 1.4 || nturn == 0){
		nturn = m->nlay;
	}

	if (nturn >= m->ncmb){
//		printf("nturn=%d\tcmb=%d\n",nturn,m->ncmb);
		return -1;
	}

	if (ks >= nturn || nturn == 0 || ke >= nturn){
//if (nturn > 45)
//printf("ks=%d\tke=%d\tnturn=%d\n",ks,ke,nturn);
//		fprintf(stderr,"The wave cannot arrives the ke layer\n");
		return -1;
	}

	if ( ks > ke ){
		l0 = ke;
		l1 = ks;
	}else{
		l0 = ks;
		l1 = ke;
	}

	p2 = p*p;
	xsum  = 0.;
	tsum  = 0.;
	tstar = 0.;
	pdiv0 = acos(p/m->p[l0]);
	psq0  = sqrt(m->p[l0]*m->p[l0] - p2);

	for (i=l0+1;i<=l1;i++){
		w = m->B[i-1];
		pdiv1  = acos(p/m->p[i]);
		psq1   = sqrt(m->p[i]*m->p[i] - p2);
		xsum  += w*(pdiv0 - pdiv1);
		dt     = w*(psq0 - psq1);
		tsum  += dt;
		if ( m->Q[i-1] < SMALL ){
//			printf ("Bad Q value in raytrace2\n");
			return -1;
		}else{
			tstar += dt/m->Q[i-1];
		}
		pdiv0  = pdiv1;
		psq0   = psq1;
	}

	if (TURN == 1){
		xtmp = xsum;
		ttmp = tsum;
		temp = tstar;
		xsum = 0.;
		tsum = 0.;
		tstar = 0.;

		for (i=l1+1;i<nturn;i++){
			w = m->B[i-1];
			pdiv1  = acos(p/m->p[i]);
			psq1   = sqrt(m->p[i]*m->p[i] - p2);
			xsum  += w*(pdiv0 - pdiv1);
			dt     = w*(psq0 - psq1);
			tsum  += dt;
			if ( m->Q[i-1] < SMALL ){
//				printf ("Bad Q value in raytrace2\n");
				return -1;
			}else{
				tstar += dt/m->Q[i-1];
			}
			pdiv0  = pdiv1;
			psq0   = psq1;
		}
		w = m->B[nturn-1];
		xsum  += w*pdiv0;
		dt     = w*psq0;
		tsum  += dt;
		tstar += dt/m->Q[nturn-1];
		
		*xc = R2D*(2.*xsum + xtmp);
		*ti = 2.*tsum + ttmp;
		bot = Rearth - m->radius[nturn-1]*pow(p/m->p[nturn-1],w);
		*tt = 2.*tstar + temp;
	}else{
		*xc = xsum*R2D;
		*ti = tsum;
		bot = m->depth[l1];
		*tt = tstar;
	}
/*	if ( (bot - m->depth[m->ncmb]) < -SMALL ){
//printf("bot=%f\tcmb=%f\n",bot,m->depth[m->ncmb]);
		return -1;
	}
*/
	return 0;
}

