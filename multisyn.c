#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"
#include "Complex.h"
#include "radiats.h"

int main(int argc, char **argv) {
  int	i,j,k,nn,ns,npt[3],npt1,npt2,error,intg,diff,src_type,filter,pt_shift;
  char	nam1[128],nam2[128],outnm[128],*ccc,*ddd,com[3]={'z','r','t'};
  int	dynamic=1,core=0;
  float	coef,rad[4][3],m0,az,*grn1,*grn2,*syn[3],*src,*pt,disp[3],mt[3][3];
  float cmpinc[3]={0,90.,90.}, cmpaz[3];
  float dt, dura, rise, tp, ts, tmp, dist, shift;
  float	*trapezoid(float, float, float, int *);
  float b[3],b1[3],b2[3],e[3],e1[3],e2[3];
  SACHEAD	hd,hd1,hd2;
#ifdef SAC_LIB
  char type[2] = {'B','P'}, proto[2] = {'B','U'};
  float	sn[30], sd[30];
  double f1, f2;
  long int order, nsects;
#endif
  void  fttq_(float *, float *, int *, int *, float *);
  int	mftm=2048, nftm;
  float	tstar=0., ftm[2048];

  /* input parameters */
  ns = 0;
  dura = 0.;
  src_type=0;
  intg=0;
  diff=0;
  filter=0;
  shift=0.;
  error = 0;
  for (i=1; !error && i < argc; i++) {
     if (argv[i][0] == '-') {
	switch(argv[i][1]) {

 	   case 'A':
	      sscanf(&argv[i][2], "%f",&az);
	      cmpaz[0] = 0.;
	      cmpaz[1] = az;
	      cmpaz[2] = az+90.;
	      if (cmpaz[2] > 360.) cmpaz[2] -= 360.;
	      break;
      
      case 'C':
         strcpy(nam2,&argv[i][2]);
         core = 1;
	      break;

 	   case 'D':
	      j = sscanf(&argv[i][2], "%f/%f",&dura,&rise);
	      if (j<2) rise = 0.5;
	      break;

#ifdef SAC_LIB
 	   case 'F':
	      filter = 1;
	      j = sscanf(&argv[i][2], "%lf/%lf/%ld",&f1,&f2,&order);
	      if (j<3) order = 4;
	      break;
#endif

 	   case 'G':
	      strcpy(nam1,&argv[i][2]);
	      break;

	   case 'I':
	      intg = 1;
	      break;

	   case 'J':
	      diff = 1;
	      break;

 	   case 'M':
	      src_type = sscanf(&argv[i][2], "%f/%f/%f/%f/%f/%f/%f",&m0,&mt[0][0],&mt[0][1],&mt[0][2],&mt[1][1],&mt[1][2],&mt[2][2]);
	      break;

 	   case 'O':
	      strcpy(outnm, &argv[i][2]);
	      break;

	   case 'P':
	      dynamic = 0;
	      break;

	   case 'Q':
	      sscanf(&argv[i][2], "%f",&tstar);
	      break;

 	   case 'S':
	      if ( (src=read_sac(&argv[i][2],&hd)) != NULL ) {
                 ns = hd.npts;
		 shift = -hd.b;
	      }
	      break;

	   default:
	      error = 1;
	      break;
	}
     }
        else error = 1;
  }

  switch (src_type) {
  case 1:
     nn = 1;
     m0 = m0*1.0e-20;
     break;
  case 3:
     nn = 2;
     sf_radiat(az-mt[0][0],mt[0][1],rad);
     m0 = m0*1.0e-15;
     break;
  case 4:
     nn = 3;
     dc_radiat(az-mt[0][0],mt[0][1],mt[0][2],rad);
     m0 = pow(10.,1.5*m0+16.1-20);
     break;
  case 7:
     nn = 4;
     mt_radiat(az,mt,rad);
     m0 = m0*1.0e-20;
     break;
  default:
     error = 1;
  }

  if (dynamic) ccc = strrchr(nam1, (int) '.') + 1;
  if (core) ddd = strrchr(nam2, (int) '.') + 1;

  if(argc == 1 || error || (dynamic && src_type == 1 && (*ccc) != 'a') ) {
     fprintf(stderr,"Usage: %s -Mmag([[/Strike/Dip]/Rake]|/Mxx/Mxy/Mxz/Myy/Myz/Mzz) -Aazimuth ([-SsrcFunctionName | -Ddura[/rise]] [-Ff1/f2[/n]] [-I | -J] -OoutName.z -GFirstCompOfGreen | -P)\n\
   Compute displacements in cm in the up, radial, and transverse (clockwise) directions produced by difference seismic sources\n\
   -M Specify source magnitude and orientation or moment-tensor\n\
      For double-couple, mag is Mw, strike/dip/rake are in A&R convention\n\
      For explosion; mag in in dyne-cm, no strike, dip, and rake needed\n\
      For single-force source; mag is in dyne, only strike and dip are needed\n\
      For moment-tensor; mag in dyne-cm, x=N,y=E,z=Down\n\
   -A Set station azimuth in degree measured from the North\n\
   -C Give the name of the first component of the core FK Green function\n\
   -D Specify the source time function as a trapezoid,\n\
      give the total duration and rise-time (0-0.5, default 0.5=triangle)\n\
   -F apply n-th order Butterworth band-pass filter, SAC lib required (off, n=4, must be < 10)\n\
   -G Give the name of the first component of the direct FK Green function\n\
   -I Integration once\n\
   -J Differentiate the synthetics\n\
   -O Output SAC file name\n\
   -P Compute static displacement, input Green functions from stdin in the form\n\
	distance Z45 R45 T45 ZDD RDD TDD ZSS RSS TSS [distance ZEX REX TEX]\n\
      The displacements will be output to stdout in the form of\n\
	distance azimuth z r t\n\
   -Q Convolve a Futterman Q operator of tstar (no)\n\
   -S Specify the SAC file name of the source time function (its sum. must be 1)\n\
   Examples:\n\
   * To compute three-component velocity at N33.5E azimuth from a Mw 4.5\n\
earthquake (strike 355, dip 80, rake -70), use:\n\
	multisyn -M4.5/355/80/-70 -D1 -A33.5 -OPAS.z -Gsp6direct_10/50.grn.0 -Csp6core_10/50.grn.0\n\
   * To compute the static displacements from the same earthquake, use:\n\
	nawk \'$1==50\' st.out | multisyn -M4.5/355/80/-70 -A33.5 -P\n\
   * To compute displacement from an explosion, use:\n\
   	multisyn -M3.3e20 -D1 -A33.5 -OPAS.z -Gsp6direct_10/50.grn.a -Csp6core_10/50.grn.a\n\
      or\n\
        multisyn -M3.3e20/1/0/0/1/0/1 -D1 -A33.5 -OPAS.z -Gsp6direct_10/50.grn.0 -Csp6core_10/50.grn.0\n\
	\n",argv[0]);
    return -1;
  }

  for(j=0; j<3; j++) disp[j] = 0.;
  for(i=0; i<nn; i++) {	/* sum up contribution from a DD, DS, SS, and EX */
//   for(i=0; i<1; i++) {	/* sum up contribution from a DD, DS, SS, and EX */
     for(j=0; j<3; j++) {
	coef = m0; if (nn>1) coef *= rad[i][j];
	if (!dynamic) {
	   if ((i==0||i==3) && j==0) scanf("%f",&dist);
	   scanf("%f",&tmp);
	   disp[j] += coef*tmp;
	   continue;
	}
        if ( (grn1=read_sac(nam1,&hd1)) == NULL ) continue;
        if (core) {
           if ( (grn2=read_sac(nam2,&hd2)) == NULL ) continue;
        }

	if ( i==0 ) {
      hd = hd1;
      dt = hd1.delta;
      npt1 = hd1.npts;
      npt[j] = npt1;
      b[j] = hd1.b;
      e[j] = hd1.e;
      if (core) {
         b1[j] = hd1.b;
         e1[j] = hd1.e;
         b2[j] = hd2.b;
         e2[j] = hd2.e;
         npt2 = hd2.npts;
         if ( b1[j]<=b2[j] ) b[j] = b1[j];
         if ( b1[j]>b2[j] ) b[j] = b2[j];
         if ( e1[j]>=e2[j] ) e[j] = e1[j];
         if ( e1[j]<e2[j] ) e[j] = e2[j];
         npt[j] = (int) ((e[j]-b[j])/dt+2);
         e[j] = b[j]+(npt[j]-1)*dt;
      }
	   tp = hd1.t1; ts = hd1.t2;
  	   syn[j]=(float *) calloc(npt[j], sizeof(float));
	   if (dura>0.&& j==0) src = trapezoid(dura,rise,dt,&ns);
	}
	else if (hd1.npts != npt1) {
	   fprintf(stderr,"number points in %s not agree with %d\n",nam1,npt1);
	}
   else if (core) {
	   if (hd2.npts != npt2) {
         fprintf(stderr,"number points in %s not agree with %d\n",nam2,npt2);
      }
   }
        for(pt=syn[j],k=0;k<npt1;k++,pt++) (*pt) += coef*grn1[k];
        free(grn1);
        (*ccc)++;
        if (core) {
           pt_shift = (int) (fabs(b1[j]-b2[j])/dt);
           for(pt=syn[j]+pt_shift,k=0;k<npt2;k++,pt++) (*pt) += coef*grn2[k];
           free(grn2);
           (*ddd)++;
        }
        
	if (*ccc == '9') {
      (*ccc)+=40;	/* explosion components start at 'a' instead of '0' */
      if (core) (*ddd)+=40;
   }
     }
  }

  if (!dynamic) {
     printf("%8.2f %8.2f %10.3e %10.3e %10.3e\n",dist,az,disp[0],disp[1],disp[2]);
     return 0;
  }

  /* convolve a source time function. integrate or filtering if needed */
#ifdef SAC_LIB
  if (filter) design(order, type, proto, 1., 1., f1, f2, (double) dt, sn, sd, &nsects);
#endif
  if (tstar>0.) fttq_(&dt,&tstar,&mftm,&nftm,ftm);
  for(j=0;j<3;j++) {
     if (ns > 0) conv(src,ns,syn[j],npt[j]);
     if (intg) cumsum(syn[j],npt[j],dt);
     if (diff) diffrt(syn[j],npt[j],dt);
#ifdef SAC_LIB
     if (filter) apply(syn[j],(long int)npt,0,sn,sd,nsects);
#endif
     if (tstar>0.) conv(ftm,nftm,syn[j],npt[j]);
  }

  /* output */
  ccc = outnm + strlen(outnm) - 1;
  if (ns == 0) strcat(outnm,"i");
//   hd.b -= shift;
//   hd.e -= shift;
  hd.t1 = tp; hd.t2 = ts;
  hd.az = az;
  for(j=0;j<3;j++) {
     *ccc = com[j];
     hd.b = b[j]-shift;
     hd.e = e[j]-shift;
     hd.npts = npt[j];
     hd.cmpinc = cmpinc[j];
     hd.cmpaz  = cmpaz[j];
     write_sac(outnm,hd,syn[j]);
  }

  return 0;

}

/* construct a trapezoid of duration dura with rise portion rise */
/* and an area of 1*dt */
float *trapezoid(float dura, float rise, float dt, int *ns) {
   int i, nr;
   float amp, *src;
   *ns = rint(dura/dt); if (*ns<2) *ns = 2;
   src = malloc((1+*ns)*sizeof(float));
   if (src == NULL) {*ns=0; return src;}
   nr = rint(rise*(*ns)); if (nr<1) nr = 1; if (2*nr>*ns) nr=*ns/2;
   amp = 1./(nr*(*ns-nr));
   for(i=0;i<nr;i++) src[i] = i*amp;
   for(;i<*ns-nr;i++) src[i] = nr*amp;
   for(;i<=*ns;i++)   src[i] = (*ns-i)*amp;
   (*ns)++;
   return src;
}
