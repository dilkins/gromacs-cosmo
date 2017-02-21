// From Numerical Recipes in C.
#include <math.h>
#include "nrutil.h"
void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;

	float *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}


void polin2(float x1a[], float x2a[], float **ya, int m, int n, float x1,
float x2, float *y, float *dy)
{
	int j;
	float *ymtmp;
	ymtmp=vector(1,m);
	for (j=1;j<=m;j++) {
		polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	polint(x1a,ymtmp,m,x1,y,dy);
	free_vector(ymtmp,1,m);
}

void polin3(float x1a[], float x2a[], float x3a[], float ***yb, int l, int m, int n, float x1, float x2, float x3, float *y, float *dy)
{
	int j,k;
	float **ya;
//	float **ymtmp;

//??	for (j=1;j<=l;j++)
	for (k=0;k<n;k++)
	{
		// First take cuts in the z direction. For this value of l, we create a square. This square
		// will have x1 and x2 coordinates, and will have an array **ya, which we can pass into polin2.
		for (i=0;i<l;i++)
		{
			for (j=0;j<m;j++)
			{
				ya[i][j] = yb[i][j][k];
			}
		}
		polin2(x1a,x2a,ya,l,m,x1,x2,&ymtmp[j],dy);
	}



	for (j=1;j<=l;j++)
	{
		ymtmp = vector(1,m);
		for (k=1;k<=m;k++)
		{
			polint(x3a,ya[j][k],n,x3,&ymtmp[k],dy);
		}
	}
	
}
