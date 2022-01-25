#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "MD_trial.h"
#include"water.h"
#include "rand.h"
#include "wanglandau.h"

extern int myrank;
void sphere(double ox, double oy, double oz, double* npos){
	double v1, v2,v3,a;

	v1 = (2.0*randd1()-1.0);
	v2 = (2.0*randd1()-1.0);
	a = v1*v1 +v2*v2;
	while(a>1.0){
		v1 = (2.0*randd1()-1.0);
		v2 = (2.0*randd1()-1.0);	
		a = v1*v1+v2*v2;
	}
	v3 = 1-2.0*a;
	a = sqrt(1-a);
	v1 = 2.0*a*v1;
	v2 = 2.0*a*v2;

	npos[0] = ox + Roh*v1;
	npos[1] = oy + Roh*v2;
	npos[2] = oz + Roh*v3;

}
void cone(double ox, double oy, double oz, double h1x, double h1y, double h1z, double* npos){
   double w1,w2,w3,v1,v2,v3,t1,t2,t3,a,bb,c,d,h;

	w1 = ox - h1x;
	w2 = oy - h1y;
	w3 = oz - h1z;
	
	a = w1*w1+w2*w2+w3*w3;
	a=sqrt(a);
	w1 = w1/a;
	w2 = w2/a;
	w3 = w3/a;
    
	v1 = (2.0*randd1()-1.0);
	v2 = (2.0*randd1()-1.0);
	a = v1*v1 +v2*v2;
	while(a>1.0){
		v1 = (2.0*randd1()-1.0);
		v2 = (2.0*randd1()-1.0);	
		a = v1*v1+v2*v2;
	}

	a = w1*w1+w2*w2-1.0;
	bb = -2.0*(v1*w1*w3+v2*w2*w3);
	c =v2*v2*w1*w1+v1*v1*w2*w2+v1*v1*w3*w3+v2*v2*w3*w3-2.0*v1*v2*w1*w2-v1*v1-v2*v2;
	d = bb*bb-4.0*a*c;

	if(d<0.0 && d>-1.0e-7)
		v3=-bb/(2.0*a);
	else if(randd1()<0.5)
		v3=(-bb+sqrt(d))/(2.0*a);
	else
		v3=(-bb-sqrt(d))/(2.0*a);

	a=sqrt(v1*v1+v2*v2+v3*v3);
	v1=v1/a;
	v2=v2/a;
	v3=v3/a;

	h=-1.0/tan(Angle*Pi/180.0);
	t1=h*w1+v1;
	t2=h*w2+v2;
	t3=h*w3+v3;
	a=sqrt(t1*t1+t2*t2+t3*t3);
	t1/=a;
	t2/=a;
	t3/=a;
	
	npos[0]=ox+Roh*t1;
	npos[1]=oy+Roh*t2;
	npos[2]=oz+Roh*t3;

}



int constrain(int n,  int flag){
	double Msum,xm,ym,zm,r;
	int i, inside;
	if(flag==0) // rotation about cm; 
		return 1;

	r = Lc0*Lc0;

	inside = 1;
	for(i = 0 ; i < N; i++){
		if(i == n) continue;

		if((Wmol[n].cx-Wmol[i].cx)*(Wmol[n].cx-Wmol[i].cx) >r ||(Wmol[n].cy-Wmol[i].cy)*(Wmol[n].cy-Wmol[i].cy) > r || (Wmol[n].cz-Wmol[i].cz)*(Wmol[n].cz-Wmol[i].cz) >r){
			inside = 0;
			break;
		}	
	}

	return inside;
}

int constrain_ini(const struct Water* Wclus, int n){
    double xm,ym,zm,r,r1;
	double Msum;
	double* ixm = malloc(sizeof(double)*N);
	double* iym = malloc(sizeof(double)*N);
	double* izm = malloc(sizeof(double)*N);
	int i;

	xm=0.0;
    ym=0.0;
	zm=0.0;

	for(i=0;i<N;i++){
		ixm[i] = Wclus[i].Ox *Mo +(Wclus[i].H1x+Wclus[i].H2x)*Mh;
		iym[i] = Wclus[i].Oy *Mo +(Wclus[i].H1y+Wclus[i].H2y)*Mh;
		izm[i] = Wclus[i].Oz *Mo +(Wclus[i].H1z+Wclus[i].H2z)*Mh;
	}
	Msum =( Mo+2*Mh);
	for(i=0;i<N;i++){
		ixm[i]/=Msum;
		iym[i]/=Msum;
		izm[i]/=Msum;

		xm +=ixm[i];
		ym +=iym[i];
		zm +=izm[i];
	}

	xm/=N;
	ym/=N;
	zm/=N;

	r = Rconstrain*Rconstrain;
	r1 = (ixm[n]-xm)*(ixm[n]-xm)+(iym[n]-ym)*(iym[n]-ym)+(izm[n]-zm)*(izm[n]-zm);

//	r1 = (Wclus[n].Ox-xm)*(Wclus[n].Ox-xm)+(Wclus[n].Oy-ym)*(Wclus[n].Oy-ym)+(Wclus[n].Oz-zm)*(Wclus[n].Oz-zm);
//	r2 = (Wclus[n].H1x-xm)*(Wclus[n].H1x-xm)+(Wclus[n].H1y-ym)*(Wclus[n].H1y-ym)+(Wclus[n].H1z-zm)*(Wclus[n].H1z-zm);
//  r3=  (Wclus[n].H2x-xm)*(Wclus[n].H2x-xm)+(Wclus[n].H2y-ym)*(Wclus[n].H2y-ym)+(Wclus[n].H2z-zm)*(Wclus[n].H2z-zm);
	free(ixm);
	free(iym);
	free(izm);

	if(r1>r)
		return 0;
	else 
		return 1;
}

void ini_water(int Read_conf, const char* filename){
		double rHH, rhh;	
	if(Read_conf==0)
	{
		int i,j,k,cnt,L;
		double r = Rconstrain;
		double x,y,z;
		double H1pos[3],H2pos[3];

		int flag;
		Wmol = malloc( N * sizeof(struct Water) );
// lattice ini. with constrain on Roo;
		L = 0;
		i = 0;
		while(i < N){
			++L;
			i = L*L*L;
		}
		cnt = -1;
		for(i=0; i< L; i++)
			for(j=0; j<L; j++)
				for(k=0; k<L; k++){
					++cnt;

					if(cnt < N){

						Wmol[cnt].Ox = x = i*Rconstrain;
						Wmol[cnt].Oy = y = j*Rconstrain;
						Wmol[cnt].Oz = z = k*Rconstrain;

						sphere(x,y,z,H1pos);
						Wmol[cnt].H1x = H1pos[0];
						Wmol[cnt].H1y = H1pos[1];
						Wmol[cnt].H1z = H1pos[2];
						assert( fabs((x-H1pos[0])*(x-H1pos[0]) + (y-H1pos[1])*(y-H1pos[1]) +(z-H1pos[2])*(z-H1pos[2]) - Roh*Roh) < 1e-7);

						cone(x,y,z,H1pos[0],H1pos[1],H1pos[2], H2pos);
						Wmol[cnt].H2x = H2pos[0];
						Wmol[cnt].H2y = H2pos[1];
						Wmol[cnt].H2z = H2pos[2];
						assert( fabs((x-H2pos[0])*(x-H2pos[0]) + (y-H2pos[1])*(y-H2pos[1]) +(z-H2pos[2])*(z-H2pos[2]) - Roh*Roh) < 1e-7);
					}
				}
		Lc0 = 1.1*Rconstrain*L; // make sure all the mol. inside the box;
	/*	do{
			for(i=0;i<N;i++)
			{
				do{	
				x = r* (2.0*randd1()-1.0);
				y = r* (2.0*randd1()-1.0);
				z = r* (2.0*randd1()-1.0);
				}while(x*x+y*y+z*z > r*r);
				Wmol[i].Ox = x;
				Wmol[i].Oy = y;
				Wmol[i].Oz = z;
				sphere(x,y,z,H1pos);
				Wmol[i].H1x = H1pos[0];
				Wmol[i].H1y = H1pos[1];
				Wmol[i].H1z = H1pos[2];
				assert( fabs((x-H1pos[0])*(x-H1pos[0]) + (y-H1pos[1])*(y-H1pos[1]) +(z-H1pos[2])*(z-H1pos[2]) - Roh*Roh) < 1e-7);
				cone(x,y,z,H1pos[0],H1pos[1],H1pos[2], H2pos);
				Wmol[i].H2x = H2pos[0];
				Wmol[i].H2y = H2pos[1];
				Wmol[i].H2z = H2pos[2];
				assert( fabs((x-H2pos[0])*(x-H2pos[0]) + (y-H2pos[1])*(y-H2pos[1]) +(z-H2pos[2])*(z-H2pos[2]) - Roh*Roh) < 1e-7);
#ifndef NDEBUG
				rHH =Roh*Roh*2 - 2*Roh*Roh*cos(Angle*Pi/180.0);
				rhh = (Wmol[i].H1x-Wmol[i].H2x)*(Wmol[i].H1x-Wmol[i].H2x) + (Wmol[i].H1y-Wmol[i].H2y)*(Wmol[i].H1y-Wmol[i].H2y) +(Wmol[i].H1z-Wmol[i].H2z)*(Wmol[i].H1z-Wmol[i].H2z);
				assert(fabs(rHH-rhh) < 1e-7);
#endif
			}
			for(i=0;i<N;i++)
				if(constrain_ini(Wmol,i)==0){
					flag=0;
					break;
				}
				else
					flag=1;
		}while(Etot(Wmol)>0.0 || flag==0); */
	}
	else
	{
	//  read conf. from file;
		FILE * ifp;
		int i;
		char dum[8];
		double x,y,z;
		Wmol = malloc( N * sizeof(struct Water) );
		
		ifp=fopen( filename,"r");
		if(ifp==NULL){
			printf("The file was not opened\n");
			exit(1);
		}
		else
		{
			for(i=0;i<N;++i)
			{
				fscanf(ifp,"%*s %lf  %lf %lf \n", &x, &y, &z );
				Wmol[i].Ox=x;
				Wmol[i].Oy=y;
				Wmol[i].Oz=z;
				fscanf(ifp,"%*s %lf  %lf %lf \n", &x, &y, &z );
				Wmol[i].H1x=x;
				Wmol[i].H1y=y;
				Wmol[i].H1z=z;
				fscanf(ifp,"%*s %lf  %lf %lf \n", &x, &y, &z );
				Wmol[i].H2x=x;
				Wmol[i].H2y=y;
				Wmol[i].H2z=z;
#ifndef NDEBUG
                rhh= (Wmol[i].Ox-Wmol[i].H2x)*(Wmol[i].Ox-Wmol[i].H2x) + (Wmol[i].Oy-Wmol[i].H2y)*(Wmol[i].Oy-Wmol[i].H2y) +(Wmol[i].Oz-Wmol[i].H2z)*(Wmol[i].Oz-Wmol[i].H2z); 
				assert( fabs(rhh- Roh*Roh) < 1e-5);
				rHH =Roh*Roh*2 - 2*Roh*Roh*cos(Angle*Pi/180.0);
				rhh = (Wmol[i].H1x-Wmol[i].H2x)*(Wmol[i].H1x-Wmol[i].H2x) + (Wmol[i].H1y-Wmol[i].H2y)*(Wmol[i].H1y-Wmol[i].H2y) +(Wmol[i].H1z-Wmol[i].H2z)*(Wmol[i].H1z-Wmol[i].H2z);
				assert(fabs(rHH-rhh) < 1e-5);
#endif
			};
			fclose(ifp);

//			for(i=0;i<N;i++)
//			constrain(Wmol,i);
		}
	}
}
double Eupdate(int i, const struct Water* iWmol){
	double deltaE=0.0;
	int j;
	double roo6, r12;
	double temp,t,Ei;
	double xm[2],ym[2],zm[2],xc[2],yc[2],zc[2];
	double Msum;
	double temp1,temp2;
	for(j=0;j<N;j++){
		if(j!=i){
			temp = (iWmol[i].Ox- iWmol[j].Ox)*(iWmol[i].Ox- iWmol[j].Ox)+ (iWmol[i].Oy- iWmol[j].Oy)* (iWmol[i].Oy- iWmol[j].Oy)+ (iWmol[i].Oz- iWmol[j].Oz)* (iWmol[i].Oz- iWmol[j].Oz);
			roo6 = temp*temp*temp;
			temp1= (Alj/roo6 -Clj)/roo6;

#ifdef SPCE
			r12 = sqrt(temp);
			temp2= Qo*Qo/r12;
		    r12 =  (iWmol[i].H1x- iWmol[j].Ox)*(iWmol[i].H1x- iWmol[j].Ox)+ (iWmol[i].H1y- iWmol[j].Oy)* (iWmol[i].H1y- iWmol[j].Oy)+ (iWmol[i].H1z- iWmol[j].Oz)* (iWmol[i].H1z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;
		    r12 =  (iWmol[i].H2x- iWmol[j].Ox)*(iWmol[i].H2x- iWmol[j].Ox)+ (iWmol[i].H2y- iWmol[j].Oy)* (iWmol[i].H2y- iWmol[j].Oy)+ (iWmol[i].H2z- iWmol[j].Oz)* (iWmol[i].H2z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H1x)*(iWmol[i].Ox- iWmol[j].H1x)+ (iWmol[i].Oy- iWmol[j].H1y)* (iWmol[i].Oy- iWmol[j].H1y)+ (iWmol[i].Oz- iWmol[j].H1z)* (iWmol[i].Oz- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H2x)*(iWmol[i].Ox- iWmol[j].H2x)+ (iWmol[i].Oy- iWmol[j].H2y)* (iWmol[i].Oy- iWmol[j].H2y)+ (iWmol[i].Oz- iWmol[j].H2z)* (iWmol[i].Oz- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
#endif

#ifdef TIP3P
			r12 = sqrt(temp);
			temp2= Qo*Qo/r12;
		    r12 =  (iWmol[i].H1x- iWmol[j].Ox)*(iWmol[i].H1x- iWmol[j].Ox)+ (iWmol[i].H1y- iWmol[j].Oy)* (iWmol[i].H1y- iWmol[j].Oy)+ (iWmol[i].H1z- iWmol[j].Oz)* (iWmol[i].H1z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;
		    r12 =  (iWmol[i].H2x- iWmol[j].Ox)*(iWmol[i].H2x- iWmol[j].Ox)+ (iWmol[i].H2y- iWmol[j].Oy)* (iWmol[i].H2y- iWmol[j].Oy)+ (iWmol[i].H2z- iWmol[j].Oz)* (iWmol[i].H2z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H1x)*(iWmol[i].Ox- iWmol[j].H1x)+ (iWmol[i].Oy- iWmol[j].H1y)* (iWmol[i].Oy- iWmol[j].H1y)+ (iWmol[i].Oz- iWmol[j].H1z)* (iWmol[i].Oz- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H2x)*(iWmol[i].Ox- iWmol[j].H2x)+ (iWmol[i].Oy- iWmol[j].H2y)* (iWmol[i].Oy- iWmol[j].H2y)+ (iWmol[i].Oz- iWmol[j].H2z)* (iWmol[i].Oz- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
#endif

#ifdef TIP4P
		  xc[0] = iWmol[i].Ox * Mo + (iWmol[i].H1x+iWmol[i].H2x)*Mh;
	 	  yc[0] = iWmol[i].Oy * Mo + (iWmol[i].H1y+iWmol[i].H2y)*Mh;
	 	  zc[0] = iWmol[i].Oz * Mo + (iWmol[i].H1z+iWmol[i].H2z)*Mh;

		  xc[1] = iWmol[j].Ox * Mo + (iWmol[j].H1x+iWmol[j].H2x)*Mh;
	 	  yc[1] = iWmol[j].Oy * Mo + (iWmol[j].H1y+iWmol[j].H2y)*Mh;
	 	  zc[1] = iWmol[j].Oz * Mo + (iWmol[j].H1z+iWmol[j].H2z)*Mh;
		  
		   Msum = (Mo+2*Mh);
		   xc[0] = xc[0]/Msum;
		   yc[0] = yc[0]/Msum;
		   zc[0] = zc[0]/Msum;		
			
	//       t = sqrt(  (iWmol[i].Ox - xc[0])*(iWmol[i].Ox - xc[0]) + (iWmol[i].Oy - yc[0])*(iWmol[i].Oy - yc[0]) + (iWmol[i].Oz - zc[0])*(iWmol[i].Oz - zc[0]) );
		   t = Rom/Roc;
		   xm[0] = iWmol[i].Ox - (iWmol[i].Ox - xc[0])*t;  // line equation in 3D;
		   ym[0] = iWmol[i].Oy - (iWmol[i].Oy - yc[0])*t;
		   zm[0] = iWmol[i].Oz - (iWmol[i].Oz - zc[0])*t;

		   xc[1] = xc[1]/Msum;
		   yc[1] = yc[1]/Msum;
		   zc[1] = zc[1]/Msum;	

	 //      t = sqrt(  (iWmol[j].Ox - xc[1])*(iWmol[j].Ox - xc[1]) + (iWmol[j].Oy - yc[1])*(iWmol[j].Oy - yc[1]) + (iWmol[j].Oz - zc[1])*(iWmol[j].Oz - zc[1]) );
	//	   t = Rom/t;
		   xm[1] = iWmol[j].Ox - (iWmol[j].Ox - xc[1])*t;  // line equation in 3D;
		   ym[1] = iWmol[j].Oy - (iWmol[j].Oy - yc[1])*t;
		   zm[1] = iWmol[j].Oz - (iWmol[j].Oz - zc[1])*t;

	//	   r12 = (xm[0]- iWmol[i].Ox)*(xm[0]- iWmol[i].Ox)+ (ym[0]- iWmol[i].Oy)* (ym[0]- iWmol[i].Oy)+ (zm[0]- iWmol[i].Oz)* (zm[0]- iWmol[i].Oz);

			r12 =  (xm[0]- xm[1])* (xm[0]- xm[1])+(ym[0]- ym[1])* (ym[0]- ym[1])+(zm[0]- zm[1])* (zm[0]- zm[1]);
			r12 = sqrt(r12);
			temp2= Qm*Qm/r12;
		    r12 =  (iWmol[i].H1x- xm[1])*(iWmol[i].H1x- xm[1])+ (iWmol[i].H1y- ym[1])* (iWmol[i].H1y- ym[1])+ (iWmol[i].H1z- zm[1])* (iWmol[i].H1z- zm[1]);
			r12 =sqrt(r12);
			temp2+= Qh*Qm/r12;
		    r12 =  (iWmol[i].H2x- xm[1])*(iWmol[i].H2x- xm[1])+ (iWmol[i].H2y- ym[1])* (iWmol[i].H2y- ym[1])+ (iWmol[i].H2z- zm[1])* (iWmol[i].H2z- zm[1]);
			r12 =sqrt(r12);
			temp2+= Qh*Qm/r12;

		    r12 =  (xm[0]- iWmol[j].H1x)*(xm[0]- iWmol[j].H1x)+ (ym[0]- iWmol[j].H1y)* (ym[0]- iWmol[j].H1y)+ (zm[0]- iWmol[j].H1z)* (zm[0]- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qm*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (xm[0]- iWmol[j].H2x)*(xm[0]- iWmol[j].H2x)+ (ym[0]- iWmol[j].H2y)* (ym[0]- iWmol[j].H2y)+ (zm[0]- iWmol[j].H2z)* (zm[0]- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qm*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

#endif
			Ei=Ep[i][j];
			Ep[i][j]=Ep[j][i]= temp1+temp2*ke2;
			deltaE += Ep[i][j] - Ei;
		}
	}
	return deltaE;
} 
double Etot(const struct Water* iWmol){
	double Elj=0.0;
	double Ec=0.0;
	int i,j;
	double roo6, r12;
	double temp,t;
	double xm[2],ym[2],zm[2],xc[2],yc[2],zc[2];
	double Msum;
	double temp1, temp2;

	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			Ep[i][j] = 0.0;

	for(i=0 ; i< N-1; i++){
		for(j=i+1; j< N;j++){
			temp = (iWmol[i].Ox- iWmol[j].Ox)*(iWmol[i].Ox- iWmol[j].Ox)+ (iWmol[i].Oy- iWmol[j].Oy)* (iWmol[i].Oy- iWmol[j].Oy)+ (iWmol[i].Oz- iWmol[j].Oz)* (iWmol[i].Oz- iWmol[j].Oz);
			roo6 = temp*temp*temp;
			temp1= (Alj/roo6 -Clj)/roo6;

#ifdef SPCE
			r12 = sqrt(temp);
			temp2= Qo*Qo/r12;
		    r12 =  (iWmol[i].H1x- iWmol[j].Ox)*(iWmol[i].H1x- iWmol[j].Ox)+ (iWmol[i].H1y- iWmol[j].Oy)* (iWmol[i].H1y- iWmol[j].Oy)+ (iWmol[i].H1z- iWmol[j].Oz)* (iWmol[i].H1z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;
		    r12 =  (iWmol[i].H2x- iWmol[j].Ox)*(iWmol[i].H2x- iWmol[j].Ox)+ (iWmol[i].H2y- iWmol[j].Oy)* (iWmol[i].H2y- iWmol[j].Oy)+ (iWmol[i].H2z- iWmol[j].Oz)* (iWmol[i].H2z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H1x)*(iWmol[i].Ox- iWmol[j].H1x)+ (iWmol[i].Oy- iWmol[j].H1y)* (iWmol[i].Oy- iWmol[j].H1y)+ (iWmol[i].Oz- iWmol[j].H1z)* (iWmol[i].Oz- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H2x)*(iWmol[i].Ox- iWmol[j].H2x)+ (iWmol[i].Oy- iWmol[j].H2y)* (iWmol[i].Oy- iWmol[j].H2y)+ (iWmol[i].Oz- iWmol[j].H2z)* (iWmol[i].Oz- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
#endif

#ifdef TIP3P
			r12 = sqrt(temp);
			temp2= Qo*Qo/r12;
		    r12 =  (iWmol[i].H1x- iWmol[j].Ox)*(iWmol[i].H1x- iWmol[j].Ox)+ (iWmol[i].H1y- iWmol[j].Oy)* (iWmol[i].H1y- iWmol[j].Oy)+ (iWmol[i].H1z- iWmol[j].Oz)* (iWmol[i].H1z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;
		    r12 =  (iWmol[i].H2x- iWmol[j].Ox)*(iWmol[i].H2x- iWmol[j].Ox)+ (iWmol[i].H2y- iWmol[j].Oy)* (iWmol[i].H2y- iWmol[j].Oy)+ (iWmol[i].H2z- iWmol[j].Oz)* (iWmol[i].H2z- iWmol[j].Oz);
			r12 =sqrt(r12);
			temp2+= Qh*Qo/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H1x)*(iWmol[i].Ox- iWmol[j].H1x)+ (iWmol[i].Oy- iWmol[j].H1y)* (iWmol[i].Oy- iWmol[j].H1y)+ (iWmol[i].Oz- iWmol[j].H1z)* (iWmol[i].Oz- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (iWmol[i].Ox- iWmol[j].H2x)*(iWmol[i].Ox- iWmol[j].H2x)+ (iWmol[i].Oy- iWmol[j].H2y)* (iWmol[i].Oy- iWmol[j].H2y)+ (iWmol[i].Oz- iWmol[j].H2z)* (iWmol[i].Oz- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qo*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
#endif

#ifdef TIP4P
		  xc[0] = iWmol[i].Ox * Mo + (iWmol[i].H1x+iWmol[i].H2x)*Mh;
	 	  yc[0] = iWmol[i].Oy * Mo + (iWmol[i].H1y+iWmol[i].H2y)*Mh;
	 	  zc[0] = iWmol[i].Oz * Mo + (iWmol[i].H1z+iWmol[i].H2z)*Mh;

		  xc[1] = iWmol[j].Ox * Mo + (iWmol[j].H1x+iWmol[j].H2x)*Mh;
	 	  yc[1] = iWmol[j].Oy * Mo + (iWmol[j].H1y+iWmol[j].H2y)*Mh;
	 	  zc[1] = iWmol[j].Oz * Mo + (iWmol[j].H1z+iWmol[j].H2z)*Mh;
		  
		   Msum = (Mo+2*Mh);
		   xc[0] = xc[0]/Msum;
		   yc[0] = yc[0]/Msum;
		   zc[0] = zc[0]/Msum;		
			
	//       t = sqrt(  (iWmol[i].Ox - xc[0])*(iWmol[i].Ox - xc[0]) + (iWmol[i].Oy - yc[0])*(iWmol[i].Oy - yc[0]) + (iWmol[i].Oz - zc[0])*(iWmol[i].Oz - zc[0]) );
		   t = Rom/Roc;		//t;
		   xm[0] = iWmol[i].Ox - (iWmol[i].Ox - xc[0])*t;  // line equation in 3D;
		   ym[0] = iWmol[i].Oy - (iWmol[i].Oy - yc[0])*t;
		   zm[0] = iWmol[i].Oz - (iWmol[i].Oz - zc[0])*t;

		   xc[1] = xc[1]/Msum;
		   yc[1] = yc[1]/Msum;
		   zc[1] = zc[1]/Msum;	

//	       t = sqrt(  (iWmol[j].Ox - xc[1])*(iWmol[j].Ox - xc[1]) + (iWmol[j].Oy - yc[1])*(iWmol[j].Oy - yc[1]) + (iWmol[j].Oz - zc[1])*(iWmol[j].Oz - zc[1]) );
//		   t = Rom/t;
		   xm[1] = iWmol[j].Ox - (iWmol[j].Ox - xc[1])*t;  // line equation in 3D;
		   ym[1] = iWmol[j].Oy - (iWmol[j].Oy - yc[1])*t;
		   zm[1] = iWmol[j].Oz - (iWmol[j].Oz - zc[1])*t;


	//	   r12 = (xm[0]- iWmol[i].Ox)*(xm[0]- iWmol[i].Ox)+ (ym[0]- iWmol[i].Oy)* (ym[0]- iWmol[i].Oy)+ (zm[0]- iWmol[i].Oz)* (zm[0]- iWmol[i].Oz);

			r12 =  (xm[0]- xm[1])* (xm[0]- xm[1])+(ym[0]- ym[1])* (ym[0]- ym[1])+(zm[0]- zm[1])* (zm[0]- zm[1]);
			r12 = sqrt(r12);
			temp2= Qm*Qm/r12;
		    r12 =  (iWmol[i].H1x- xm[1])*(iWmol[i].H1x- xm[1])+ (iWmol[i].H1y- ym[1])* (iWmol[i].H1y- ym[1])+ (iWmol[i].H1z- zm[1])* (iWmol[i].H1z- zm[1]);
			r12 =sqrt(r12);
			temp2+= Qh*Qm/r12;
		    r12 =  (iWmol[i].H2x- xm[1])*(iWmol[i].H2x- xm[1])+ (iWmol[i].H2y- ym[1])* (iWmol[i].H2y- ym[1])+ (iWmol[i].H2z- zm[1])* (iWmol[i].H2z- zm[1]);
			r12 =sqrt(r12);
			temp2+= Qh*Qm/r12;

		    r12 =  (xm[0]- iWmol[j].H1x)*(xm[0]- iWmol[j].H1x)+ (ym[0]- iWmol[j].H1y)* (ym[0]- iWmol[j].H1y)+ (zm[0]- iWmol[j].H1z)* (zm[0]- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qm*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H1x)*(iWmol[i].H1x- iWmol[j].H1x)+ (iWmol[i].H1y- iWmol[j].H1y)* (iWmol[i].H1y- iWmol[j].H1y)+ (iWmol[i].H1z- iWmol[j].H1z)* (iWmol[i].H1z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H1x)*(iWmol[i].H2x- iWmol[j].H1x)+ (iWmol[i].H2y- iWmol[j].H1y)* (iWmol[i].H2y- iWmol[j].H1y)+ (iWmol[i].H2z- iWmol[j].H1z)* (iWmol[i].H2z- iWmol[j].H1z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

		    r12 =  (xm[0]- iWmol[j].H2x)*(xm[0]- iWmol[j].H2x)+ (ym[0]- iWmol[j].H2y)* (ym[0]- iWmol[j].H2y)+ (zm[0]- iWmol[j].H2z)* (zm[0]- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qm*Qh/r12;
			 r12 =  (iWmol[i].H1x- iWmol[j].H2x)*(iWmol[i].H1x- iWmol[j].H2x)+ (iWmol[i].H1y- iWmol[j].H2y)* (iWmol[i].H1y- iWmol[j].H2y)+ (iWmol[i].H1z- iWmol[j].H2z)* (iWmol[i].H1z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;
			 r12 =  (iWmol[i].H2x- iWmol[j].H2x)*(iWmol[i].H2x- iWmol[j].H2x)+ (iWmol[i].H2y- iWmol[j].H2y)* (iWmol[i].H2y- iWmol[j].H2y)+ (iWmol[i].H2z- iWmol[j].H2z)* (iWmol[i].H2z- iWmol[j].H2z);
			r12 =sqrt(r12);
			temp2+= Qh*Qh/r12;

#endif
			Elj += temp1;
			Ec += temp2;
			Ep[i][j]=Ep[j][i]= temp1+temp2*ke2;
	//		Epold[i][j]=Epold[j][i]=Ep[i][j];
		}
	}
	Ec *= ke2;
	return Elj+Ec;
}
void Vol(){
	double Lnew, RatL, RevL, Ei, Ef, Ev,  rHH, rhh;
	double AXX,AXY,AXZ,AYX,AYY,AYZ,AZX,AZY,AZZ;
	int i, j, iti;

	for(i=0; i<N; i++)
		Wold[i] = Wmol[i];


	Ei = currEtot;
	 iti = (int) ((Ei*invN-WLD1min)*invdWLD1);
	Lnew = Lc0 + (2.0*randd1() - 1.0)*Lmax;
	++attemptV[iti];

	if(Lnew > 0){
		RatL = Lc0/Lnew;
		RevL = 1.0/RatL;
		ini_quat();
		for(i = 0 ;i < N; i++){

			Wmol[i].cx  *= RevL;
			Wmol[i].cy *= RevL;
			Wmol[i].cz  *= RevL;

				AXX = quat[i].q0 *quat[i].q0 +quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			   AXY = 2.0 * (quat[i].q1*quat[i].q2+ quat[i].q0 *quat[i].q3);
			   AXZ = 2.0 * (quat[i].q1*quat[i].q3- quat[i].q0 *quat[i].q2);
				AYX = 2.0 * (quat[i].q1*quat[i].q2- quat[i].q0 *quat[i].q3);
				AYY = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 +quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
				AYZ = 2.0 * (quat[i].q2*quat[i].q3+ quat[i].q0 *quat[i].q1);
				AZX = 2.0 * (quat[i].q1*quat[i].q3+ quat[i].q0 *quat[i].q2);
				AZY = 2.0 * (quat[i].q2*quat[i].q3- quat[i].q0 *quat[i].q1);
				AZZ = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 +quat[i].q3*quat[i].q3;

			   Wmol[i].Ox = (Wmol[i].cx -0.0207858*Ul * AZX) ;
			   Wmol[i].Oy = (Wmol[i].cy-0.0207858*Ul  *  AZY) ;
			   Wmol[i].Oz = (Wmol[i].cz-0.0207858*Ul  *  AZZ) ;

			   Wmol[i].H1x = (Wmol[i].cx + 0.239997*Ul  * AYX +0.1649726*Ul   * AZX) ;
			   Wmol[i].H1y = (Wmol[i].cy +0.239997*Ul  * AYY +0.1649726 *Ul  * AZY) ;
			   Wmol[i].H1z = (Wmol[i].cz +0.239997*Ul  * AYZ +0.1649726 *Ul  * AZZ) ;

			   Wmol[i].H2x = (Wmol[i].cx-0.239997*Ul  * AYX +0.1649726*Ul   * AZX) ;
			   Wmol[i].H2y = (Wmol[i].cy-0.239997*Ul  * AYY +0.1649726 *Ul  * AZY) ;
			   Wmol[i].H2z = (Wmol[i].cz-0.239997*Ul  * AYZ +0.1649726*Ul   * AZZ) ;	


		}

		Ef = Etot(Wmol);

		Ev = 3*N*log(RatL) + Pres/136.26*(Lnew*Lnew*Lnew - Lc0*Lc0*Lc0);

				if(WangLandau(Ei,Ef, Ev,0,0,0,-2 )==1)
				{
					++acceptV[iti];
					currEtot = Ef;	
					Lc0 = Lnew;
					Lc[iti] = Lnew;
					
					for(i=0;i<N;i++)
						for(j=0;j<N;j++)
							Epold[i][j] = Epold[j][i]  = Ep[i][j];

				}
				else
				{
					//reject - return monomer to old position
					Lc[iti] = Lc0;
					for(j=0; j<N; j++)
						Wmol[j] = Wold[j];
					currEtot = Ei;
					for(i=0;i<N;i++)
						for(j=0;j<N;j++)
							Ep[i][j] = Ep[j][i]  = Epold[i][j];

				};		

	}
}



void Diff(){
	int i=0;
	struct Water WmolOld; 
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
    int j, axis,iti,flag;
	double XMsum,YMsum,ZMsum,Msum,xc,yc,zc,temp1,temp2,tx,ty,tz;
	double alpha,beta,gamma,ca,cb,cr,sa,sb,sr;
	double cosD, sinD,dr,df;
	double RM[4][4];
	double x,y,z;
	double rand;
    double c,s,u,rx,ry,rz,eta1,eta2,etasq;
	double Dfi, Dri;
	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
	 Ei = currEtot;
	
	 iti = (int) ((Ei*invN-WLD1min)*invdWLD1);
	Dfi = Df[iti];
	//Dri = Dr[iti];

 //diffuse the ith monomer
	WmolOld=Wmol[i];
	tx = ty = tz =0.0;
	df=dr=0.0;
//	if( randd1() < ProbD){// Diffusion;
			flag = 1;

     		tx = Dfi*(2.0*randd1()-1.0);
			ty = Dfi*(2.0*randd1()-1.0);
			tz = Dfi*(2.0*randd1()-1.0);
			

			Wmol[i].Ox += tx;
			Wmol[i].Oy += ty;
			Wmol[i].Oz += tz;

			Wmol[i].H1x += tx;
			Wmol[i].H1y += ty;
			Wmol[i].H1z += tz;

			Wmol[i].H2x += tx;
			Wmol[i].H2y += ty;
			Wmol[i].H2z += tz;

			Wmol[i].cx += tx;
			Wmol[i].cy += ty;
			Wmol[i].cz += tz;
/*
		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
			rx = 2.0*eta1*sqrt(1.0-etasq);
			ry = 2.0*eta2*sqrt(1.0-etasq);
			rz = 1.0 - 2.0*etasq;

			df= Dfi*randd1();
            tx = df*rx;
			ty = df*ry;
			tz = df*rz;


			Wmol[i].Ox += tx;
			Wmol[i].Oy += ty;
			Wmol[i].Oz += tz;

			Wmol[i].H1x += tx;
			Wmol[i].H1y += ty;
			Wmol[i].H1z += tz;

			Wmol[i].H2x += tx;
			Wmol[i].H2y += ty;
			Wmol[i].H2z += tz;
*/
		//	attemptdiff +=1;


	//	}

		if(constrain(i, flag)==1){
			 Ef=Eupdate(i,Wmol)+Ei;
		  
			if(WangLandau(Ei,Ef, dr,tx,ty,tz,flag )==1)
			{
				currEtot = Ef;	

					for(j=0;j<N;j++){
						if(j!=i){
							Epold[i][j]=Epold[j][i]=Ep[i][j];
						}
					}

			}
			else
			{
				//reject - return monomer to old position
				Wmol[i]=WmolOld;
				currEtot = Ei;

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

			};		


		}else{
					//reject - return monomer to old position
               WangLandau(Ei,Ei,dr,tx,ty,tz, flag  ); // update old state;
				Wmol[i]=WmolOld;
				currEtot = Ei;	

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

		}
}

void Rot(){
	int i=0;
	struct Water WmolOld; 
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
    int j, axis,iti,flag;
	double XMsum,YMsum,ZMsum,Msum,xc,yc,zc,temp1,temp2,tx,ty,tz;
	double alpha,beta,gamma,ca,cb,cr,sa,sb,sr;
	double cosD, sinD,dr,df;
	double RM[4][4];
	double x,y,z;
	double rand;
    double c,s,u,rx,ry,rz,eta1,eta2,etasq;
	double Dfi, Dri;
	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
	 Ei = currEtot;
	
	 iti = (int) ((Ei*invN-WLD1min)*invdWLD1);
//	Dfi = Df[iti];
	Dri = Dr[iti];

 //diffuse the ith monomer
	WmolOld=Wmol[i];
	tx = ty = tz =0.0;
	df=dr=0.0;

//		else{	// Rotation trial about center of mass;
			flag =0;
			//XMsum=YMsum=ZMsum=0.0;
			xc = Wmol[i].cx;
			yc = Wmol[i].cy;
			zc = Wmol[i].cz;


/*	       dr=(2.0*randd1()-1.0)*Dri;
		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
			rx = 2.0*eta1*sqrt(1.0-etasq);
			ry = 2.0*eta2*sqrt(1.0-etasq);
			rz = 1.0 - 2.0*etasq;

			//angle constants
			c=cos(dr);
			s=sin(dr);
			u=1.0-c;		  
	
			x = Wmol[i].Ox-xc;
			y = Wmol[i].Oy-yc;
			z = Wmol[i].Oz-zc;

			//Perform the rotation on the selected segment
			Wmol[i].Ox = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].Oy = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].Oz = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = Wmol[i].H1x-xc;
			y = Wmol[i].H1y-yc;
			z = Wmol[i].H1z-zc;
			Wmol[i].H1x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].H1y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].H1z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = Wmol[i].H2x-xc;
			y = Wmol[i].H2y-yc;
			z = Wmol[i].H2z-zc;
			Wmol[i].H2x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].H2y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].H2z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

*/

			dr=(2.0*randd1()-1.0)*Dri;
			cosD = cos(dr);
			sinD = sin(dr);
			axis =  (int)(3.0*randd1() )+1;
			switch(axis){
				case(1): // y-z plane;
					ty = Wmol[i].Oy - yc;
					tz = Wmol[i].Oz - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].Oy = temp1 + yc;
					Wmol[i].Oz = temp2 + zc;

					ty = Wmol[i].H1y - yc;
					tz = Wmol[i].H1z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].H1y = temp1 + yc;
					Wmol[i].H1z = temp2 + zc;

					ty = Wmol[i].H2y - yc;
					tz = Wmol[i].H2z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].H2y = temp1 + yc;
					Wmol[i].H2z = temp2 + zc;
					break;
				case(2): // x-z plane
					tz = Wmol[i].Oz - zc;
					tx = Wmol[i].Ox - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].Oz = temp1 + zc;
					Wmol[i].Ox = temp2 + xc;

					tz = Wmol[i].H1z - zc;
					tx = Wmol[i].H1x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].H1z = temp1 + zc;
					Wmol[i].H1x = temp2 + xc;

					tz = Wmol[i].H2z - zc;
					tx = Wmol[i].H2x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].H2z = temp1 + zc;
					Wmol[i].H2x = temp2 + xc;
					break;
				case(3): // x-y plane;
					tx = Wmol[i].Ox - xc;
					ty = Wmol[i].Oy- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].Ox = temp1 + xc;
					Wmol[i].Oy = temp2 + yc;

					tx = Wmol[i].H1x - xc;
					ty = Wmol[i].H1y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].H1x = temp1 + xc;
					Wmol[i].H1y = temp2 + yc;

					tx = Wmol[i].H2x - xc;
					ty = Wmol[i].H2y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].H2x = temp1 + xc;
					Wmol[i].H2y = temp2 + yc;
					break;
			};  
		
	//	attemptrot +=1;
			
//		};
//		}while(constrain(Wmol, i)!=1);
    //Calculate the final energy

//		fti=(int) ((Ef*invN-WLD1min)*invdWLD1);

		if(constrain(i, flag)==1){
			 Ef=Eupdate(i,Wmol)+Ei;
		  
			if(WangLandau(Ei,Ef, dr,tx,ty,tz,flag )==1)
			{
				currEtot = Ef;	

					for(j=0;j<N;j++){
						if(j!=i){
							Epold[i][j]=Epold[j][i]=Ep[i][j];
						}
					}

			}
			else
			{
				//reject - return monomer to old position
				Wmol[i]=WmolOld;
				currEtot = Ei;

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

			};		


		}else{
					//reject - return monomer to old position
               WangLandau(Ei,Ei,dr,tx,ty,tz, flag  ); // update old state;
				Wmol[i]=WmolOld;
				currEtot = Ei;	

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

		}
}



void DiffandRot(){
	int i=0;
	struct Water WmolOld; 
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
    int j, axis,flag,iti;
	double XMsum,YMsum,ZMsum,Msum,xc,yc,zc,temp1,temp2,tx,ty,tz;
	double alpha,beta,gamma,ca,cb,cr,sa,sb,sr;
	double cosD, sinD,dr,df;
	double RM[4][4];
	double x,y,z;
    double c,s,u,rx,ry,rz,eta1,eta2,etasq;
	double Dfi, Dri;
	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
//	Ei=Etot();
//	assert(abs(currEtot-Etot(Wmol) ) < 1e-7);
	 Ei =currEtot;

	 iti = (int) ((Ei*invN-WLD1min)*invdWLD1);
	Dfi = Df[iti];
	Dri = Dr[iti];

	//diffuse the ith monomer
	WmolOld=Wmol[i];
	tx=ty=tz=dr=0.0;

//	do{
		if( randd1() < ProbD){// Diffusion;
			flag = 1;

			tx = Dfi*(2.0*randd1()-1.0);
			ty = Dfi*(2.0*randd1()-1.0);
			tz = Dfi*(2.0*randd1()-1.0);
			

			Wmol[i].Ox += tx;
			Wmol[i].Oy += ty;
			Wmol[i].Oz += tz;

			Wmol[i].H1x += tx;
			Wmol[i].H1y += ty;
			Wmol[i].H1z += tz;

			Wmol[i].H2x += tx;
			Wmol[i].H2y += ty;
			Wmol[i].H2z += tz;

			Wmol[i].cx += tx;
			Wmol[i].cy += ty;
			Wmol[i].cz += tz;
/*
		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
			rx = 2.0*eta1*sqrt(1.0-etasq);
			ry = 2.0*eta2*sqrt(1.0-etasq);
			rz = 1.0 - 2.0*etasq;
 
			df= Dfi*randd1();
            tx = df*rx;
			ty = df*ry;
			tz = df*rz;

			Wmol[i].Ox += tx;
			Wmol[i].Oy += ty;
			Wmol[i].Oz += tz;

			Wmol[i].H1x += tx;
			Wmol[i].H1y += ty;
			Wmol[i].H1z += tz;

			Wmol[i].H2x += tx;
			Wmol[i].H2y += ty;
			Wmol[i].H2z += tz;

*/


		}
		else{	// Rotation trial about center of mass;
			flag=0;
			//XMsum=YMsum=ZMsum=0.0;

			xc = Wmol[i].cx;
			yc = Wmol[i].cy;
			zc = Wmol[i].cz;


/*	       dr=(2.0*randd1()-1.0)*Pi*Dri;
		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
			rx = 2.0*eta1*sqrt(1.0-etasq);
			ry = 2.0*eta2*sqrt(1.0-etasq);
			rz = 1.0 - 2.0*etasq;

			//angle constants
			c=cos(dr);
			s=sin(dr);
			u=1.0-c;		  
	
			x = Wmol[i].Ox-xc;
			y = Wmol[i].Oy-yc;
			z = Wmol[i].Oz-zc;

			//Perform the rotation on the selected segment
			Wmol[i].Ox = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].Oy = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].Oz = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = Wmol[i].H1x-xc;
			y = Wmol[i].H1y-yc;
			z = Wmol[i].H1z-zc;
			Wmol[i].H1x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].H1y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].H1z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = Wmol[i].H2x-xc;
			y = Wmol[i].H2y-yc;
			z = Wmol[i].H2z-zc;
			Wmol[i].H2x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			Wmol[i].H2y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			Wmol[i].H2z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

*/

			dr=(2.0*randd1()-1.0)*Dri;
			cosD = cos(dr);
			sinD = sin(dr);
			axis =  (int)(3.0*randd1() )+1;
			switch(axis){
				case(1): // y-z plane;
					ty = Wmol[i].Oy - yc;
					tz = Wmol[i].Oz - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].Oy = temp1 + yc;
					Wmol[i].Oz = temp2 + zc;

					ty = Wmol[i].H1y - yc;
					tz = Wmol[i].H1z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].H1y = temp1 + yc;
					Wmol[i].H1z = temp2 + zc;

					ty = Wmol[i].H2y - yc;
					tz = Wmol[i].H2z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					Wmol[i].H2y = temp1 + yc;
					Wmol[i].H2z = temp2 + zc;
					break;
				case(2): // x-z plane
					tz = Wmol[i].Oz - zc;
					tx = Wmol[i].Ox - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].Oz = temp1 + zc;
					Wmol[i].Ox = temp2 + xc;

					tz = Wmol[i].H1z - zc;
					tx = Wmol[i].H1x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].H1z = temp1 + zc;
					Wmol[i].H1x = temp2 + xc;

					tz = Wmol[i].H2z - zc;
					tx = Wmol[i].H2x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					Wmol[i].H2z = temp1 + zc;
					Wmol[i].H2x = temp2 + xc;
					break;
				case(3): // x-y plane;
					tx = Wmol[i].Ox - xc;
					ty = Wmol[i].Oy- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].Ox = temp1 + xc;
					Wmol[i].Oy = temp2 + yc;

					tx = Wmol[i].H1x - xc;
					ty = Wmol[i].H1y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].H1x = temp1 + xc;
					Wmol[i].H1y = temp2 + yc;

					tx = Wmol[i].H2x - xc;
					ty = Wmol[i].H2y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					Wmol[i].H2x = temp1 + xc;
					Wmol[i].H2y = temp2 + yc;
					break;
			};  
			

			
		};
//		}while(constrain(Wmol, i)!=1);
    //Calculate the final energy
		if(constrain(i, flag)==1){
			 Ef=Eupdate(i,Wmol)+Ei;
		  
			if(WangLandau(Ei,Ef,dr,tx, ty, tz, flag)==1)
			{
				currEtot = Ef;	

					for(j=0;j<N;j++){
						if(j!=i){
							Epold[i][j]=Epold[j][i]=Ep[i][j];
						}
					}

			}
			else
			{
				//reject - return monomer to old position
				Wmol[i]=WmolOld;
				currEtot = Ei;

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

			};		
		}else{
					//reject - return monomer to old position
               WangLandau(Ei,Ei,dr,tx,ty,tz,flag); // update old state;
				Wmol[i]=WmolOld;
				currEtot = Ei;	

					for(j=0;j<N;j++){
						if(j!=i){
							Ep[i][j]=Ep[j][i]=Epold[i][j];
						}
					}

		}
}


void reset_XYZ(){
	int i,j;
	double cx=0.0,cy=0.0,cz=0.0;
    double eng, Msum, xm, ym, zm,r ;
	//Find the approx. center of mass for each component
	for(i=0;i<N;i++)
    {
      cx+=Wmol[i].Ox;
      cy+=Wmol[i].Oy;
      cz+=Wmol[i].Oz;
    };

	cx*=invN;
	cy*=invN;
	cz*=invN;

	//fprintf(stderr,"%g\t%g\t%g\n",cmass(),currEtot,currnonbondE);
	for(i=0;i<N;i++)
    {
		Wmol[i].Ox=Wmol[i].Ox-cx;
		Wmol[i].Oy=Wmol[i].Oy-cy;
		Wmol[i].Oz=Wmol[i].Oz-cz;
		
		Wmol[i].H1x=Wmol[i].H1x-cx;
		Wmol[i].H1y=Wmol[i].H1y-cy;
		Wmol[i].H1z=Wmol[i].H1z-cz;

		Wmol[i].H2x=Wmol[i].H2x-cx;
		Wmol[i].H2y=Wmol[i].H2y-cy;
		Wmol[i].H2z=Wmol[i].H2z-cz;

    };
	

	for(i=0;i<N;i++){
		Wmol[i].cx = (Wmol[i].Ox *Mo +(Wmol[i].H1x+Wmol[i].H2x)*Mh ) /(2*Mh+Mo) ;
		Wmol[i].cy = ( Wmol[i].Oy *Mo +(Wmol[i].H1y+Wmol[i].H2y)*Mh) /(2*Mh+Mo) ;
		Wmol[i].cz = (Wmol[i].Oz *Mo +(Wmol[i].H1z+Wmol[i].H2z)*Mh) /(2*Mh+Mo) ;
	}

	 r = Lc0*Lc0;
	for(i = 0; i< N-1; i++)
		for(j = i+1; j<N;j++)
			if((Wmol[j].cx-Wmol[i].cx)*(Wmol[j].cx-Wmol[i].cx) >r ||(Wmol[j].cy-Wmol[i].cy)*(Wmol[j].cy-Wmol[i].cy) > r || (Wmol[j].cz-Wmol[i].cz)*(Wmol[j].cz-Wmol[i].cz) >r)
				exit(119);
	


	//Recalculate the overall energies
	eng=Etot(Wmol);
	//printf("myrank=%d, delta=%g\n", myrank, eng-currEtot);
	assert(fabs(eng-currEtot) < 1e-5);

}

void count_Hbond(){
	int i,j;
	for(i=0; i<N-1 ;i++)
		for(j=i+1; j<N; j++){
			
		}
}


void orderqs(void){

}
void grcalculate(int pbin){

}
double coredensity(void){
	return 0.0;
}
double Rgyr2(void){
	return 0.0;
}
double EEdist(void){
	return 0.0;
}
void initialize(const char* filename){
	int i,j;

   	//for relatively big jumps
	D=0.2*Pi;  DD=0.1;
	D0=0.02*Pi;DD0=0.01;

	Ep=(double**)malloc(N*sizeof(double));
	for(i=0;i<N;i++)
		Ep[i]=(double*)malloc(N*sizeof(double));

	Epold=(double**)malloc(N*sizeof(double));
	for(i=0;i<N;i++)
		Epold[i]=(double*)malloc(N*sizeof(double));


	//precalculations for chain length
	invN=1.0/(1.0*N);
	rootinvN = 1.0/sqrt(1.0*N);
	if(filename == NULL)
		ini_water(0, NULL);
	else
		ini_water(1, filename);

	Wold = malloc(sizeof(struct Water)*N);

	for(i=0;i<N;i++){
		Wmol[i].cx = (Wmol[i].Ox *Mo +(Wmol[i].H1x+Wmol[i].H2x)*Mh ) /(2*Mh+Mo) ;
		Wmol[i].cy = ( Wmol[i].Oy *Mo +(Wmol[i].H1y+Wmol[i].H2y)*Mh) /(2*Mh+Mo) ;
		Wmol[i].cz = (Wmol[i].Oz *Mo +(Wmol[i].H1z+Wmol[i].H2z)*Mh) /(2*Mh+Mo) ;
	}

	currEtot=Etot(Wmol);

	
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			Epold[i][j] = Epold[j][i]  = Ep[i][j];
}

void writeinc(int frame, double energy)
{
  int i,j;
  char s[512];
  FILE *ofp;
  double cx=0.0,cy=0.0,cz=0.0;

  for(i=0;i<N;++i)
    {
      cx+=Wmol[i].Ox;
      cy+=Wmol[i].Oy;
      cz+=Wmol[i].Oz;
    };

  cx*=invN;
  cy*=invN;
  cz*=invN;

  sprintf(s,"snap_%d.inc",frame);
  ofp=fopen(s,"w");
  
  fprintf(ofp,"//N=%d, Etot=%18.10e\n",N,currEtot*invN);
  fprintf(ofp,"#declare polychain=union{\n");

  for(i=0;i<N-1;++i)
    {
      fprintf(ofp,"sphere{<%g,%g,%g>, Rs texture{stex }}\n", (Wmol[i].Ox), (Wmol[i].Oy), (Wmol[i].Oz));
      fprintf(ofp,"cylinder{<%g,%g,%g>,<%g,%g,%g>, Rc texture{ctex }}\n",
	      (Wmol[i].Ox), (Wmol[i].Oy), (Wmol[i].Oz),  (Wmol[i+1].Ox), (Wmol[i+1].Oy), (Wmol[i+1].Oz) );
    };
  fprintf(ofp,"sphere{<%g,%g,%g>, Rs texture{stex}}\n", (Wmol[N-1].Ox), (Wmol[N-1].Oy), (Wmol[N-1].Oz));
  
  fprintf(ofp,"finish{ambient 0.25 diffuse 0.75 specular 0.3}\ntranslate -<%g,%g,%g>\n}\n",cx,cy,cz);
  //fprintf(ofp,"finish{ambient 0.25 diffuse 0.75 specular 0.3}\n}\n");
  fflush(ofp);
  fclose(ofp);
}



//similar to writeinc, but writes out a list of the x,y,z monomer positions
void write_mol2(int frame, const struct Water* iWmol)
{
  int i,j;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"snap_%d.mol2",frame);
  ofp=fopen(s,"w");

  fprintf(ofp,"@<TRIPOS>MOLECULE\n");
  fprintf(ofp,"WATER CLUSTER\n");
  fprintf(ofp," %d %d\n",N*3,N*2);
  fprintf(ofp,"SMALL\n");
  fprintf(ofp,"NO_CHARGES\n");
  fprintf(ofp,"@<TRIPOS>ATOM\n");

  for(i=0;i<N*3;i+=3)
    {
      fprintf(ofp,"%d MON%d %g %g %g O.spc\n",i+1,i+1, (iWmol[i/3].Ox), (iWmol[i/3].Oy), (iWmol[i/3].Oz) );
      fprintf(ofp,"%d MON%d %g %g %g H.spc\n",i+2,i+2, (iWmol[i/3].H1x), (iWmol[i/3].H1y), (iWmol[i/3].H1z) );
      fprintf(ofp,"%d MON%d %g %g %g H.spc\n",i+3,i+3, (iWmol[i/3].H2x), (iWmol[i/3].H2y), (iWmol[i/3].H2z) );
		//	&(Wmol[i].H1x),&(Wmol[i].H1y),&(Wmol[i].H1z),&(Wmol[i].H2x),&(Wmol[i].H2y),&(Wmol[i].H2z));
    
    };
  fprintf(ofp,"@<TRIPOS>BOND\n");

  j=0;
  for(i=0;i<N;++i)
    {
	  j++;
      fprintf(ofp,"%d %d %d %d\n",j,i*3+1,i*3+2,1);
	  j++;
      fprintf(ofp,"%d %d %d %d\n",j,i*3+1,i*3+3,1);    
    };

  fclose(ofp);

}

//Calculates thermodynamic quantities and writes modification-factor labeled files
void thermoqs()
{
	int i;	
	double T,U,Z,C,F,S;
	double centerE,lambda,Bw,Nfree;
	
	FILE *therm_op;
	char s1[512];
	
	//Opening File containing all thermodynamic quantities (in following order)
	//	T	U	Cv	freeE  Entropy	Rgyr2  EEdist
	sprintf(s1,"therm.dat");
	therm_op=fopen(s1,"w");
	
	if( (therm_op==NULL) )
	{
		fprintf(stderr, "\nHey, this file ( in thermoqs() ) could not be opened!\n\n");
		exit(1);
	};
	
	//Normalization for free energy
	Nfree = 0.0;
	
	//Main Temperature Loop
	for(T=TTi;T<=TTf+dTT;T=T+dTT)
	{			
		U = 0.0;	// Initializing the average energy <E>	
		Z = 0.0;	// Initializing the Partition Function 		
		C = 0.0;	// Initialize the Specific Heat
		F = 0.0;	// Initialize the free energy
		S = 0.0;	// Initialize the entropy
		Bw = 0.0;	//Boltzmann weight
		lambda = -1.0e300;	//Normalization shift (max value of DOS considering T)
		centerE = 0.0;	// Taking the center of the energy bin
		
		//Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
		for(i=0;i<D1BINS;i++)
		{
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
				
			if( (lambda < ((wllng[i]) - 1.0*(centerE)/(Rg*T))) )// && (wllng[i] > 0.0) ) 
				lambda = ((wllng[i]) - 1.0*(centerE)/(Rg*T));  
			
		};
		
		//Central Loop for calculating thermodynamic properties from the DOS
		for(i=0;i<D1BINS;i++)
		{			
			//Taking the center of the bin
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
			
			//Boltzmann Factor
			Bw =  exp( wllng[i]  - (centerE)/(Rg*T) - lambda );
		
			//Partition Function
			Z = Z + Bw;
				
			//Average Energy <E>
			U = U + (centerE)*Bw;		   		   
				
			//Average Energy Squared <E^2>
			C = C + (centerE)*(centerE)*Bw;								
		};
		
		//Internal Energy
		U = U/Z;	
		//Normalization of free energy
		if(T == TTi)
		{
			Nfree = -(Rg*T)*( lambda + log(Z) ) - U;
		};
		//Specific Heat
		C = ( (C/Z) - (U*U) ) / (Rg*T*T);	
		//Free Energy
		F = -(Rg*T)*( lambda + log(Z) ) - Nfree;
		//Entropy
		S = (U - F)/T;
			
		fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\n",T,U*invN,C*invN+3*Rg,F*invN,S*invN);	
	};	
		
	fclose(therm_op);
}

