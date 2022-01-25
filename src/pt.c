#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "pt.h"
#include "water.h"
#include "rand.h"

void ini_T(double Ti, double Tf, int nT){
	int i;
	if(Ti*Tf > 1e-5){
		pT = malloc(sizeof(double)*nT);
		if(pT != NULL){
			pT[0] = log(Ti);
			pT[nT-1]=log(Tf);
			for(i=1;i<nT-1;i++){
				pT[i] = pT[0] +(pT[nT-1]-pT[0])*i/(nT-1);
			}
			for(i=0;i<nT;i++)
				pT[i]=exp(pT[i]);
		}
		else
			exit(1);
	}else{
		exit(1);
	}
}

void ini_sys(int nT,const char* filename){
       int i;
	   attemps=malloc(sizeof(double)*nT);
	   accepts=malloc(sizeof(double)*nT);
	   pE= malloc(sizeof(double)*nT);
        pW= malloc(sizeof(struct Water*)*nT);
		if(pW!=NULL){
			for(i=0;i<nT;i++){
				attemps[i]=accepts[i]=0.0;
				initialize(filename);
				pW[i] = Wmol;
				pE[i] = Etot(Wmol);
			}
			attemps[nT-1]=-1;
  
		}else
			exit(1);
}

void swap(int nT){
	int i;
	double delta, etemp;
	struct Water* wtemp;

	for(i=0;i<nT-1;i++){
		attemps[i] +=1;
		delta = (pE[i]-pE[i+1])*(1/(Rg*pT[i+1])-1/(Rg*pT[i]));
		if(delta <=0 || randd1() < exp(-delta)){
		   wtemp = pW[i+1];
		   pW[i+1] = pW[i];
		   pW[i] = wtemp;

		   etemp = pE[i+1];
		   pE[i+1] = pE[i];
		   pE[i] = etemp;

		   accepts[i] +=1;
		}
	}

}
int Metropolis(double Ei, double Ef, int n)
{
  double R;

  if(Ef<Ei)
    {
      //accept
//      currEtot=Ef;
  //    currnonbondE=Efastf;
      return 1;
    }
  else
    {
	  if(n==-1)
		  R=exp(-(Ef-Ei)/(Rg*T));
	  else
		R=exp(-(Ef-Ei)/(Rg*pT[n]));

      if(randd1()<R)
	{
	  //accept
//	  currEtot=Ef;
//	  currnonbondE=Efastf;
	  return 1;
	}
      else
	{
	  //reject
	  return 0;
	};
    };
}

void mcdiff(int n){
//	assert(n>=0 && n<sizeof(pT)/sizeof(double));
	struct Water* iWmol = pW[n];

	int i=0;
	struct Water WmolOld; 
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
    int j, axis,flag;
	double XMsum,YMsum,ZMsum,Msum,xc,yc,zc,temp1,temp2,tx,ty,tz;
	double cosD,sinD,dtheta;
	double alpha,beta,gamma,ca,cb,cr,sa,sb,sr;
	double RM[4][4];
	double x,y,z;
	double rand;
   double rHH, rhh;
    double c,s,u,rx,ry,rz,eta1,eta2,etasq;
	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
//	Ei=Etot(iWmol);
    Ei=pE[n];

	assert(fabs(Ei-Etot(iWmol)) < 1e-7);
	//diffuse the ith monomer
	WmolOld=iWmol[i];

//	do{
		if(randd1() < ProbD){// Diffusion;
			flag = 1;
/*			tx = DD*(2.0*randd1()-1.0);
			ty = DD*(2.0*randd1()-1.0);
			tz = DD*(2.0*randd1()-1.0);
			

			iWmol[i].Ox += tx;
			iWmol[i].Oy += ty;
			iWmol[i].Oz += tz;

			iWmol[i].H1x += tx;
			iWmol[i].H1y += ty;
			iWmol[i].H1z += tz;

			iWmol[i].H2x += tx;
			iWmol[i].H2y += ty;
			iWmol[i].H2z += tz;
*/
		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
			rx = 2.0*eta1*sqrt(1.0-etasq);
			ry = 2.0*eta2*sqrt(1.0-etasq);
			rz = 1.0 - 2.0*etasq;

			temp1= DD*randd1();
            tx = temp1*rx;
			ty = temp1*ry;
			tz = temp1*rz;

			iWmol[i].Ox += tx;
			iWmol[i].Oy += ty;
			iWmol[i].Oz += tz;

			iWmol[i].H1x += tx;
			iWmol[i].H1y += ty;
			iWmol[i].H1z += tz;

			iWmol[i].H2x += tx;
			iWmol[i].H2y += ty;
			iWmol[i].H2z += tz;


		}
		else{	// Rotation trial about center of mass;
			flag = 0;
		//	XMsum=YMsum=ZMsum=0.0;

			XMsum = iWmol[i].Ox * Mo + (iWmol[i].H1x+iWmol[i].H2x)*Mh;
			YMsum = iWmol[i].Oy * Mo + (iWmol[i].H1y+iWmol[i].H2y)*Mh;
			ZMsum = iWmol[i].Oz * Mo + (iWmol[i].H1z+iWmol[i].H2z)*Mh;

			Msum = (Mo+2*Mh);
			xc = XMsum/Msum;
			yc = YMsum/Msum;
			zc = ZMsum/Msum;

	       dtheta=(2.0*randd1()-1.0)*Pi*D;
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
			c=cos(dtheta);
			s=sin(dtheta);
			u=1.0-c;		  
	
			//Perform the rotation on the selected segment
			x = iWmol[i].Ox-xc;
			y = iWmol[i].Oy-yc;
			z = iWmol[i].Oz-zc;
			iWmol[i].Ox = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			iWmol[i].Oy = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			iWmol[i].Oz = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = iWmol[i].H1x-xc;
			y = iWmol[i].H1y-yc;
			z = iWmol[i].H1z-zc;
			iWmol[i].H1x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			iWmol[i].H1y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			iWmol[i].H1z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	

			x = iWmol[i].H2x-xc;
			y = iWmol[i].H2y-yc;
			z = iWmol[i].H2z-zc;
			iWmol[i].H2x = (u*rx*rx + c)*x     + (u*ry*rx - s*rz)*y  + (u*rz*rx + ry*s)*z +xc;
			iWmol[i].H2y = (u*rx*ry + rz*s)*x  + (u*ry*ry + c)*y     + (u*rz*ry - rx*s)*z +yc;
			iWmol[i].H2z = (u*rx*rz - ry*s)*x  + (u*ry*rz + rx*s)*y  + (u*rz*rz + c)*z +zc;	
  /*   		alpha=(2.0*randd1()-1.0)*Pi*D;
			beta=(2.0*randd1()-1.0)*Pi*D;
			gamma=(2.0*randd1()-1.0)*Pi*D;

	//		temp1 = alpha*alpha;
	//		temp2 = temp1*temp1;
			ca=  cos(alpha);   // 1-temp1/2.0 + temp2/24.0;
			sa =  sin(alpha); //alpha*(1 - temp1/6.0 + temp2/120.0);

	//		temp1 = beta*beta;
	//		temp2 = temp1*temp1;
			cb= cos(beta);  // 1-temp1/2.0 + temp2/24.0;
			sb = sin(beta);   //beta*(1 - temp1/6.0 + temp2/120.0);

	//		temp1 = gamma*gamma;
	//		temp2 = temp1*temp1;
			cr= cos(gamma); //1-temp1/2.0 + temp2/24.0;
			sr = sin(gamma); //gamma*(1 - temp1/6.0 + temp2/120.0);

			RM[1][1] = ca*cb;
			RM[1][2] = -sr*sb*ca-sa*cr;
			RM[1][3] = -cr*sb*ca+sa*sr;

			RM[2][1] = -sa*cb;
			RM[2][2] = sa*sr*sb-ca*cr;
			RM[2][3] = sa*cr*sb+ca*sr;

			RM[3][1] = sb;
			RM[3][2] = cb*sr;
			RM[3][3] = cb*cr;


			x = iWmol[i].Ox-xc;
			y = iWmol[i].Oy-yc;
			z = iWmol[i].Oz-zc;
			iWmol[i].Ox = x*RM[1][1] + y*RM[2][1] + z*RM[3][1]+xc;
			iWmol[i].Oy = x*RM[1][2] + y*RM[2][2] + z*RM[3][2]+yc ;
			iWmol[i].Oz = x*RM[1][3] + y*RM[2][3] + z*RM[3][3]+zc ;

			x = iWmol[i].H1x-xc;
			y = iWmol[i].H1y-yc;
			z = iWmol[i].H1z-zc;
			iWmol[i].H1x = x*RM[1][1] + y*RM[2][1] + z*RM[3][1]+xc;
			iWmol[i].H1y = x*RM[1][2] + y*RM[2][2] + z*RM[3][2]+yc ;
			iWmol[i].H1z = x*RM[1][3] + y*RM[2][3] + z*RM[3][3]+zc ;

			x = iWmol[i].H2x-xc;
			y = iWmol[i].H2y-yc;
			z = iWmol[i].H2z-zc;
			iWmol[i].H2x = x*RM[1][1] + y*RM[2][1] + z*RM[3][1]+xc;
			iWmol[i].H2y = x*RM[1][2] + y*RM[2][2] + z*RM[3][2]+yc ;
			iWmol[i].H2z = x*RM[1][3] + y*RM[2][3] + z*RM[3][3]+zc ;
			*/
/*

			dtheta=(2.0*randd1()-1.0)*Pi*D;
//			temp1 = dtheta*dtheta;
//			temp2 = temp1*temp1;
			cosD = cos(dtheta);  //1-temp1/2.0 + temp2/24.0;
			sinD =  sin(dtheta); //dtheta*(1 - temp1/6.0 + temp2/120.0);

			axis =  (int)(3.0*randd1() )+1;
			switch(axis){
				case(1): // y-z plane;
					ty = iWmol[i].Oy - yc;
					tz = iWmol[i].Oz - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					iWmol[i].Oy = temp1 + yc;
					iWmol[i].Oz = temp2 + zc;

					ty = iWmol[i].H1y - yc;
					tz = iWmol[i].H1z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					iWmol[i].H1y = temp1 + yc;
					iWmol[i].H1z = temp2 + zc;

					ty = iWmol[i].H2y - yc;
					tz = iWmol[i].H2z - zc;
					temp1= cosD*ty-sinD*tz;
					temp2 = sinD*ty+cosD*tz;
					iWmol[i].H2y = temp1 + yc;
					iWmol[i].H2z = temp2 + zc;
					break;
				case(2): // x-z plane
					tz = iWmol[i].Oz - zc;
					tx = iWmol[i].Ox - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					iWmol[i].Oz = temp1 + zc;
					iWmol[i].Ox = temp2 + xc;

					tz = iWmol[i].H1z - zc;
					tx = iWmol[i].H1x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					iWmol[i].H1z = temp1 + zc;
					iWmol[i].H1x = temp2 + xc;

					tz = iWmol[i].H2z - zc;
					tx = iWmol[i].H2x - xc;
					temp1= cosD*tz-sinD*tx;
					temp2 = sinD*tz+cosD*tx;
					iWmol[i].H2z = temp1 + zc;
					iWmol[i].H2x = temp2 + xc;
					break;
				case(3): // x-y plane;
					tx = iWmol[i].Ox - xc;
					ty = iWmol[i].Oy- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					iWmol[i].Ox = temp1 + xc;
					iWmol[i].Oy = temp2 + yc;

					tx = iWmol[i].H1x - xc;
					ty = iWmol[i].H1y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					iWmol[i].H1x = temp1 + xc;
					iWmol[i].H1y = temp2 + yc;

					tx = iWmol[i].H2x - xc;
					ty = iWmol[i].H2y- yc;
					temp1= cosD*tx-sinD*ty;
					temp2 = sinD*tx+cosD*ty;
					iWmol[i].H2x = temp1 + xc;
					iWmol[i].H2y = temp2 + yc;
					break;
			}; */

			 attemptrot+=1;
		}; 

#ifndef NDEBUG
				rHH =Roh*Roh*2 - 2*Roh*Roh*cos(Angle*Pi/180.0);
				rhh = (iWmol[i].H1x-iWmol[i].H2x)*(iWmol[i].H1x-iWmol[i].H2x) + (iWmol[i].H1y-iWmol[i].H2y)*(iWmol[i].H1y-iWmol[i].H2y) +(iWmol[i].H1z-iWmol[i].H2z)*(iWmol[i].H1z-iWmol[i].H2z);
				assert(fabs(rHH-rhh) < 1e-5);
#endif


//		}while(constrain(iWmol, i) ==0);
		if(constrain(i,  flag)==1){
			//Calculate the final energy
			Ef=Etot(iWmol);
		  
			if(Metropolis(Ei,Ef,n)==1)
			{
				pE[n] = Ef;

			}
			else
			{
				//reject - return monomer to old position
				iWmol[i]=WmolOld;
				pE[n] = Ei;
			};		
		}else{
					//reject - return monomer to old position
				iWmol[i]=WmolOld;
				pE[n] = Ei;	
		}

}

void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP){
	int mcs,i,j,k,n;
	 FILE *ofp;
	 FILE * ofp2;

	double * avgE=malloc(sizeof(double)*nT);
	double * avgE2=malloc(sizeof(double)*nT);
	double * c=malloc(sizeof(double)*nT);

	double* Et = malloc(sizeof(double)*SAMPS);
	double GS = 0.0;
    ofp = fopen("stat.dat","w");
	for(n=0;n<nT;n++){
		avgE[n]=0.0;
		avgE2[n]=0.0;
	}


	for(mcs=0;mcs<DROPI;mcs++){
		for(n=0;n<nT;n++)
			for(i =0; i< N; i++)
				mcdiff(n);
		if(mcs%10==0)
			swap(nT);
	}

	  for(i=0;i<SAMPS;++i)
		{
			for(j=0;j<SEP;++j){
				for(n=0;n<nT;n++)
					for(k=0;k<N;k++)
						mcdiff(n);
				if(i%10==0)
					swap(nT);
			}
//			Et[i] = pE[0];
			for(n=0;n<nT;n++){							
				avgE[n] +=pE[n];
				avgE2[n] +=pE[n]*pE[n];
			}
			if(pE[0] < GS){
				write_mol2(0,pW[0]);
				GS=pE[0];	
			}
			
		};
	for(n=0;n<nT;n++){
	  avgE[n]/=1.0*SAMPS;
      avgE2[n]/=1.0*SAMPS;
	  c[n]= (avgE2[n] - avgE[n]*avgE[n])/(Rg*pT[n]*pT[n]) ;
	}

	for(n=0;n<nT;n++){
      fprintf(ofp,"%g\t%g\t%g\n",pT[n],avgE[n]*invN,c[n]*invN+3*Rg);
	}	
//	for(i=0;i<SAMPS;i++)
//			fprintf(ofp,"%d\t%g\n",i,Et[i]);
	  fflush(ofp);

	 ofp2 = fopen("misc.dat","w");
	 fprintf(ofp2,"Ground state energy:%g\n",GS);
	// fprintf(ofp2,"Avg_diff:%g  Avg_rot:%g Sample size:%d\n",1.0*acceptdiff/attemptdiff, 1.0*acceptrot/attemptrot,  (int)SAMPS);
	 fprintf(ofp2,"swap prob:\n");
	 for(n=0;n<nT-1;n++){
	 fprintf(ofp2,"%g\t%g\n",pT[n], 1.0*accepts[n]/attemps[n]);		
	 }
}