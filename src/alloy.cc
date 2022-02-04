#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define  DEFINE_GLOBALS
#include "alloy.hpp"
#include "pt.hpp"
#include "rand.hpp"
#include "wanglandau.hpp"

extern int myrank, nprocs;

void initialize(int state){

	//maximum rotational angle;
	D=1.0*Pi; DD = 0.2;

	attd=accd=0;

	//initialize spin configuration;
	ini_conf(state);
	ini_coupling();
	ini_alloy();
//	write_mol2(-1);
	invN=1.0/(N*N*N);
	

}

void ini_coupling(){
	int i,j,k, shell, ti, tj;
	double tmp, Jr;
	FILE* fop;
	for(i=0; i<5; i++)
		for(j=0; j<5; j++)
			for(k=0; k<15; k++)
				J[i][j][k] = 1.0;

	NS[0] = 12; NS[1] = 6;
	NS[2] = 24; NS[3] = 12;
	NS[4] = 24; NS[5] = 8;
	NS[6] = 48; NS[7] = 6;
	NS[8] = 36; NS[9] = 24;
// read from input file;
	fop = fopen("source_code/data/J.dat","r");
   if(fop==NULL){
			printf("Jfit.dat file was not opened\n");
			exit(1);
   }
	for(i=0;i<150; i++) // 10 shell, each with 15 J;
    {
		fscanf(fop,"%d %lg %d %d %lg\n",&shell, &tmp, &ti, &tj, &Jr);
		J[ti-1][tj-1][shell-1] = 2*Jr;
		J[tj-1][ti-1][shell-1] = 2*Jr;
    };

	fclose(fop);	


}

void ini_alloy(){

	int i,j,k,t,t2,N5;
	int cnt[5] = {0};
// initialize alloy atom species with equal probabilities.
	N5 = N*N*N/5;
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			for(k=0; k<N; k++){		
				t = (int)(randd1()*5);
				Type[i*N_2+j*N+k] = t;
				cnt[t]++;
			}
		for(t = 0 ; t < 5; t++){
			while(cnt[t] > N5){

				do{
					i = (int) (randd1()*N);
					j = (int) (randd1()*N);
					k = (int)(randd1()*N);	
				}while(Type[i*N_2+j*N+k] != t);

				do{
					t2 = (int)(randd1()*5);
				}while(t2 <= t);
				Type[i*N_2+j*N+k] = t2;

				cnt[t]--;
				cnt[t2]++;
			}

			while(cnt[t] < N5){

				do{
					i = (int) (randd1()*N);
					j = (int) (randd1()*N);
					k = (int)(randd1()*N);	
				}while(Type[i*N_2+j*N+k] <= t);
				t2 = Type[i*N_2+j*N+k]; 				
				Type[i*N_2+j*N+k] = t;

				cnt[t]++;
				cnt[t2]--;
			}
		}
	assert(cnt[0] == N5);assert(cnt[1] == N5);assert(cnt[2] == N5);assert(cnt[3] == N5);		
/////////////////////////////////////////
}

void ini_conf(int state){
	
	int i,j,k,t,t2,N4;
	double eta1,eta2,etasq;
	double x, y,z;


	if(state==ORDER)
	{
// ordered initial state;ordered initial postion;
		for(i=0; i<N; i++)
			for(j=0; j<N; j++)
				for(k=0; k<N; k++){
					S[i*N_2+j*N+k] = 1;
					S[N_3+i*N_2+j*N+k] = 0;
					S[2*N_3+i*N_2+j*N+k] = 0;

				}

	}
	else
	{
//disordered inital state;  
		for(i=0; i<N; i++)
			for(j=0; j<N; j++)
				for(k=0; k<N; k++){
			
				   do{
						eta1 = 1.0 - 2.0*randd1();
						eta2 = 1.0 - 2.0*randd1();
						etasq = eta1*eta1 + eta2*eta2;
				   }while(etasq >1.0);
					//these are the new unit vectors
					x = 2.0*eta1*sqrt(1.0-etasq);
					y = 2.0*eta2*sqrt(1.0-etasq);
					z = 1.0 - 2.0*etasq;				

					S[i*N_2+j*N+k] = x;
					S[N_3+i*N_2+j*N+k] = y;
					S[2*N_3+i*N_2+j*N+k] = z;

			}
	
	}
	
}
inline void noffset(int i, int j, int k, int offi, int offj, int offk, int*nn, int cnt){
	int si, sj,sk;
	si = i + offi;
	if(si < 0) si += N;
	if(si >= N) si -= N;
	sj = j + offj;
	if(sj < 0) sj += N;
	if(sj >= N) sj -= N;
	sk = k + offk;
	if(sk < 0) sk += N;
	if(sk >= N) sk -= N;
	nn[cnt++] = si; nn[cnt++] = sj; nn[cnt++] =sk;
	
	si = i - offi;
	if(si < 0) si += N;
	if(si >= N) si -= N;
	sj = j - offj;
	if(sj < 0) sj += N;
	if(sj >= N) sj -= N;
	sk = k - offk;
	if(sk < 0) sk += N;
	if(sk >= N) sk -= N;
	nn[cnt++] = si; nn[cnt++] = sj; nn[cnt++] =sk;
}
inline void  neighbor(int i, int j, int k, int* nn){
	int ii, jj, kk, cnt=0;
// 1st 12;
	noffset(i,j,k,1,0,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,1,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,0,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-1,0, nn, cnt); cnt+=6;
	noffset(i,j,k,1,0,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,0,1,-1, nn, cnt); cnt+=6;
//2nd 6;
	noffset(i,j,k,1,1,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-1,1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,1,1, nn, cnt); cnt+=6;
//3rd 24;
	noffset(i,j,k,1,1,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,1,0,1, nn, cnt); cnt+=6;
	noffset(i,j,k,0,1,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,1,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-2,1, nn, cnt); cnt+=6;	
	noffset(i,j,k,-2,1,1, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,-1,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,0,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,2,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,-1,2,0, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,0,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,-1,2, nn, cnt); cnt+=6;	
//4th 12
	noffset(i,j,k,2,0,0, nn, cnt); cnt+=6;
	noffset(i,j,k,0,2,0, nn, cnt); cnt+=6;
	noffset(i,j,k,0,0,2, nn, cnt); cnt+=6;
	noffset(i,j,k,2,-2,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,0,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,0,2,-2, nn, cnt); cnt+=6;	
// 5th 24
	noffset(i,j,k,2,-1,1, nn, cnt); cnt+=6;
	noffset(i,j,k,2,1,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,2,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,-1,2,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-1,2, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,1,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,-2,1, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,1,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,-2,2,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,2,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,-2,1,2, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-2,2, nn, cnt); cnt+=6;	
// 6th 8
	noffset(i,j,k,1,1,1, nn, cnt); cnt+=6;
	noffset(i,j,k,3,-1,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,3,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,-1,3, nn, cnt); cnt+=6;
// 7th 48
	noffset(i,j,k,2,1,0, nn, cnt); cnt+=6;
	noffset(i,j,k,0,2,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,0,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,0,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,2,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,1,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,-2,0, nn, cnt); cnt+=6;
	noffset(i,j,k,3,0,-2, nn, cnt); cnt+=6;	
	noffset(i,j,k,-2,3,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,3,-2, nn, cnt); cnt+=6;	
	noffset(i,j,k,-2,0,3, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,-2,3, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,-1,-2, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,-2,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,3,-2, nn, cnt); cnt+=6;	
	noffset(i,j,k,-2,3,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,-2,3, nn, cnt); cnt+=6;
	noffset(i,j,k,-2,-1,3, nn, cnt); cnt+=6;
	noffset(i,j,k,3,-1,0, nn, cnt); cnt+=6;
	noffset(i,j,k,3,0,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,-1,3,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,3,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,-1,0,3, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,-1,3, nn, cnt); cnt+=6;
// 8th 6
	noffset(i,j,k,2,-2,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,2,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,-2,2,2, nn, cnt); cnt+=6;	
// 9th 36
/*	noffset(i,j,k,2,2,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,2,-1,2, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,2,2, nn, cnt); cnt+=6;
	noffset(i,j,k,3,-2,1, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,1,-2, nn, cnt); cnt+=6;
	noffset(i,j,k,1,3,-2, nn, cnt); cnt+=6;	
	noffset(i,j,k,-2,3,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-2,3, nn, cnt); cnt+=6;
	noffset(i,j,k,-2,1,3, nn, cnt); cnt+=6;
	noffset(i,j,k,2,2,-3, nn, cnt); cnt+=6;
	noffset(i,j,k,2,-3,2, nn, cnt); cnt+=6;
	noffset(i,j,k,-3,2,2, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,0,0, nn, cnt); cnt+=6;
	noffset(i,j,k,0,3,0, nn, cnt); cnt+=6;
	noffset(i,j,k,0,0,3, nn, cnt); cnt+=6;
	noffset(i,j,k,3,-3,0, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,0,-3, nn, cnt); cnt+=6;	
	noffset(i,j,k,0,3,-3, nn, cnt); cnt+=6;	
// 10th 24
	noffset(i,j,k,3,1,-1, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,-1,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,3,-1, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,3,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-1,3, nn, cnt); cnt+=6;
	noffset(i,j,k,-1,1,3, nn, cnt); cnt+=6;
	noffset(i,j,k,3,-3,1, nn, cnt); cnt+=6;	
	noffset(i,j,k,3,1,-3, nn, cnt); cnt+=6;
	noffset(i,j,k,-3,3,1, nn, cnt); cnt+=6;
	noffset(i,j,k,1,3,-3, nn, cnt); cnt+=6;
	noffset(i,j,k,-3,1,3, nn, cnt); cnt+=6;
	noffset(i,j,k,1,-3,3, nn, cnt); cnt+=6;	*/
}

double Etot(){
	int i,j,k;
	double E = 0.0;

	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			for(k=0; k<N; k++)
				E += 0.5*Esite(i,j,k);

		return E;

}

double Esite(int i, int j, int k){

	int ii,cnt,ti,tj,x,y,z, shell;
	double sx, sy, sz;
	double E = 0.0;
	int nn[600]; // 10 shells contains 200 neighbors for fcc; 


		sx = S[i*N_2+j*N+k] ; 
		sy = S[N_3+i*N_2+j*N+k] ; 
		sz = S[2*N_3+i*N_2+j*N+k] ; 

		ti = Type[i*N_2+j*N+k];

		neighbor(i,j,k,nn);
		cnt = 0; 
		for(shell =0; shell < SH; shell++){
			for(ii=0;ii<NS[shell];ii++){
				x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
				tj = Type[x*N_2+y*N+z];
				E += -J[ti][tj][shell]*( sx*S[x*N_2+y*N+z]
					  +sy*S[N_3+x*N_2+y*N+z]
					  +sz*S[2*N_3+x*N_2+y*N+z] );
			}
		}


		return E;
}


void Rot(){
	int i, j, k, x, y, z, ii, jj, kk;
	double dr, eta1, eta2, etasq, deltaE, E1, E2;
	double rx, ry, rz, c, s, u, sx,sy,sz,im,jm;
	double drx, dry, drz, J;

//	for(i=0; i<N; i++)
//		for(j=0; j<N; j++)
//			for(k=0; k<N; k++){ 
		i = (int) (randd1()*N);
		j = (int) (randd1()*N);
		k = (int)(randd1()*N);				
			   dr=(2.0*randd1()-1.0)*D;
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
		//		c=cos(dr);
		//		s=sin(dr);
		//		u=1.0-c;		  
		
				sx = S[i*N_2+j*N+k];
				sy = S[N_3+i*N_2+j*N+k];
				sz = S[2*N_3+i*N_2+j*N+k];

				E1 = Esite(i,j,k);
				//Perform the rotation on the selected spin
				S[i*N_2+j*N+k] = rx; //(u*rx*rx + c)*sx     + (u*ry*rx - s*rz)*sy  + (u*rz*rx + ry*s)*sz ;
				S[N_3+i*N_2+j*N+k] = ry; //(u*rx*ry + rz*s)*sx  + (u*ry*ry + c)*sy     + (u*rz*ry - rx*s)*sz ;
				S[2*N_3+i*N_2+j*N+k] = rz; //(u*rx*rz - ry*s)*sx  + (u*ry*rz + rx*s)*sy  + (u*rz*rz + c)*sz ;	

				//evaluate energy change
				E2 = Esite(i,j,k);
				deltaE = E2 - E1;

				attd++;
				if(WangLandau(currEtot, currEtot+deltaE) == 1){// accept
					currEtot += deltaE; accd++;
				}else{//reject
					S[i*N_2+j*N+k] = sx;
					S[N_3+i*N_2+j*N+k] = sy;
					S[2*N_3+i*N_2+j*N+k] = sz;
				}

}

void wolff(int i, int j, int k, double rx, double ry, double rz){

        double sx,sy,sz,pi,pj, delta;
        int ti,tj, cnt, ii, shell, x, y, z, nn[600];

        cluster[i*N_2+j*N+k] = true;
        sx = S[i*N_2+j*N+k];
        sy = S[N_3+i*N_2+j*N+k];
        sz = S[2*N_3+i*N_2+j*N+k];
        pi = (sx*rx + sy*ry + sz*rz);
        ti = Type[i*N_2+j*N+k];

        S[i*N_2+j*N+k] = sx - 2*pi*rx;
        S[N_3+i*N_2+j*N+k] = sy - 2*pi*ry;
        S[2*N_3+i*N_2+j*N+k] = sz - 2*pi*rz;

                neighbor(i,j,k,nn);
                cnt = 0;
                for(shell =0; shell < SH; shell++){
                        for(ii=0;ii<NS[shell];ii++){
                                x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
                                if(!cluster[x*N_2+y*N+z]){
                                        tj = Type[x*N_2+y*N+z];
                                        pj = (S[x*N_2+y*N+z]*rx + S[N_3+x*N_2+y*N+z]*ry + S[2*N_3+x*N_2+y*N+z]*rz);
                                        delta = -2*J[ti][tj][shell]/(T_scale*pT)*pi*pj;
                                        if(delta < 0 && randd1() < (1-exp(delta)) )     
                                                wolff(x,y,z,rx,ry,rz);
                                }
                        }
                }
}


void write_pos(){
  char s[512];
  FILE *ofp;
  sprintf(s,"compos%d.dat",myrank);
  ofp=fopen(s,"w");
	int i,j,k;
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			for(k=0;k<N;k++){
				fprintf(ofp, "%d\t%d\t%d\t%d\t \n",i+k, i+j, j+k,Type[i*N_2+j*N+k]);
			}
	fclose(ofp);
}

void write_mol2(int frame)
{
  int i,j,k,cnt;
  char s[512];
  double M[6];
  FILE *ofp;

  Mag(M);
  sprintf(s,"snap_rank%d_fr%d.mol2",myrank,frame-10*myrank);
  ofp=fopen(s,"w");

  fprintf(ofp,"@<TRIPOS>MOLECULE\n");
  fprintf(ofp,"E=%g M=%g M1=%g M2=%g M3=%g M4=%g M5=%g\n", currEtot,M[5],M[0],M[1],M[2],M[3],M[4]);
  fprintf(ofp," %d %d\n",N*N*N*2,N*N*N);
  fprintf(ofp,"SMALL\n");
  fprintf(ofp,"NO_CHARGES\n");
  fprintf(ofp,"@<TRIPOS>ATOM\n");



  cnt = 0;
  for(i=0;i<N;i++)
	  for(j=0; j<N;j++)
		  for(k=0;k<N;k++)
    {
		cnt++;
      fprintf(ofp,"%d H %d %d %d \n",cnt,i+k, i+j, j+k);
      cnt++;
	  fprintf(ofp,"%d %c %g %g %g \n",cnt, atom[i*N_2+j*N+k],i+k+S[i*N_2+j*N+k], i+j+S[N_3+i*N_2+j*N+k], j+k+S[2*N_3+i*N_2+j*N+k] );

		  };
  fprintf(ofp,"@<TRIPOS>BOND\n");


  for(i=0;i<N*N*N;++i)
    {
      fprintf(ofp,"%d %d %d %d\n",i+1,i*2+1,i*2+2,1);
    };

  fclose(ofp);


}


//Calculates thermodynamic quantities and writes modification-factor labeled files
void Mag(double* M){
	int i, j, k, t;
	double Mx[6], My[6], Mz[6];

	for(i=0;i<6;i++){
		Mx[i]=My[i]=Mz[i] = 0.0;
	}
	for( i=0; i<N; i++)
		for( j=0; j<N; j++)
			for( k=0; k<N; k++){
				t = Type[i*N_2+j*N+k];
				Mx[t] += S[i*N_2+j*N+k];
				My[t] += S[N_3+i*N_2+j*N+k];
				Mz[t] += S[2*N_3+i*N_2+j*N+k];
				Mx[5] += S[i*N_2+j*N+k];
				My[5] += S[N_3+i*N_2+j*N+k];
				Mz[5] += S[2*N_3+i*N_2+j*N+k];
			}
	for(i = 0; i<6; i++){
		M[i]= sqrt(Mx[i]*Mx[i] +My[i]*My[i]+Mz[i]*Mz[i]);
		if(i < 5 )
			M[i] *= invN*5;
		else 
			M[i] *= invN;
	}
}




