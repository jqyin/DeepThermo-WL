#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <map>
#include "mpi.h"
#include "random.h"
#include "pt.hpp"
#include "alloy.hpp"
#include "rand.hpp"

//#define RESTART
//#define Time_Series
//#define HIST


extern int myrank, nprocs;

void ini_T(double Ti, double Tf, int nT){
	int i;

	if(Ti*Tf > 1e-5){
		T = (double*)malloc(sizeof(double)*nT);
		if(T != NULL){
			T[0] = log(Ti);
			T[nT-1]=log(Tf);
			for(i=1;i<nT-1;i++){
				T[i] = T[0] +(T[nT-1]-T[0])*i/(nT-1);
			}
			for(i=0;i<nT;i++)
				T[i]=exp(T[i]);
		}
		else
			exit(1);
	}else{
		exit(1);
	}

	pT = T[myrank];

}

void ini_sys(){

		N_2 = N*N;
		N_3 = N*N*N;
       int i,j,k;
	   unsigned mem_s = 3*N_3*sizeof(double);

        S= (double*)malloc(mem_s);
		Type=(int*)malloc(sizeof(int)*N_3);
		cluster = (bool*)malloc(sizeof(bool)*N_3);
		atom = (char*) malloc(sizeof(char)*N_3);
		att=acc=0;
		init_genrand(myrank*rand() );
 
		if(S!=NULL){
				initialize(DISORDER);
 
		}else
			exit(1);
//sync ini conf;
		MPI_Bcast(S, 3*N_3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(Type, N_3, MPI_INT, 0, MPI_COMM_WORLD);
	  for(i=0;i<N;i++)
		  for(j=0; j<N;j++)
			  for(k=0;k<N;k++){
				if(Type[i*N_2+j*N+k] == 0 ) // cr
					atom[i*N_2+j*N+k] = 'O';
				else if(Type[i*N_2+j*N+k] == 1 ) 
					atom[i*N_2+j*N+k] = 'N';
				else if(Type[i*N_2+j*N+k] == 2 )
					atom[i*N_2+j*N+k] = 'C';
				else if(Type[i*N_2+j*N+k] == 3 )
					atom[i*N_2+j*N+k] = 'S';
				else 
					atom[i*N_2+j*N+k] = 'P';
			  }
		currEtot = Etot();
#ifdef RESTART
		read_state();
#endif
}

void swap(bool even){
	int i,j,k, itag=123;
	double delta, Ei, Ej, Ti, Tj;
	bool flag;
	MPI_Status istatus;
	unsigned mem_s = N_3*sizeof(double);
	double* bufi = (double*) malloc(3*mem_s);
	double* bufj = (double*) malloc(3*mem_s);

	if(even){

			if(myrank%2 == 0){
				Ei = pE;
				MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank+1, itag,
										&Ej,1,MPI_DOUBLE, myrank+1,itag,
										MPI_COMM_WORLD,&istatus);
				delta = (Ei-Ej)*(1/(T_scale*T[myrank+1])-1/(T_scale*T[myrank]));
				if(delta <=0 || genrand_real2() < exp(-delta))
					flag = true;
				else
					flag = false;
				MPI_Send(&flag, 1, MPI_LOGICAL,myrank+1,itag, MPI_COMM_WORLD);
				att++;
				if(flag){
					memcpy(&bufi[0], S, 3*mem_s);
					MPI_Sendrecv(bufi, 3*N_3, MPI_DOUBLE, myrank+1, itag,
											bufj,3*N_3,MPI_DOUBLE, myrank+1,itag,
											MPI_COMM_WORLD,&istatus);

					memcpy(S, &bufj[0], 3*mem_s);

					pE = Ej;
					acc ++;
				}
			}else{
				Ei = pE;
				MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank-1, itag,
										&Ej,1,MPI_DOUBLE, myrank-1,itag,
										MPI_COMM_WORLD,&istatus);
				MPI_Recv(&flag, 1, MPI_LOGICAL,myrank-1,itag, MPI_COMM_WORLD, &istatus);	
				if(flag){
					memcpy(&bufi[0], S, 3*mem_s);
					MPI_Sendrecv(bufi, 3*N_3, MPI_DOUBLE, myrank-1, itag,
											bufj,3*N_3,MPI_DOUBLE, myrank-1,itag,
											MPI_COMM_WORLD,&istatus);

					memcpy(S, &bufj[0], 3*mem_s);

					pE = Ej;

				}
			}	

	}else{

			if(myrank%2 != 0 && myrank!= nprocs-1){
				Ei = pE;
				MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank+1, itag,
										&Ej,1,MPI_DOUBLE, myrank+1,itag,
										MPI_COMM_WORLD,&istatus);
				delta = (Ei-Ej)*(1/(T_scale*T[myrank+1])-1/(T_scale*T[myrank]));
				if(delta <=0 || genrand_real2() < exp(-delta))
					flag = true;
				else
					flag = false;
				MPI_Send(&flag, 1, MPI_LOGICAL,myrank+1,itag, MPI_COMM_WORLD);
				att ++;
				if(flag){
					memcpy(&bufi[0], S, 3*mem_s);
					MPI_Sendrecv(bufi, 3*N_3, MPI_DOUBLE, myrank+1, itag,
											bufj,3*N_3,MPI_DOUBLE, myrank+1,itag,
											MPI_COMM_WORLD,&istatus);

					memcpy(S, &bufj[0], 3*mem_s);


					pE = Ej;


					acc++;
				}
			}else if(myrank != 0 && myrank != nprocs-1){
				Ei = pE;
				MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank-1, itag,
										&Ej,1,MPI_DOUBLE, myrank-1,itag,
										MPI_COMM_WORLD,&istatus);
				MPI_Recv(&flag, 1, MPI_LOGICAL,myrank-1,itag, MPI_COMM_WORLD, &istatus);	
				if(flag){

					memcpy(&bufi[0], S, 3*mem_s);

					MPI_Sendrecv(bufi, 3*N_3, MPI_DOUBLE, myrank-1, itag,
											bufj,3*N_3,MPI_DOUBLE, myrank-1,itag,
											MPI_COMM_WORLD,&istatus);

					memcpy(S, &bufj[0], 3*mem_s);

					pE = Ej;


				}
			}		

	}


	//MPI_Barrier(MPI_COMM_WORLD);

	free(bufi);
	free(bufj);

}
int Metropolis(double Ei, double Ef)
{
  double R;

  if(Ef<Ei)
    {
      //accept
      return 1;
    }
  else
    {

		R=exp(-(Ef-Ei)/(T_scale*pT));

      if(randd1()<R)
	{
	  //accept
	  return 1;
	}
      else
	{
	  //reject
	  return 0;
	};
    };
}

void mchybrid(){

	int cnt, i, j,k;
	double eta1, eta2, etasq, rx, ry, rz;
/*	for(cnt=0; cnt<N*N*N; cnt++){
			Rot();
	}*/
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			for(k = 0; k < N; k++)
				cluster[i*N_2+j*N+k] = false;

		i = (int) (randd1()*N);
		j = (int) (randd1()*N);
		k = (int)(randd1()*N);		
	do{
		eta1 = 1.0 - 2.0*randd1();
		eta2 = 1.0 - 2.0*randd1();
		etasq = eta1*eta1 + eta2*eta2;
	 }while(etasq >1.0);
	//these are the new unit vectors
	rx = 2.0*eta1*sqrt(1.0-etasq);
	ry = 2.0*eta2*sqrt(1.0-etasq);
	rz = 1.0 - 2.0*etasq;	

	wolff(i,j,k, rx,ry,rz);
	currEtot = Etot();
}
void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int irun){
	int mcs,i,j,k,n, t;
	 FILE *ofp, *ofp1;
	 FILE * ofp2;
	double Mx[6], My[6], Mz[6]; // 0 - cr; 1 - Mn, 2-fe; 3-co; 4-ni; 5-total;
	double M[6], M2[6], M4[6];
	double avgE, avgE2, avgM[6], avgM2[6], avgM4[6];
	double c, x[6], bc[6];
	double* Et = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Mt = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Ea = (double*)malloc(sizeof(double)*nT);
	double* Ma = (double*)malloc(sizeof(double)*nT*6);
	double* ca = (double*)malloc(sizeof(double)*nT);
	double* xa = (double*)malloc(sizeof(double)*nT*6);
	double* bca = (double*)malloc(sizeof(double)*nT*6);
	int* attempts=(int*)malloc(sizeof(int)*nT);
	int* accepts=(int*)malloc(sizeof(int)*nT);
	double var =0.0;
	char s[512];
	bool flag = true, stop = false;
	time_t t1,t2;
	std::map<double, int> histE, histM, histM2, histM4;
	avgE = avgE2 = 0.0;
	for(i=0; i<6;i++){
		avgM[i] = avgM2[i] = avgM4[i] =0.0;
	}
 
	if(myrank ==0)
		t1 = time(NULL);

#ifndef RESTART
	for(mcs=0;mcs<DROPI;mcs++){
	//	for(n=0;n<nT;n++)
		mchybrid();
		if(mcs%5==0){ 
			swap(flag);
			if(flag) 
				flag = false;
			else
				flag = true;
		}
#ifdef Time_Series
		Et[mcs] = pE;

				Mx[1] = My[1] = Mz[1] =0.0;
				for( i=0; i<N; i++)
					for( j=0; j<N; j++)
						for( k=0; k<N; k++){
							t = Type[i*N_2+j*N+k];
							Mx[t] += S[i*N_2+j*N+k];
							My[t] += S[N_3+i*N_2+j*N+k];
							Mz[t] += S[2*N_3+i*N_2+j*N+k];
						}
				Mt[mcs]= sqrt(Mx[1]*Mx[1] +My[1]*My[1]+Mz[1]*Mz[1])*invN;
#endif
	}

#endif

	mcs = 0;
#ifdef RESTART
	
		sprintf(s, "mc%d.input", myrank);
		ofp1 = fopen(s,"r");
		if(ofp1 == NULL){
			printf("mc.input no found\n");
			exit(-3);
		}
		fscanf(ofp1,"%d\t%lg\t%lg\t%lg\t%lg\t%lg",
							&mcs, &avgE,&avgE2,&avgM[0],&avgM2[0],&avgM4[0]);
		fclose(ofp1);
	

#endif

	  while(mcs < SAMPS){

			for(i=0;i<SEP;++i){
			//	for(n=0;n<nT;n++)
					mchybrid();
					if( int(mcs*SEP+i)%5==0 ){	
					swap(flag);
					if(flag) 
						flag = false;
					else
						flag = true;
				}
			}

		//	for(n=0;n<nT;n++){							
				avgE +=pE;
				avgE2 +=pE*pE;
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
					M2[i] = M[i]*M[i];
					M4[i] = M2[i]*M2[i];	

					avgM[i] += M[i];
					avgM2[i] += M2[i];
					avgM4[i] += M4[i];
				}
/////////////////////////////////////////////////////////
			if(mcs % 50000 == 0){ //checkpoint

				printf("myrank =%d, pE = %g, Eng = %g\n", myrank, pE, Etot());
				if(myrank==0){
					t2 = time(NULL);
					if(1.0*(t2-t1)/3600 > 23.4)
						stop = true;
				}
				MPI_Bcast(&stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
				if(stop){

						write_state();
						sprintf(s, "mc%d.input", myrank);
						ofp1 = fopen(s, "w");
						fprintf(ofp1, "%d\t%g\t%g\t%g\t%g\t%g",
							mcs, avgE,avgE2,avgM[0],avgM2[0],avgM4[0]);
						fclose(ofp1);


					MPI_Finalize();
					exit(123);
				}
			}
///////////////////////////////////////////////////////////////////

#ifdef HIST
				++histE[pE];
				++histM[M[1]];
				++histM2[M2[1]];
				++histM4[M4[1]];
#endif

#ifdef Time_Series
				Mt[mcs] = M[1];
	//		}
			Et[mcs] = pE;
#endif
			mcs++;
		};

//	for(n=0;n<nT;n++){
	  avgE/=1.0*SAMPS;
      avgE2/=1.0*SAMPS;
	  for(i=0; i<6;i++){
		  avgM[i]/=1.0*SAMPS;
		  avgM2[i]/=1.0*SAMPS;
		  avgM4[i]/=1.0*SAMPS;
	  }
	  c = (avgE2 - avgE*avgE)/(T_scale*T_scale*pT*pT) ;
	  for(i=0;i<6;i++){
		  if(i < 5)
			 x[i] = (avgM2[i] - avgM[i]*avgM[i])*N*N*N/5/pT/T_scale;
		  else
			 x[i] = (avgM2[i] - avgM[i]*avgM[i])*N*N*N/pT/T_scale;

		  bc[i] = 1 - avgM4[i] / (avgM2[i]*avgM2[i]) / 3.0;	 
	  }
//	}


#ifdef Time_Series
	if(myrank < nprocs){
		sprintf(s, "EMt%g.dat", pT);
		ofp1 = fopen(s,"w");
		for(i=0;i<(SAMPS+DROPI);i++)
				fprintf(ofp1,"%d\t%g\t%g\n",i,Et[i],Mt[i]);
		fclose(ofp1);
	}
#endif

#ifdef HIST
	if(myrank < nprocs){
		std::map<double,int>::iterator m_it;

		sprintf(s, "histE%g.dat", pT);
		ofp1 = fopen(s,"w");
		for(m_it = histE.begin(); m_it!= histE.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);

		sprintf(s, "histM%g.dat", pT);
		ofp1 = fopen(s,"w");
		for(m_it = histM.begin(); m_it!= histM.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);

		sprintf(s, "histM2%g.dat", pT);
		ofp1 = fopen(s,"w");
		for(m_it = histM2.begin(); m_it!= histM2.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);

		sprintf(s, "histM4%g.dat", pT);
		ofp1 = fopen(s,"w");
		for(m_it = histM4.begin(); m_it!= histM4.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);
	}
#endif

	  MPI_Gather(&avgE, 1, MPI_DOUBLE, Ea, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	  MPI_Gather(avgM, 6, MPI_DOUBLE, Ma, 6, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	  MPI_Gather(&c, 1, MPI_DOUBLE, ca, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	  MPI_Gather(x, 6, MPI_DOUBLE, xa, 6, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	  MPI_Gather(bc, 6, MPI_DOUBLE, bca, 6, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	  MPI_Gather(&att, 1, MPI_INT, attempts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	  MPI_Gather(&acc, 1, MPI_INT, accepts, 1, MPI_INT, 0 , MPI_COMM_WORLD);

	  if(myrank == 0){
		  sprintf(s, "stat%d.dat", irun);
		  ofp = fopen(s,"w");
		for(n=0;n<nT;n++){
		  fprintf(ofp,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			  T[n],Ea[n]*invN,ca[n]*invN, 
			  Ma[n*6+5], xa[n*6+5], bca[n*6+5],
			  Ma[n*6], xa[n*6], bca[n*6],
			  Ma[n*6+1], xa[n*6+1], bca[n*6+1],
			  Ma[n*6+2], xa[n*6+2], bca[n*6+2],
			  Ma[n*6+3], xa[n*6+3], bca[n*6+3],
			  Ma[n*6+4], xa[n*6+4], bca[n*6+4]);
		}	
		  fclose(ofp);



		  sprintf(s, "misc%d.dat", irun);
		 ofp2 = fopen(s,"w");
	//	 fprintf(ofp2,"Avg_diff:%g  Avg_rot:%g Sample size:%d\n",1.0*acceptdiff/attemptdiff, 1.0*acceptrot/attemptrot,  (int)SAMPS);
		 fprintf(ofp2,"swap prob:\n");
		 for(n=0;n<nT-1;n++){
		 fprintf(ofp2,"%g\t%g\n",T[n], 1.0*accepts[n]/attempts[n]);		
		 }
	  }
	  MPI_Gather(&attd, 1, MPI_INT, attempts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	  MPI_Gather(&accd, 1, MPI_INT, accepts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	  if(myrank == 0){
		 fprintf(ofp2,"Rot prob:\n");
		 for(n=0;n<nT;n++){
		 fprintf(ofp2,"%g\t%g\n",T[n], 1.0*accepts[n]/attempts[n]);		
		 }
	  }
	 free(Ea);
	 free(Ma);
	 free(ca);
	 free(xa);
	 free(bca);
	 free(Et);
	 free(Mt);
	 free(attempts);
	 free(accepts);
}

void freePT(){
	free(S);
	free(Type);
	free(cluster);
	free(T);

}


void write_state(){

	char s[512];
	sprintf(s, "state%d.input", myrank);
	FILE* fp = fopen(s, "wb");

	fwrite(S, sizeof(double), 3*N_3, fp);

	fwrite(&pE, sizeof(double), 1, fp);
	fclose(fp);
	
	if(myrank == 0){
		fp = fopen("type.input", "wb");	
		fwrite(Type, sizeof(int), N_3, fp);
		fclose(fp);
	}

}
void read_state(){

	char s[512];
	double eng;
	sprintf(s, "state%d.input", myrank);
	FILE* fp = fopen(s, "rb");
	if(fp == NULL){
		printf("state file no found\n");
		exit(-1);
	}
	fread(S, sizeof(double), 3*N_3, fp);

	fread(&eng,sizeof(double), 1, fp);
	pE = Etot();
	if(fabs(eng-pE)>1e-5){
		printf("state file corrupted!\n");
		exit(-2);
	}

	fclose(fp);

	fp = fopen("type.input", "wb");	
	fread(Type, sizeof(int), N_3, fp);
	fclose(fp);

}
