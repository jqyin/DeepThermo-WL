#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <chrono>
#include <map>
#include "mpi.h"
#include "random.h"
#include "pt.hpp"
#include "alloy.hpp"
#include "rand.hpp"
//#define RESTART

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
        int i;

	Atom=(short*)malloc(sizeof(short)*N_3);
	Atomo=(short*)malloc(sizeof(short)*N_3);
	//Atom=(uint8_t*)malloc(sizeof(uint8_t)*N_3);
	cluster = (bool*)malloc(sizeof(bool)*N_3);
	att=acc=0;
	init_genrand(myrank*rand());
 
	initialize();
#ifdef RESTART
        if(Restart == 1)
	        read_state();
#endif
}

void swap(bool even){
	int i,j,k, itag=123;
	double delta, Ei, Ej, Ti, Tj;
	bool flag;
	MPI_Status istatus;
	unsigned mem_u = N_3*sizeof(short);
	short* buf_ui = (short*) malloc(mem_u);
	short* buf_uj = (short*) malloc(mem_u);

	if(even){
		if(myrank%2 == 0){
			Ei = currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank+1, itag,
									&Ej,1,MPI_DOUBLE, myrank+1,itag,
									MPI_COMM_WORLD,&istatus);
			delta = (Ei-Ej)*E_scale*(1/(T_scale*T[myrank+1])-1/(T_scale*T[myrank]));
			if(delta <=0 || genrand_real2() < exp(-delta))
				flag = true;
			else
				flag = false;
			MPI_Send(&flag, 1, MPI_LOGICAL,myrank+1,itag, MPI_COMM_WORLD);
			att++;
			if(flag){
				memcpy(&buf_ui[0], Atom, mem_u);
				MPI_Sendrecv(buf_ui, N_3, MPI_SHORT, myrank+1, itag,
										buf_uj,N_3,MPI_SHORT, myrank+1,itag,
										MPI_COMM_WORLD,&istatus);

				memcpy(Atom, &buf_uj[0], mem_u);

				currEtot = Ej;
				acc ++;
			}
		}else{
			Ei = currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank-1, itag,
									&Ej,1,MPI_DOUBLE, myrank-1,itag,
									MPI_COMM_WORLD,&istatus);
			MPI_Recv(&flag, 1, MPI_LOGICAL,myrank-1,itag, MPI_COMM_WORLD, &istatus);	
			if(flag){
				memcpy(&buf_ui[0], Atom, mem_u);
				MPI_Sendrecv(buf_ui, N_3, MPI_SHORT, myrank-1, itag,
										buf_uj,N_3,MPI_SHORT, myrank-1,itag,
										MPI_COMM_WORLD,&istatus);

				memcpy(Atom, &buf_uj[0], mem_u);

				currEtot = Ej;

			}
		}	

	}else{
		if(myrank%2 != 0 && myrank!= nprocs-1){
			Ei = currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank+1, itag,
									&Ej,1,MPI_DOUBLE, myrank+1,itag,
									MPI_COMM_WORLD,&istatus);
			delta = (Ei-Ej)*E_scale*(1/(T_scale*T[myrank+1])-1/(T_scale*T[myrank]));
			if(delta <=0 || genrand_real2() < exp(-delta))
				flag = true;
			else
				flag = false;
			MPI_Send(&flag, 1, MPI_LOGICAL,myrank+1,itag, MPI_COMM_WORLD);
			att ++;
			if(flag){
				memcpy(&buf_ui[0], Atom, mem_u);
                                MPI_Sendrecv(buf_ui, N_3, MPI_SHORT, myrank+1, itag,
                                                                                  buf_uj,N_3,MPI_SHORT, myrank+1,itag,
                                                                                  MPI_COMM_WORLD,&istatus);

                                memcpy(Atom, &buf_uj[0], mem_u);

				currEtot = Ej;

				acc++;
			}
		}else if(myrank != 0 && myrank != nprocs-1){
			Ei = currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, myrank-1, itag,
									&Ej,1,MPI_DOUBLE, myrank-1,itag,
									MPI_COMM_WORLD,&istatus);
			MPI_Recv(&flag, 1, MPI_LOGICAL,myrank-1,itag, MPI_COMM_WORLD, &istatus);	
			if(flag){
				memcpy(&buf_ui[0], Atom, mem_u);
                                MPI_Sendrecv(buf_ui, N_3, MPI_SHORT, myrank-1, itag,
                                                                                 buf_uj,N_3,MPI_SHORT, myrank-1,itag,
                                                                                 MPI_COMM_WORLD,&istatus);

                                memcpy(Atom, &buf_uj[0], mem_u);

				currEtot = Ej;


			}
		}		

	}


	ini_W();
	//MPI_Barrier(MPI_COMM_WORLD);

	free(buf_ui);
	free(buf_uj);

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

      R=exp(-(Ef-Ei)*E_scale/(T_scale*pT));

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

void mchybrid(SamplingMode mode){

	int cnt, i, j,k;
	double eta1, eta2, etasq, rx, ry, rz;
	for(cnt=0; cnt<N*N*N; cnt++){
		BondSwap(mode);
	}
	//vae_update(mode);
/*	for(i = 0; i < N; i++)
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
*/
}
/*
void correlation(double* avgScorr){
	int cnt,i,j,k,shell,ii;
	int nn[MAX_NEIGHBORS*3]; // 10 shells contains 200 neighbors for fcc; 
	double sx,sy,sz,Scorr, spin_corr[SH];
	int x,y,z;
	for(shell=0; shell < SH; shell++)
		spin_corr[shell] = 0.0;
	for(i = 0; i < N; i++)for(j = 0; j < N; j++)for(k = 0; k < N; k++){
		neighbor(i,j,k,nn);
		cnt = 0; 
		sx = S[i*N_2+j*N+k];	
		sy = S[N_3+i*N_2+j*N+k];	
		sz = S[2*N_3+i*N_2+j*N+k];	
		for(shell =0; shell < SH; shell++){
			Scorr = 0.0;
			for(ii=0;ii<NS[shell];ii++){
				x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
				Scorr += sx*S[x*N_2+y*N+z]+
						    sy*S[N_3+x*N_2+y*N+z]+
						    sz*S[2*N_3+x*N_2+y*N+z];
			}
			spin_corr[shell] += Scorr/NS[shell];
		}
	}
	for(shell=0; shell < SH; shell++){
		spin_corr[shell] /= N_3;
		avgScorr[shell] += spin_corr[shell];
	}
}*/

void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int irun, SamplingMode mode){
	int mcs,i,j,k,n, t;
	FILE *ofp, *ofp1, *ofp2, *ofp_time;
	double avgE, avgE2, M, avgM, avgM2, avgM4;
	double c, x, bc, avgScorr[SH];
	double* Et = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Mt = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Ea = (double*)malloc(sizeof(double)*nT);
	double* Ma = (double*)malloc(sizeof(double)*nT);
	double* ca = (double*)malloc(sizeof(double)*nT);
	double* xa = (double*)malloc(sizeof(double)*nT);
	double* bca = (double*)malloc(sizeof(double)*nT);
	double* Scorra = (double*)malloc(sizeof(double)*nT*SH);
	int* attempts=(int*)malloc(sizeof(int)*nT);
	int* accepts=(int*)malloc(sizeof(int)*nT);
	double var =0.0;
	char s[512];
	bool flag = true, stop = false;
	time_t t1,t2;
	std::chrono::high_resolution_clock::time_point t3,t4;
	std::map<double, int> histE, histM, histM2, histM4;
	avgE = avgE2 = 0.0;
	avgM = avgM2 = avgM4 =0.0;
	for(i=0; i<SH; i++)
		avgScorr[i] = 0.0;
	if(myrank ==0)
		t1 = time(NULL);

//#ifndef RESTART
        if(Restart == 0){
	        for(mcs=0;mcs<DROPI;mcs++){
		        mchybrid(mode);
		        if(mcs%2==0){ 
			        swap(flag);
			        if(flag) 
				        flag = false;
			        else
				        flag = true;
		        }
#ifdef Time_Series
		        Et[mcs] = currEtot;
	        	M = L1(); 
		        Mt[mcs] = M[NE];
#endif
	        }
        }

//#endif

	mcs = 0;
//#ifdef RESTART
        if(Restart == 1){
	        sprintf(s, "mc%d.input", myrank);
	        ofp1 = fopen(s,"rb");
	        if(ofp1 == NULL){
		        printf("mc.input no found\n");
		        exit(-3);
	        }
		fread(&mcs, sizeof(int), 1, ofp1);
		fread(&avgE, sizeof(double), 1, ofp1);
		fread(&avgE2, sizeof(double), 1, ofp1);
		fread(&avgM, sizeof(double), 1, ofp1);
		fread(&avgM2, sizeof(double), 1, ofp1);
		fread(&avgM4, sizeof(double), 1, ofp1);
	        fclose(ofp1);
        }
//#endif

	while(mcs < SAMPS){
		for(i=0;i<SEP;++i){
			/*if(mcs == 1){
				if(myrank == 0){
                			ofp_time= fopen("infer.dat","w");
                			//t3 = time(NULL);
					t3 = std::chrono::high_resolution_clock::now();
        			}
			}*/
			mchybrid(mode);
			/*if(mcs == 1){
				if(myrank == 0){
					//t4 = time(NULL);
					t4 = std::chrono::high_resolution_clock::now();
					std::chrono::duration<double, std::milli> ms_double = t4 - t3;
					fprintf(ofp_time, "infers/ms: %f\n", 1.0*N*N*N/ms_double.count());
					fflush(ofp_time);
				}

			}*/
			if( int(mcs*SEP+i)%2==0 ){	
				swap(flag);
				if(flag) 
					flag = false;
				else
					flag = true;
			}
		}

		avgE +=currEtot;
		avgE2 +=currEtot*currEtot;
		M = L1();
		avgM += M;
		avgM2 += M*M;
		avgM4 += M*M;
/////////////////////////////////////////////////////////
		if(mcs % CHPT_STEPS == 0){ //checkpoint

			printf("myrank =%d, currEtot = %g, Eng = %g\n", myrank, currEtot, Etot());
			if(myrank==0){
				t2 = time(NULL);
				if(1.0*(t2-t1)/3600 >  0.9*TIMER)
					stop = true;
			}
			MPI_Bcast(&stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
			if(stop){

				write_state();
				sprintf(s, "mc%d.input", myrank);
				ofp1 = fopen(s, "wb");
				mcs++; 
				fwrite(&mcs, sizeof(int), 1, ofp1);
				fwrite(&avgE, sizeof(double), 1, ofp1);
				fwrite(&avgE2, sizeof(double), 1, ofp1);
				fwrite(&avgM, sizeof(double), 1, ofp1);
				fwrite(&avgM2, sizeof(double),1, ofp1);
				fwrite(&avgM4, sizeof(double),1, ofp1);
				fclose(ofp1);
				printf("avgE2-avgE*avgE %g\n",avgE2/mcs-avgE*avgE/mcs/mcs);
				MPI_Finalize();
				exit(123);
			}
		}
///////////////////////////////////////////////////////////////////

#ifdef HIST
		++histE[currEtot];
		++histM[M[NE]];
#endif

#ifdef Time_Series
		Mt[mcs] = M[NE];
		Et[mcs] = currEtot;
#endif
		mcs++;
	};
	// output spin configuration;
	if(myrank == 0)
		write_pos(); 
	//write_xyz(SAMPS);
	avgE/=1.0*SAMPS;
        avgE2/=1.0*SAMPS;
	avgM/=1.0*SAMPS;
        avgM2/=1.0*SAMPS;
        avgM4/=1.0*SAMPS;
	c = (avgE2 - avgE*avgE)*E_scale*E_scale/(T_scale*T_scale*pT*pT) ;
	//for(i=0;i<SH;i++)
	//  	avgScorr[i]/=1.0*SAMPS;
	x = (avgM2 - avgM*avgM)*N_3/pT/T_scale;
        bc = 1 - avgM4/avgM2/avgM2/3;

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
	}
#endif

	MPI_Gather(&avgE, 1, MPI_DOUBLE, Ea, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&avgM, 1, MPI_DOUBLE, Ma, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&c, 1, MPI_DOUBLE, ca, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	//MPI_Gather(avgScorr, SH, MPI_DOUBLE, Scorra, SH, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&x, 1, MPI_DOUBLE, xa, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&bc, 1, MPI_DOUBLE, bca, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&att, 1, MPI_INT, attempts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	MPI_Gather(&acc, 1, MPI_INT, accepts, 1, MPI_INT, 0 , MPI_COMM_WORLD);

	if(myrank == 0){
		sprintf(s, "stat%d.dat", irun);
		ofp = fopen(s,"w");
		for(n=0;n<nT;n++){
		  fprintf(ofp,"%g\t%g\t%g\t",
			  T[n],Ea[n]*invN,ca[n]*invN
			 /* Ma[n*5+3], xa[n*5+3], bca[n*5+3]*/);
		  fprintf(ofp, "%g\t%g\t%g\t",Ma[n+i],xa[n*+i],bca[n+i]);
		  fprintf(ofp,"\n");
		}	
		fclose(ofp);
		  
		/*sprintf(s, "Scorr%d.dat", irun);
		ofp = fopen(s,"w");
		fprintf(ofp, "r(a)\\T(K)\t");
		for(n=0;n<nT;n++){
			fprintf(ofp, "%g\t", T[n]);
		}
		fprintf(ofp, "\n");
		for(i=0; i<SH;i++){
			fprintf(ofp, "%g\t",Dist[i]);
			for(n=0;n<nT;n++){
				fprintf(ofp, "%g\t", Scorra[n*SH+i]);
			}
			fprintf(ofp, "\n");
		}*/ 

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
		fflush(ofp2);
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
	int t;
	//free(Atom);
	free(cluster);
	free(T);
	for(t=0; t<NE; t++)
		free(inputPos[t]);
}


void write_state(){

	char s[512];
	sprintf(s, "state%d.input", myrank);
	FILE* fp = fopen(s, "wb");

	fwrite(Atom, sizeof(short), N_3, fp);
	fwrite(NT, sizeof(int), NE, fp);
	fwrite(&currEtot, sizeof(double), 1, fp);
	fclose(fp);

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
	fread(Atom, sizeof(short), N_3, fp);
	fread(NT, sizeof(int), NE, fp);
	fread(&eng,sizeof(double), 1, fp);
        ini_apos();
	currEtot = Etot();
	if(fabs(eng-currEtot)>1e-5){
		printf("state file corrupted!\n");
		exit(-2);
	}

	fclose(fp);

}
