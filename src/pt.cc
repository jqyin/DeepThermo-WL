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


void ini_T(double Ti, double Tf, int nT){
	int i;

	if(Ti*Tf > 1e-5){
		ptState.T = (double*)malloc(sizeof(double)*nT);
		if(ptState.T != NULL){
			ptState.T[0] = log(Ti);
			ptState.T[nT-1]=log(Tf);
			for(i=1;i<nT-1;i++){
				ptState.T[i] = ptState.T[0] +(ptState.T[nT-1]-ptState.T[0])*i/(nT-1);
			}
			for(i=0;i<nT;i++)
				ptState.T[i]=exp(ptState.T[i]);
		}
		else
			exit(1);
	}else{
		exit(1);
	}

	ptState.pT = ptState.T[mpiState.myrank];

}

void ini_sys(){

	alloyState.N_2 = alloyState.N*alloyState.N;
	alloyState.N_3 = alloyState.N*alloyState.N*alloyState.N;
        int i;

	alloyState.Atom=(short*)malloc(sizeof(short)*alloyState.N_3);
	alloyState.Atomo=(short*)malloc(sizeof(short)*alloyState.N_3);
	ptState.att=ptState.acc=0;
	init_genrand(mpiState.myrank*rand());
 
	initialize();
#ifdef RESTART
        if(ptState.Restart == 1)
	        read_state();
#endif
}

void swap(bool even){
	int i,j,k, itag=123;
	double delta, Ei, Ej, Ti, Tj;
	bool flag;
	MPI_Status istatus;
	unsigned mem_u = alloyState.N_3*sizeof(short);
	short* buf_ui = (short*) malloc(mem_u);
	short* buf_uj = (short*) malloc(mem_u);

	if(even){
		if(mpiState.myrank%2 == 0){
			Ei = alloyState.currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, mpiState.myrank+1, itag,
									&Ej,1,MPI_DOUBLE, mpiState.myrank+1,itag,
									MPI_COMM_WORLD,&istatus);
			delta = (Ei-Ej)*E_scale*(1/(T_scale*ptState.T[mpiState.myrank+1])-1/(T_scale*ptState.T[mpiState.myrank]));
			if(delta <=0 || genrand_real2() < exp(-delta))
				flag = true;
			else
				flag = false;
			MPI_Send(&flag, 1, MPI_LOGICAL,mpiState.myrank+1,itag, MPI_COMM_WORLD);
			ptState.att++;
			if(flag){
				memcpy(&buf_ui[0], alloyState.Atom, mem_u);
				MPI_Sendrecv(buf_ui, alloyState.N_3, MPI_SHORT, mpiState.myrank+1, itag,
										buf_uj,alloyState.N_3,MPI_SHORT, mpiState.myrank+1,itag,
										MPI_COMM_WORLD,&istatus);

				memcpy(alloyState.Atom, &buf_uj[0], mem_u);

				alloyState.currEtot = Ej;
				ptState.acc ++;
			}
		}else{
			Ei = alloyState.currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, mpiState.myrank-1, itag,
									&Ej,1,MPI_DOUBLE, mpiState.myrank-1,itag,
									MPI_COMM_WORLD,&istatus);
			MPI_Recv(&flag, 1, MPI_LOGICAL,mpiState.myrank-1,itag, MPI_COMM_WORLD, &istatus);	
			if(flag){
				memcpy(&buf_ui[0], alloyState.Atom, mem_u);
				MPI_Sendrecv(buf_ui, alloyState.N_3, MPI_SHORT, mpiState.myrank-1, itag,
										buf_uj,alloyState.N_3,MPI_SHORT, mpiState.myrank-1,itag,
										MPI_COMM_WORLD,&istatus);

				memcpy(alloyState.Atom, &buf_uj[0], mem_u);

				alloyState.currEtot = Ej;

			}
		}	

	}else{
		if(mpiState.myrank%2 != 0 && mpiState.myrank!= mpiState.nprocs-1){
			Ei = alloyState.currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, mpiState.myrank+1, itag,
									&Ej,1,MPI_DOUBLE, mpiState.myrank+1,itag,
									MPI_COMM_WORLD,&istatus);
			delta = (Ei-Ej)*E_scale*(1/(T_scale*ptState.T[mpiState.myrank+1])-1/(T_scale*ptState.T[mpiState.myrank]));
			if(delta <=0 || genrand_real2() < exp(-delta))
				flag = true;
			else
				flag = false;
			MPI_Send(&flag, 1, MPI_LOGICAL,mpiState.myrank+1,itag, MPI_COMM_WORLD);
			ptState.att++;
			if(flag){
				memcpy(&buf_ui[0], alloyState.Atom, mem_u);
                                MPI_Sendrecv(buf_ui, alloyState.N_3, MPI_SHORT, mpiState.myrank+1, itag,
                                                                                  buf_uj,alloyState.N_3,MPI_SHORT, mpiState.myrank+1,itag,
                                                                                  MPI_COMM_WORLD,&istatus);

                                memcpy(alloyState.Atom, &buf_uj[0], mem_u);

				alloyState.currEtot = Ej;

				ptState.acc++;
			}
		}else if(mpiState.myrank != 0 && mpiState.myrank != mpiState.nprocs-1){
			Ei = alloyState.currEtot;
			MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, mpiState.myrank-1, itag,
									&Ej,1,MPI_DOUBLE, mpiState.myrank-1,itag,
									MPI_COMM_WORLD,&istatus);
			MPI_Recv(&flag, 1, MPI_LOGICAL,mpiState.myrank-1,itag, MPI_COMM_WORLD, &istatus);	
			if(flag){
				memcpy(&buf_ui[0], alloyState.Atom, mem_u);
                                MPI_Sendrecv(buf_ui, alloyState.N_3, MPI_SHORT, mpiState.myrank-1, itag,
                                                                                 buf_uj,alloyState.N_3,MPI_SHORT, mpiState.myrank-1,itag,
                                                                                 MPI_COMM_WORLD,&istatus);

                                memcpy(alloyState.Atom, &buf_uj[0], mem_u);

				alloyState.currEtot = Ej;


			}
		}		

	}


	ini_W();

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

      R=exp(-(Ef-Ei)*E_scale/(T_scale*ptState.pT));

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
	for(cnt=0; cnt<alloyState.N*alloyState.N*alloyState.N; cnt++){
		BondSwap(mode);
	}
}

void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int irun, SamplingMode mode){
	int mcs,i,j,k,n, t;
	FILE *ofp, *ofp1, *ofp2, *ofp_time;
	double avgE, avgE2, M, avgM, avgM2, avgM4;
	double c, x, bc;
	double* Et = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Mt = (double*)malloc(sizeof(double)*(SAMPS+DROPI));
	double* Ea = (double*)malloc(sizeof(double)*nT);
	double* Ma = (double*)malloc(sizeof(double)*nT);
	double* ca = (double*)malloc(sizeof(double)*nT);
	double* xa = (double*)malloc(sizeof(double)*nT);
	double* bca = (double*)malloc(sizeof(double)*nT);
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
	if(mpiState.myrank ==0)
		t1 = time(NULL);

//#ifndef RESTART
        if(ptState.Restart == 0){
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
        if(ptState.Restart == 1){
	        sprintf(s, "mc%d.input", mpiState.myrank);
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
				if(mpiState.myrank == 0){
                			ofp_time= fopen("infer.dat","w");
                			//t3 = time(NULL);
					t3 = std::chrono::high_resolution_clock::now();
        			}
			}*/
			mchybrid(mode);
			/*if(mcs == 1){
				if(mpiState.myrank == 0){
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

		avgE +=alloyState.currEtot;
		avgE2 +=alloyState.currEtot*alloyState.currEtot;
		M = L1();
		avgM += M;
		avgM2 += M*M;
		avgM4 += M*M;
/////////////////////////////////////////////////////////
		if(mcs % CHPT_STEPS == 0){ //checkpoint

			printf("mpiState.myrank =%d, currEtot = %g, Eng = %g\n", mpiState.myrank, alloyState.currEtot, Etot());
			if(mpiState.myrank==0){
				t2 = time(NULL);
				if(1.0*(t2-t1)/3600 >  0.9*TIMER)
					stop = true;
			}
			MPI_Bcast(&stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
			if(stop){

				write_state();
				sprintf(s, "mc%d.input", mpiState.myrank);
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
	if(mpiState.myrank == 0)
		write_pos(); 
	avgE/=1.0*SAMPS;
        avgE2/=1.0*SAMPS;
	avgM/=1.0*SAMPS;
        avgM2/=1.0*SAMPS;
        avgM4/=1.0*SAMPS;
	c = (avgE2 - avgE*avgE)*E_scale*E_scale/(T_scale*T_scale*ptState.pT*ptState.pT) ;
	x = (avgM2 - avgM*avgM)*alloyState.N_3/ptState.pT/T_scale;
        bc = 1 - avgM4/avgM2/avgM2/3;

#ifdef Time_Series
	if(mpiState.myrank < mpiState.nprocs){
		sprintf(s, "EMt%g.dat", ptState.pT);
		ofp1 = fopen(s,"w");
		for(i=0;i<(SAMPS+DROPI);i++)
				fprintf(ofp1,"%d\t%g\t%g\n",i,Et[i],Mt[i]);
		fclose(ofp1);
	}
#endif

#ifdef HIST
	if(mpiState.myrank < mpiState.nprocs){
		std::map<double,int>::iterator m_it;

		sprintf(s, "histE%g.dat", ptState.pT);
		ofp1 = fopen(s,"w");
		for(m_it = histE.begin(); m_it!= histE.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);

		sprintf(s, "histM%g.dat", ptState.pT);
		ofp1 = fopen(s,"w");
		for(m_it = histM.begin(); m_it!= histM.end(); m_it++)
				fprintf(ofp1,"%g\t%d\n",(*m_it).first,(*m_it).second );
		fclose(ofp1);
	}
#endif

	MPI_Gather(&avgE, 1, MPI_DOUBLE, Ea, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&avgM, 1, MPI_DOUBLE, Ma, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&c, 1, MPI_DOUBLE, ca, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&x, 1, MPI_DOUBLE, xa, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&bc, 1, MPI_DOUBLE, bca, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Gather(&ptState.att, 1, MPI_INT, attempts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	MPI_Gather(&ptState.acc, 1, MPI_INT, accepts, 1, MPI_INT, 0 , MPI_COMM_WORLD);

	if(mpiState.myrank == 0){
		sprintf(s, "stat%d.dat", irun);
		ofp = fopen(s,"w");
		for(n=0;n<nT;n++){
		  fprintf(ofp,"%g\t%g\t%g\t",
			  ptState.T[n],Ea[n]*alloyState.invN,ca[n]*alloyState.invN
			 /* Ma[n*5+3], xa[n*5+3], bca[n*5+3]*/);
		  fprintf(ofp, "%g\t%g\t%g\t",Ma[n+i],xa[n*+i],bca[n+i]);
		  fprintf(ofp,"\n");
		}	
		fclose(ofp);
		  
		sprintf(s, "misc%d.dat", irun);
		ofp2 = fopen(s,"w");
		fprintf(ofp2,"swap prob:\n");
		for(n=0;n<nT-1;n++){
			fprintf(ofp2,"%g\t%g\n",ptState.T[n], 1.0*accepts[n]/attempts[n]);		
		}
	}
	MPI_Gather(&alloyState.attd, 1, MPI_INT, attempts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	MPI_Gather(&alloyState.accd, 1, MPI_INT, accepts, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	if(mpiState.myrank == 0){
		fprintf(ofp2,"Rot prob:\n");
		for(n=0;n<nT;n++){
			fprintf(ofp2,"%g\t%g\n",ptState.T[n], 1.0*accepts[n]/attempts[n]);		
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
	free(ptState.T);
	for(t=0; t<NE; t++)
		free(alloyState.inputPos[t]);
}


void write_state(){

	char s[512];
	sprintf(s, "state%d.input", mpiState.myrank);
	FILE* fp = fopen(s, "wb");

	fwrite(alloyState.Atom, sizeof(short), alloyState.N_3, fp);
	fwrite(alloyState.NT, sizeof(int), NE, fp);
	fwrite(&alloyState.currEtot, sizeof(double), 1, fp);
	fclose(fp);

}
void read_state(){

	char s[512];
	double eng;
	sprintf(s, "state%d.input", mpiState.myrank);
	FILE* fp = fopen(s, "rb");
	if(fp == NULL){
		printf("state file no found\n");
		exit(-1);
	}
	fread(alloyState.Atom, sizeof(short), alloyState.N_3, fp);
	fread(alloyState.NT, sizeof(int), NE, fp);
	fread(&eng,sizeof(double), 1, fp);
        ini_apos();
	alloyState.currEtot = Etot();
	if(fabs(eng-alloyState.currEtot)>1e-5){
		printf("state file corrupted!\n");
		exit(-2);
	}

	fclose(fp);

}
