#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <set>

#define  DEFINE_GLOBALS
#include "model.hpp"
#include "alloy.hpp"
#include "rand.hpp"
#include "pt.hpp"
#include "wanglandau.hpp"

extern int myrank, nprocs;
void initialize(){

	//maximum rotational angle;
	D=1.0*Pi; DD = 0.2;
	attd=accd=0;

	//initialize spin configuration;
	//ini_conf(state);

	ini_coupling();
//	write_mol2(-1);
	invN=1.0/(N*N*N);
	
}

void ini_coupling(){
	int i,j,k,ii,cnt,shell,ti[NE*(NE-1)/2 + NE], tj[NE*(NE-1)/2 + NE],x,y,z;
	double r, Jr[NE*(NE-1)/2 + NE],Jrold;
	FILE* fop;
	for(i=0; i<NE; i++)
		for(j=0; j<NE; j++)
			for(k=0; k<SH; k++)
				J[i][j][k] = 0.0;

// read from input file;
	fop = fopen("coupling.input","r");
        if(fop==NULL){
		printf("coupling.input file was not opened\n");
		exit(1);
        }
	//for(i=0;i<100; i++) // 10 shell, each with 10 J;
	shell = 0; Nneighbors=0; Dist[0] = 0.0;
	while( shell < SH)
        {
		fscanf(fop, "%lg", &r);
	//	for(i=0;i<(NE*(NE-1)/2 + NE);i++){
			i=0;
			fscanf(fop, "%lg %d %d", &Jr[i], &ti[i], &tj[i]);
	//	}
		fscanf(fop, "%d %d %d\n", &x, &y, &z);
		  
		if(Nneighbors == 0){
			Jrold = Jr[0];
		}else{
			// shell by coupling 
			//if(fabs(Jrold -Jr[0]) > 1e-10){
			// shell by distance
			if(fabs(r- Dist[shell]) > 1e-6){
				shell++;
				Jrold = Jr[0];
				if(shell == SH) break;
			}
		}
		nlist[Nneighbors].x = x;
		nlist[Nneighbors].y = y;
		nlist[Nneighbors].z = z;
	//	for(i=0;i<(NE*(NE-1)/2 + NE);i++)
			i = 0; 
			J[ti[i]-1][tj[i]-1][shell]=J[tj[i]-1][ti[i]-1][shell]= 2*Jr[i];
		Dist[shell] = r;  
		NS[shell]++; // number of neighbors within each shell
		Nneighbors++; // total number of neighbors
		
        };

	fclose(fop);	

	if(myrank==0){
		for(i=0;i<SH;i++)
			fprintf(stderr,"%d\t%d\t%lf\t%lf\n",i,NS[i],Dist[i],J[0][0][i]);
	}
	inputPos = (int**)malloc(sizeof(int*)*N_3); 
	for(i=0; i<N_3; i++){
		inputPos[i] = (int*)malloc(sizeof(int)*Nneighbors);
	}


}

void shuffle(int *array, size_t n)
{
        size_t i;
        for (i = n-1; i > 0; i--) 
        {
          size_t j = rand() % (i+1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
}

void ini_apos(){
	int idx, i, j, k, x, y, z, cnt, ii, shell;
	int nn[MAX_NEIGHBORS*3], cntt;  

// generate neighbor list for neural network inputsa
	for(i=0; i<N; i++)for(j=0; j<N; j++)for(k=0; k<N; k++){
		neighbor(i,j,k,nn);
		cnt = cntt = 0; 
		idx = i*N_2+j*N+k;
		for(shell =0; shell < SH; shell++)
			for(ii=0;ii<NS[shell];ii++){
				x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
				inputPos[idx][cntt++] = x*N_2+y*N+z;
		}			
	} 
	
}
void ini_alloy(int state){

	int i,j,k,t,t2,Ni,ii,shell,x,y,z;
	int cnt,cnti,cntt[NE];
	double p;
	int * list = (int*) malloc(sizeof(int)*N*N*N);
// initialize alloy atom species with equal probabilities.
	FILE* fop = fopen("composition.input","r");
        if(fop==NULL){
		printf("composition.input file was not opened\n");
		exit(1);
        }
	cnt = 0;
	for(t=0;t<NE;t++){
		fscanf(fop, "%lg", &p); 
		Ni = (int)(N*N*N*p);
		cnti=0;
		while(cnti < Ni){
			list[cnt++] = t;
			cnti++;
		}
		NT[t] = 0;
	}
	while(cnt < N*N*N){
		t = (int)(randd1()*NE);
		list[cnt++] = t;
	}
        if(state == 0)
		shuffle(list, N*N*N);
	cnt = 0;	
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			for(k=0; k<N; k++){		
				//t = (int)(randd1()*NE);
				t = list[cnt++];
				Atom[i*N_2+j*N+k] = t;
				//cnt[t]++;
				NT[t]++;
	}
	free(list);
        fclose(fop);
	for(t = 0 ; t < NE; t++)
		fprintf(stderr, "%s:%d\n",element[t+1],NT[t]);



/*	for(t = 0 ; t < NE; t++){
			while(cnt[t] > Ni[t]){

				do{
					i = (int) (randd1()*N);
					j = (int) (randd1()*N);
					k = (int)(randd1()*N);	
				}while(Atom[i*N_2+j*N+k] != t);

				do{
					t2 = (int)(randd1()*NE);
				}while(t2 <= t);
				Atom[i*N_2+j*N+k] = t2;

				cnt[t]--;
				cnt[t2]++;
			}

			while(cnt[t] < Ni[t]){

				do{
					i = (int) (randd1()*N);
					j = (int) (randd1()*N);
					k = (int)(randd1()*N);	
				}while(Atom[i*N_2+j*N+k] <= t);
				t2 = Atom[i*N_2+j*N+k]; 				
				Atom[i*N_2+j*N+k] = t;

				cnt[t]++;
				cnt[t2]--;
			}
		}
	for(t=0;t<NE;t++)
		assert(cnt[t] == Ni[t]);*/
/////////////////////////////////////////
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
	
/*	si = i - offi;
	if(si < 0) si += N;
	if(si >= N) si -= N;
	sj = j - offj;
	if(sj < 0) sj += N;
	if(sj >= N) sj -= N;
	sk = k - offk;
	if(sk < 0) sk += N;
	if(sk >= N) sk -= N;
	nn[cnt++] = si; nn[cnt++] = sj; nn[cnt++] =sk;
*/}
void  neighbor(int i, int j, int k, int* nn){
	int ii,cnt=0;
	for(ii=0;ii<Nneighbors;ii++){
		noffset(i,j,k,nlist[ii].x,nlist[ii].y,nlist[ii].z,nn,cnt);
		cnt += 3;
	}
}

double Etot(){
	double E;
	E = 0.0;

	#pragma omp parallel for reduction(+: E) num_threads(NE)
	for(int t = 0; t < NE; t++){
		tensorflow::TensorShape data_shape({NT[t], Nneighbors*(NE-1)});
        	tensorflow::Tensor data(tensorflow::DT_UINT8, data_shape);
        	auto data_ = data.flat<std::uint8_t>().data();
                int cnt=0;
		for(int i=0; i<N_3; i++){
			if(Atom[i] == t){
				for(int j = 0; j < Nneighbors; j++){//num of neighbors  
					for(int k=0; k<NE-1; k++){//num of representation 
                                                //printf("i%d,j%d,input%d\n",i,j,inputPos[i][j]);
						data_[cnt++] = encode[Atom[inputPos[i][j]]][k]; 
                                        }
				}
			}
		}
		std::vector<tensorflow::Tensor> outputs;
                models[t].Predict(data, outputs);
		for(int i = 0; i < NT[t]; i++){	
			E += outputs[0].flat<float>().data()[i];
		}
	}
	return E; 

}
/*
double Esite(int i, int j, int k){

	int ii,cnt,ti,tj,x,y,z, shell;
	double sx, sy, sz;
	double Elocal = 0.0;
	int nn[MAX_NEIGHBORS*3]; // 10 shells contains 200 neighbors for fcc; 

	ti = Atom[i*N_2+j*N+k];

	neighbor(i,j,k,nn);
	cnt = 0; 
	for(shell =0; shell < SH; shell++){
		for(ii=0;ii<NS[shell];ii++){
			x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
			tj = Atom[x*N_2+y*N+z];
		}
	}
	
	tensorflow::TensorShape data_shape({1, 6});
        tensorflow::Tensor data(tensorflow::DT_FLOAT, data_shape);
        
        auto data_ = data.flat<float>().data();
        data_[0] = 0;
	data_[1] = 0;
  	data_[2] = 1;
  	data_[3] = 1;
  	data_[4] = 0;
  	data_[5] = 0;

	tensor_dict feed_dict = {
      		{"input_liz:0", data},
  	};
	std::vector<tensorflow::Tensor> outputs;
	TF_CHECK_OK(
		sess->Run(feed_dict, {"sequential/BiasAdd_2:0"}, {}, &outputs));
	Elocal = outputs[0].flat<float>().data()[0];

        //std::cout<< Elocal << std::endl; 
	return Elocal;
}*/
/*
void wolff(int i, int j, int k, double rx, double ry, double rz){

	double sx,sy,sz,pi,pj, delta;
	int ti,tj, cnt, ii, shell, x, y, z, nn[MAX_NEIGHBORS*3];

	cluster[i*N_2+j*N+k] = true;
	sx = S[i*N_2+j*N+k];
	sy = S[N_3+i*N_2+j*N+k];
	sz = S[2*N_3+i*N_2+j*N+k];
	pi = (sx*rx + sy*ry + sz*rz);
	ti = Atom[i*N_2+j*N+k];

	S[i*N_2+j*N+k] = sx - 2*pi*rx;
	S[N_3+i*N_2+j*N+k] = sy - 2*pi*ry;
	S[2*N_3+i*N_2+j*N+k] = sz - 2*pi*rz;

		neighbor(i,j,k,nn);
		cnt = 0; 
		for(shell =0; shell < SH; shell++){
			for(ii=0;ii<NS[shell];ii++){
				x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
				if(!cluster[x*N_2+y*N+z]){
					tj = Atom[x*N_2+y*N+z];
					pj = (S[x*N_2+y*N+z]*rx + S[N_3+x*N_2+y*N+z]*ry + S[2*N_3+x*N_2+y*N+z]*rz);
					delta = -2*J[ti][tj][shell]/(T_scale*pT)*pi*pj;
					if(delta < 0 && randd1() < (1-exp(delta)) )	
						wolff(x,y,z,rx,ry,rz);
				}
			}
		}
}
*/
void Rot(){
	int i, j, k, x, y, z, ii, jj, kk, n, it, shell;
	int ai, aj, cnt, cntE; 
	double deltaE, dE, E1, E2;
        //int nn[MAX_NEIGHBORS*3];

	i = (int) (randd1()*N_3);
	n = (int) (randd1()*NS[0]);
	//neighbor(i/N_2,(i/N)%N,i%N,nn);			
        //x = nn[n*3]; y = nn[n*3+1]; z = nn[n*3+2];
	//j = x*N_2+y*N+z;
	j = inputPos[i][n];	
	E1 = currEtot;
	ai = Atom[i];
	aj = Atom[j];
	Atom[i] = aj;
	Atom[j] = ai;	
	E2 = Etot();

	deltaE = E2 - E1;
	attd++;
        if(WangLandau(currEtot, currEtot+deltaE) == 1){// accept
        	currEtot += deltaE; accd++;
        }else{//reject
		Atom[i] = ai;
		Atom[j] = aj;
	}
	std::cout << "E: " << currEtot << std::endl;
}

void write_pos(){
  char s[512];
  FILE *ofp;
  //sprintf(s,"compos%d.dat",myrank);
  ofp=fopen("compos.dat","w");
	int i,j,k;
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			for(k=0;k<N;k++){
				fprintf(ofp, "%d\t%d\t%d\t%d\t \n",i+k, i+j, j+k,Atom[i*N_2+j*N+k]);
			}
  fclose(ofp);
}

void write_mol2(int frame)
{
  int i,j,k,cnt;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"snap_%d_%d.mol2",frame, myrank);
  ofp=fopen(s,"w");

  fprintf(ofp,"@<TRIPOS>MOLECULE\n");
  fprintf(ofp,"Eng = %g\n", currEtot/N/N/N);
  fprintf(ofp," %d %d\n",N*N*N,N*N*N);
  fprintf(ofp,"SMALL\n");
  fprintf(ofp,"NO_CHARGES\n");
  fprintf(ofp,"@<TRIPOS>ATOM\n");

  cnt = 0;
  for(i=0;i<N;i++)
	  for(j=0; j<N;j++)
		  for(k=0;k<N;k++)
  {
	fprintf(ofp,"%d %s %d %d %d \n",cnt++, element[Atom[i*N_2+j*N+k]+1], i+k, i+j, j+k );

  };

  fclose(ofp);

}


// HEA order parameter 
void O(double* M){
	int i, j, k, x, y, z, cnt, shell, t, it, ai, aj, neighbors;
        int nn[MAX_NEIGHBORS*3];
	double Occ[NE], op[NE][NE];

	for(t = 0; t < NE; t++){
		Occ[t] = 1.0*NT[t]/N_3;
		for(it = 0; it < NE; it++)
			op[t][it] = 0.0;
	}
	neighbors = 0;
	for(shell =0; shell < O_SH; shell++)
		neighbors += NS[shell];
		
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			for(k=0; k<N; k++){
				neighbor(i,j,k,nn);
				cnt = 0; ai = Atom[i*N_2+j*N+k];
				for(shell =0; shell < O_SH; shell++)for(it=0;it<NS[shell];it++){
					x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
					aj = Atom[x*N_2+y*N+z];
					op[ai][aj] += 1.0/neighbors; 
				}
				//op[ai] += fabs( 1.0*M/((cnt+1)/3) - Occ[ai]);
	}
	//op[NE] = 0.0;
	for(t = 0; t < NE; t++){
		for(it = 0; it < NE; it++){
			op[t][it] /= N_3; //NT[t];
			op[t][it] = 1 - op[t][it]/(Occ[t]*Occ[it]);
		}
		//op[t] /= (NT[t]*2*Occ[t]*(1-Occ[t]));
		//op[NE] += op[t];
	}
	M[0] = op[0][2]; //Mo-Ta
	M[1] = op[1][2]; //Nb-Ta
	M[2] = op[1][3]; //Nb-W
	M[3] = op[2][3]; //Ta-W
	M[NE] = 0.0;
	for(t = 0; t < NE; t++){
		//M[t] = op[0][t];  
		M[NE] += fabs(M[t]);
	}	
	M[NE] /= NE; 	
	//op[NE] /= NE;			
	
}



