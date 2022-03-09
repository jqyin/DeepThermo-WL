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
#include "client.hpp"
#include "wanglandau.hpp"

extern std::string encoder_name, decoder_name;
extern int myrank, nprocs;
void initialize(){
	//maximum rotational angle;
	D=1.0*Pi; DD = 0.2;
	attd=accd=0;

	ini_coupling();
	invN=1.0/(N*N*N);
		
	SHIFT = N-1; 
	VAE_D = 2*(N-1)+ SHIFT + 1; 
	PAD = int((ceil(1.0*VAE_D/16)*16 - VAE_D)/2);
        VAE_D += PAD*2; 
	//z[0]=z[1]=z[2]=0.0;
	
	inputConfig = (float*) malloc(sizeof(float)*VAE_D*VAE_D*VAE_D*NE);
	memset(inputConfig,0,sizeof(float)*VAE_D*VAE_D*VAE_D*NE);
}

void ini_coupling(){

        int i,j,k,t,ii,cnt,shell,ti[NE*(NE-1)/2 + NE], tj[NE*(NE-1)/2 + NE],x,y,z;
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

        shell = 0; Nneighbors=0;
        for(i=0;i<SH;i++){
                NS[i] = 0;
                Dist[i] = 0.0;
        }
        while( shell < SH)
        {
                fscanf(fop, "%lg", &r);
                for(i=0;i<(NE*(NE-1)/2 + NE);i++){
                        fscanf(fop, "%lg %d %d", &Jr[i], &ti[i], &tj[i]);
                }
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
                for(i=0;i<(NE*(NE-1)/2 + NE);i++)
                        J[ti[i]-1][tj[i]-1][shell]=J[tj[i]-1][ti[i]-1][shell]= Jr[i];
                Dist[shell] = r;
                NS[shell]++; // number of neighbors within each shell
                Nneighbors++; // total number of neighbors
        }

        fclose(fop);

        if(myrank==0){
                for(i=0;i<SH;i++)
                        fprintf(stderr,"%d\t%d\t%lf\t%lf\n",i,NS[i],Dist[i],J[0][1][i]);
        }
        inputPos = (int**)malloc(sizeof(int*)*N_3);
        for(i=0; i<N_3; i++){
                inputPos[i] = (int*)malloc(sizeof(int)*Nneighbors);
        }
	for(t=0;t<NE;t++)
		elist[t] = (int*)malloc(sizeof(int)*N_3);
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

void ini_W(){
        int i, j, k, idx, cnt, shell, ii, ai, aj;
        double Wsum = 0.0;
        memset(W, 0, sizeof(int)*NE*NE*SH);
        // init W matrix
        for(i=0; i<N; i++)for(j=0; j<N; j++)for(k=0; k<N; k++){
                idx = i*N_2+j*N+k;
                ai = Atom[idx];
                cnt = 0;
                for(shell =0; shell < SH; shell++){
                        for(ii=0;ii<NS[shell];ii++){
                                aj = Atom[inputPos[idx][cnt++]];
                                if(ai <= aj)
                                        W[ai][aj][shell]++;
                                else
                                        W[aj][ai][shell]++;
                        }
                }
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
	ini_W();	
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

void updateWsite(int ai, int aj, int pi, int pj){
        int shell, j, id, a;
        int cnt=0;
        for(shell = 0; shell < SH; shell++){
                for(j=0; j<NS[shell]; j++){
                        id = inputPos[pi][cnt++];
                        if(id != pj){
                                a = Atom[id];
                                if(ai <= a)
                                        W[ai][a][shell] -= 2;
                                else
                                        W[a][ai][shell] -= 2;
                                if(aj <= a)
                                        W[aj][a][shell] += 2;
                                else
                                        W[a][aj][shell] += 2;
                        }
                }
        }

}

void updateW(int ai, int aj, int idxi, int idxj){
        //updateWsite(Atom[idxi],idxi,idxj,false);
        //updateWsite(Atom[idxj],idxi,idxj,true);
        //updateWsite(Atom[idxj],idxj,idxi,false);
        //updateWsite(Atom[idxi],idxj,idxi,true);
        updateWsite(ai, aj, idxi, idxj);
        updateWsite(aj, ai, idxj, idxi);
}


double Etot(){
	double E;
	E = 0.0;

#ifdef DL_MODEL

#ifdef Backend_TF
	tensorflow::TensorShape data_shape({1, 1, SH*NE*(NE-1)/2});
        tensorflow::Tensor data(tensorflow::DT_HALF, data_shape);
        auto data_ = data.flat<Eigen::half>().data();
        //tensorflow::Tensor data(tensorflow::DT_FLOAT, data_shape);
        //auto data_ = data.flat<float>().data();
        int cnt=0;
       	for(int shell = 0; shell < SH; shell++){
                for(int i =0; i < NE-1; i++){
                        for(int j =i+1; j < NE; j++){
				//data_[cnt++] = (invN*W[i][j][shell]/NS[shell]);
				data_[cnt++] = (Eigen::half)(invN*W[i][j][shell]/NS[shell]);
			}
		}
	}

	std::vector<tensorflow::Tensor> outputs;
        model.Predict(data, outputs);
	E = outputs[0].flat<float>().data()[0] + mlp_intercept;
        E = E/E_scale; 
#else
        float data[SH*NE*(NE-1)/2];
	std::string data_key, eng_key;
	std::string rankID = std::to_string(myrank);
        int cnt=0;
       	for(int shell = 0; shell < SH; shell++){
                for(int i =0; i < NE-1; i++){
                        for(int j =i+1; j < NE; j++){
				data[cnt++] = (invN*W[i][j][shell]/NS[shell]);
			}
		}
	}

						
	data_key = "data_"+rankID;
	eng_key = "eng_"+rankID;

	//if (MC == 5&& myrank ==0)
	//	t1 = std::chrono::high_resolution_clock::now();
	(*SRclient).put_tensor(data_key, data, {1,1,SH*NE*(NE-1)/2}, SmartRedis::TensorType::flt,
				SmartRedis::MemoryLayout::contiguous);
	//if (MC == 5&& myrank ==0)
	//	t2 = std::chrono::high_resolution_clock::now();
        //std::cout << model_name << std::endl;
	(*SRclient).run_model(model_name+rankID, {data_key}, {eng_key});
	//if (MC == 5 && myrank ==0)
	//	t3 = std::chrono::high_resolution_clock::now();
	//std::vector<float> eng;
        float eng;
	(*SRclient).unpack_tensor(eng_key, &eng, {1}, SmartRedis::TensorType::flt,
                                SmartRedis::MemoryLayout::contiguous);
	//if (MC == 5 && myrank ==0)
		//t4 = std::chrono::high_resolution_clock::now();
	E = eng + mlp_intercept;
        E = E/E_scale;

	/*if (MC == 5 && myrank == 0){
		std::chrono::duration<double, std::milli> ms_double = t2 - t1;
		std::cout << "put_tensor: " << ms_double.count() << std::endl;
		ms_double = t3-t2; 
		std::cout << "run_model: " << ms_double.count() << std::endl;
		ms_double = t4-t3; 
		std::cout << "unpack_tensor: " << ms_double.count() << std::endl;
	}*/
#endif


#else // regression model 
       	for(int shell = 0; shell < SH; shell++){
                for(int i =0; i < NE-1; i++){
                        for(int j =i+1; j < NE; j++){
                                E += (1.0*W[i][j][shell]/NS[shell]/N_3)*J[i][j][shell];
			}
		}
	}
	E += reglin_intercept;
#endif
	return E*N_3; 

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
void Rot(int mode){
	int i, j, k, x, y, z, ii, jj, kk, n, it, shell;
	int ai, aj, cnt, cntE; 
	int Wo[NE][NE][SH];
	double deltaE, dE, E1, E2;
        //int nn[MAX_NEIGHBORS*3];

        i = (int) (randd1()*N_3);
        n = (int) (randd1()*Nneighbors);
        j = inputPos[i][n];
    	if(Atom[i] != Atom[j]){
        	E1 = currEtot;
        	ai = Atom[i];
        	aj = Atom[j];
        	memcpy(Wo, W, sizeof(int)*NE*NE*SH);
        	updateW(ai, aj, i, j);
        	Atom[i] = aj;
        	Atom[j] = ai;
        	E2 = Etot();

        	deltaE = E2 - E1;
        	//attd++;
		if(mode == 0){//metropolis
        		if(Metropolis(currEtot, currEtot+deltaE) == 1){// accept
        			currEtot += deltaE; //accd++;
        		}else{//reject
                		Atom[i] = ai;
                		Atom[j] = aj;
                		memcpy(W, Wo, sizeof(int)*NE*NE*SH);
        		}
		}
		else{//wanglandau
        		if(WangLandau(currEtot, currEtot+deltaE, mode) == 1){// accept
                		currEtot += deltaE; //accd++;
        		}else{//reject
                		Atom[i] = ai;
                		Atom[j] = aj;
                		memcpy(W, Wo, sizeof(int)*NE*NE*SH);
        		}

			//std::cout << "E: " << currEtot << std::endl;
		}
    	}

}

void encode(float* z){
	int i,j,k,t; 
	std::string z_key, config_key;
	std::string rankID = std::to_string(myrank);
	memset(inputConfig,0,sizeof(float)*VAE_D*VAE_D*VAE_D*NE);
  	for(i=0;i<N;i++)for(j=0; j<N;j++)for(k=0;k<N;k++){
		int idx = (j-i+k+SHIFT+PAD)*VAE_D*VAE_D*NE + (k-j+i+SHIFT+PAD)*VAE_D*NE + (j+i-k+SHIFT+PAD)*NE;
		for(t =0; t <NE; t++)
			if(Atom[i*N_2+j*N+k]==t){
				inputConfig[idx+t] = 1;
				break;
			}
	}
	z_key = "output_z_"+rankID;
	config_key = "input_config_"+rankID;
        (*SRclient).put_tensor(config_key, inputConfig, {1, VAE_D, VAE_D, VAE_D, NE}, SmartRedis::TensorType::flt,
                                SmartRedis::MemoryLayout::contiguous);
	(*SRclient).run_model(encoder_name+rankID, {config_key}, {z_key});
	(*SRclient).unpack_tensor(z_key, z, {3}, SmartRedis::TensorType::flt,SmartRedis::MemoryLayout::contiguous);
}

void walk(float* npos){
	// shhere
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

	npos[0] += Z_R*v1;
	npos[1] += Z_R*v2;
	npos[2] += Z_R*v3; 
	
	// line
	/*a = randd1()*2 - 1;
	npos[0] += Z_R*a;
	npos[1] += Z_R*a;
	npos[2] += Z_R*a;*/
	

}

int arg_max(std::vector<float> const& vec){
	return static_cast<int>(std::distance(vec.begin(), std::max_element(vec.begin(), vec.end())));
}

void decode(float* z){
  	int i, j, k, t, tt, id, tid;	
	int sum[NE] = {0};
	int sumo[NE] = {0};
	float output[VAE_D*VAE_D*VAE_D*NE] = {0};
	std::vector<float> vec(NE); 	
	std::string z_key, config_key;
	std::string rankID = std::to_string(myrank);
	
	z_key = "input_z_"+rankID;
	config_key = "output_config_"+rankID;
        (*SRclient).put_tensor(z_key, z, {1,3}, SmartRedis::TensorType::flt,
                                SmartRedis::MemoryLayout::contiguous);
	(*SRclient).run_model(decoder_name+rankID, {z_key}, {config_key});
	(*SRclient).unpack_tensor(config_key, &output, {VAE_D*VAE_D*VAE_D*NE}, SmartRedis::TensorType::flt,SmartRedis::MemoryLayout::contiguous);

	for(t=0;t<NE;t++)
		memset(elist[t],-1,sizeof(int)*N_3);
  	for(i=0;i<N;i++)for(j=0; j<N;j++)for(k=0;k<N;k++){
		int idx = (j-i+k+SHIFT+PAD)*VAE_D*VAE_D*NE + (k-j+i+SHIFT+PAD)*VAE_D*NE + (j+i-k+SHIFT+PAD)*NE;
		for(t=0; t<NE; t++)
			vec[t] = output[idx+t];
		t = arg_max(vec);
        	Atom[i*N_2+j*N+k] = t;
		elist[t][sum[t]++] = i*N_2+j*N+k;
	}
	for(t=0;t<NE;t++)
		sumo[t] = sum[t];
	// fix concentration, better to take prob into account 
	for(t=0;t<NE;t++){
		while(sum[t] > NT[t]){
			id = int(randd1()*sumo[t]);
			if(elist[t][id] != -1){
				for(tt=0;tt<NE;tt++)
					if(tt!=t && sum[tt] < NT[tt]){
						tid = tt;
						break;
					}
				Atom[elist[t][id]] = tid;
				elist[t][id] = -1;
				sum[t]--;
				sum[tid]++;
			} 	
		}
	}
	for(t=0;t<NE;t++)
		assert(sum[t] == NT[t]);
	ini_W();
	 
}

void vae_update(int mode){
	float z[3] = {0};
	double E1, E2, deltaE;
	int Wo[NE][NE][SH];
        memcpy(Atomo, Atom, sizeof(short)*N_3);
        memcpy(Wo, W, sizeof(int)*NE*NE*SH);
	// encode the current config to latent space;
	encode(z);
        /*if(myrank == 0 && TotalSweeps%100 == 0){
		printf("step %d before: z=(%f,%f,%f), E=%f\n", TotalSweeps, z[0], z[1], z[2], currEtot);
		//write_xyz(TotalSweeps);
	}*/
	// random walk in latent space; 
	walk(z);
	// decode the data point to real space; 
        decode(z);
	// measure the energy of new configuration;
	E1 = currEtot;
	E2 = Etot(); 
        /*if(myrank == 0 && TotalSweeps%100==0){
		printf("step %d after: z=(%f,%f,%f), E=%f\n", TotalSweeps, z[0], z[1], z[2], E2);
		//write_xyz(TotalSweeps+1);
	}*/
	// update according to WL 	
	deltaE = E2 - E1;
        attd++;
	if(mode == 0){//metropolis
        	if(Metropolis(currEtot, currEtot+deltaE) == 1){// accept
        		currEtot += deltaE; accd++;
        	}else{//reject
        		memcpy(Atom, Atomo, sizeof(short)*N_3);
        		memcpy(W, Wo, sizeof(int)*NE*NE*SH);
		}
	}else{//wanglandau
        	if(WangLandau(currEtot, currEtot+deltaE, mode) == 1){// accept
               		currEtot += deltaE; accd++;
        	}else{//reject
        		memcpy(Atom, Atomo, sizeof(short)*N_3);
        		memcpy(W, Wo, sizeof(int)*NE*NE*SH);
        	}
	}
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

void write_xyz(int frame)
{
  int i,j,k,cnt, shell;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"snap_%d_%d.xyz",frame, myrank);
  ofp=fopen(s,"ab+");
  fprintf(ofp, "%d\n", N_3); 
  fprintf(ofp,"Eng = %.8lg  MC_step = %d\n", currEtot, frame);

  //sum = 0.0;
  for(shell = 0; shell < SH; shell++)
        for(i =0; i < NE-1; i++)
                for(j =i+1; j < NE; j++){
                        fprintf(ofp, "%.8lg\t",1.0*W[i][j][shell]/NS[shell]/N_3);
                        //sum += 1.0*W[i][j][shell]/NS[shell]/N_3;
                }
  fprintf(ofp, "\n");


  for(i=0;i<N;i++)
          for(j=0; j<N;j++)
                  for(k=0;k<N;k++)
                        fprintf(ofp,"%s %d %d %d \n",element[Atom[i*N_2+j*N+k]+1], j-i+k, k-j+i, j+i-k);
  fclose(ofp);

}

void thermoqs()
{
        int i;
        double T,U,Z,C,F,S;
        double centerE,lambda,Bw,Nfree;

        FILE *therm_op;
        char s1[512];

        //Opening File containing all thermodynamic quantities (in following order)
        //      T       U       Cv      freeE  Entropy  Rgyr2  EEdist
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
                U = 0.0;        // Initializing the average energy <E>
                Z = 0.0;        // Initializing the Partition Function
                C = 0.0;        // Initialize the Specific Heat
                F = 0.0;        // Initialize the free energy
                S = 0.0;        // Initialize the entropy
                Bw = 0.0;       //Boltzmann weight
                lambda = -1.0e300;      //Normalization shift (max value of DOS considering T)
                centerE = 0.0;  // Taking the center of the energy bin
               //Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
                for(i=0;i<D1BINS;i++)
                {
                        centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;

                        if( (lambda < ((wllng[i]) - 1.0*(centerE*E_scale)/(T*T_scale))) )// && (wllng[i] > 0.0) )
                                lambda = ((wllng[i]) - 1.0*(centerE*E_scale)/(T*T_scale));

                };

                //Central Loop for calculating thermodynamic properties from the DOS
                for(i=0;i<D1BINS;i++)
                {
                        //Taking the center of the bin
                        centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;

                        //Boltzmann Factor
                        Bw =  exp( wllng[i]  - (centerE*E_scale)/(T*T_scale) - lambda );

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
                        Nfree = -(T*T_scale)/E_scale*( lambda + log(Z) ) - U;
                };
                //Specific Heat
                C = ( (C/Z) - (U*U) )*E_scale*E_scale / (T*T*T_scale*T_scale);
                //Free Energy
                F = -(T*T_scale)/E_scale*( lambda + log(Z) ) - Nfree;
                //Entropy
                S = (U - F)*E_scale/(T*T_scale);

                fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\n",T,U*invN,C*invN,F*invN,S*invN);
        };

        fclose(therm_op);
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



