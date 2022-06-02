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

#ifndef TF_BACKEND
#include "client.hpp"
#endif

extern std::string encoder_name, decoder_name;

void initialize(){
	alloyState.attd=alloyState.accd=0;

	ini_coupling();
	alloyState.invN=1.0/(alloyState.N*alloyState.N*alloyState.N);
		
	alloyState.SHIFT = alloyState.N-1; 
	alloyState.VAE_D = 2*(alloyState.N-1)+ alloyState.SHIFT + 1; 
	alloyState.PAD = int((ceil(1.0*alloyState.VAE_D/16)*16 - alloyState.VAE_D)/2);
        alloyState.VAE_D += alloyState.PAD*2; 
	
	alloyState.inputConfig = (float*) malloc(sizeof(float)*alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE);
	memset(alloyState.inputConfig,0,sizeof(float)*alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE);
}

void ini_coupling(){

        int i,j,k,t,ii,cnt,shell,ti[NE*(NE-1)/2 + NE], tj[NE*(NE-1)/2 + NE],x,y,z;
        double r, Jr[NE*(NE-1)/2 + NE],Jrold;
        FILE* fop;

        for(i=0; i<NE; i++)
                for(j=0; j<NE; j++)
                        for(k=0; k<SH; k++)
                                alloyState.J[i][j][k] = 0.0;

// read from input file;
        fop = fopen("coupling.input","r");
        if(fop==NULL){
                printf("coupling.input file was not opened\n");
                exit(1);
        }

        shell = 0; alloyState.Nneighbors=0;
        for(i=0;i<SH;i++){
                alloyState.NS[i] = 0;
                alloyState.Dist[i] = 0.0;
        }
        while( shell < SH)
        {
                fscanf(fop, "%lg", &r);
                for(i=0;i<(NE*(NE-1)/2 + NE);i++){
                        fscanf(fop, "%lg %d %d", &Jr[i], &ti[i], &tj[i]);
                }
                fscanf(fop, "%d %d %d\n", &x, &y, &z);
                if(alloyState.Nneighbors == 0){
                        Jrold = Jr[0];
                }else{
                        // shell by distance
                        if(fabs(r- alloyState.Dist[shell]) > 1e-6){
                                shell++;
                                Jrold = Jr[0];
                                if(shell == SH) break;
                        }
                }
                alloyState.nlist[alloyState.Nneighbors].x = x;
                alloyState.nlist[alloyState.Nneighbors].y = y;
                alloyState.nlist[alloyState.Nneighbors].z = z;
                for(i=0;i<(NE*(NE-1)/2 + NE);i++)
                        alloyState.J[ti[i]-1][tj[i]-1][shell]=alloyState.J[tj[i]-1][ti[i]-1][shell]= Jr[i];
                alloyState.Dist[shell] = r;
                alloyState.NS[shell]++; // number of neighbors within each shell
                alloyState.Nneighbors++; // total number of neighbors
        }

        fclose(fop);

        if(mpiState.myrank==0){
                for(i=0;i<SH;i++)
                        fprintf(stderr,"%d\t%d\t%lf\t%lf\n",i,alloyState.NS[i],alloyState.Dist[i],alloyState.J[0][1][i]);
        }
        alloyState.inputPos = (int**)malloc(sizeof(int*)*alloyState.N_3);
        for(i=0; i<alloyState.N_3; i++){
                alloyState.inputPos[i] = (int*)malloc(sizeof(int)*alloyState.Nneighbors);
        }
	for(t=0;t<NE;t++)
		alloyState.elist[t] = (int*)malloc(sizeof(int)*alloyState.N_3);
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
        memset(alloyState.W, 0, sizeof(int)*NE*NE*SH);
        // init W matrix
        for(i=0; i<alloyState.N; i++)for(j=0; j<alloyState.N; j++)for(k=0; k<alloyState.N; k++){
                idx = i*alloyState.N_2+j*alloyState.N+k;
                ai = alloyState.Atom[idx];
                cnt = 0;
                for(shell =0; shell < SH; shell++){
                        for(ii=0;ii<alloyState.NS[shell];ii++){
                                aj = alloyState.Atom[alloyState.inputPos[idx][cnt++]];
                                if(ai <= aj)
                                        alloyState.W[ai][aj][shell]++;
                                else
                                        alloyState.W[aj][ai][shell]++;
                        }
                }
        }
}

void ini_apos(){
	int idx, i, j, k, x, y, z, cnt, ii, shell;
	int nn[MAX_NEIGHBORS*3], cntt;  

// generate neighbor list for neural network inputs
	for(i=0; i<alloyState.N; i++)for(j=0; j<alloyState.N; j++)for(k=0; k<alloyState.N; k++){
		neighbor(i,j,k,nn);
		cnt = cntt = 0; 
		idx = i*alloyState.N_2+j*alloyState.N+k;
		for(shell =0; shell < SH; shell++)
			for(ii=0;ii<alloyState.NS[shell];ii++){
				x = nn[cnt++]; y = nn[cnt++]; z = nn[cnt++];
				alloyState.inputPos[idx][cntt++] = x*alloyState.N_2+y*alloyState.N+z;
		}			
	} 
	ini_W();	
}
void ini_alloy(int state){

	int i,j,k,t,t2,Ni,ii,shell,x,y,z;
	int cnt,cnti,cntt[NE];
	double p;
	int * list = (int*) malloc(sizeof(int)*alloyState.N*alloyState.N*alloyState.N);
// initialize alloy atom species with equal probabilities.
	FILE* fop = fopen("composition.input","r");
        if(fop==NULL){
		printf("composition.input file was not opened\n");
		exit(1);
        }
	cnt = 0;
	for(t=0;t<NE;t++){
		fscanf(fop, "%lg", &p); 
		Ni = (int)(alloyState.N*alloyState.N*alloyState.N*p);
		cnti=0;
		while(cnti < Ni){
			list[cnt++] = t;
			cnti++;
		}
		alloyState.NT[t] = 0;
	}
	while(cnt < alloyState.N*alloyState.N*alloyState.N){
		t = (int)(randd1()*NE);
		list[cnt++] = t;
	}
        if(state == 0)
		shuffle(list, alloyState.N*alloyState.N*alloyState.N);
	cnt = 0;	
	for(i=0; i<alloyState.N; i++)
		for(j=0; j<alloyState.N; j++)
			for(k=0; k<alloyState.N; k++){		
				t = list[cnt++];
				alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k] = t;
				alloyState.NT[t]++;
	}
	free(list);
        fclose(fop);
	for(t = 0 ; t < NE; t++)
		fprintf(stderr, "%s:%d\n",alloyState.element[t+1],alloyState.NT[t]);

}

inline void noffset(int i, int j, int k, int offi, int offj, int offk, int*nn, int cnt){
	int si, sj,sk;
	si = i + offi;
	if(si < 0) si += alloyState.N;
	if(si >= alloyState.N) si -= alloyState.N;
	sj = j + offj;
	if(sj < 0) sj += alloyState.N;
	if(sj >= alloyState.N) sj -= alloyState.N;
	sk = k + offk;
	if(sk < 0) sk += alloyState.N;
	if(sk >= alloyState.N) sk -= alloyState.N;
	nn[cnt++] = si; nn[cnt++] = sj; nn[cnt++] =sk;
	
}
void  neighbor(int i, int j, int k, int* nn){
	int ii,cnt=0;
	for(ii=0;ii<alloyState.Nneighbors;ii++){
		noffset(i,j,k,alloyState.nlist[ii].x,alloyState.nlist[ii].y,alloyState.nlist[ii].z,nn,cnt);
		cnt += 3;
	}
}

void updateWsite(int ai, int aj, int pi, int pj){
        int shell, j, id, a;
        int cnt=0;
        for(shell = 0; shell < SH; shell++){
                for(j=0; j<alloyState.NS[shell]; j++){
                        id = alloyState.inputPos[pi][cnt++];
                        if(id != pj){
                                a = alloyState.Atom[id];
                                if(ai <= a)
                                        alloyState.W[ai][a][shell] -= 2;
                                else
                                        alloyState.W[a][ai][shell] -= 2;
                                if(aj <= a)
                                        alloyState.W[aj][a][shell] += 2;
                                else
                                        alloyState.W[a][aj][shell] += 2;
                        }
                }
        }

}

void updateW(int ai, int aj, int idxi, int idxj){
        updateWsite(ai, aj, idxi, idxj);
        updateWsite(aj, ai, idxj, idxi);
}


double Etot(){
	double E;
	E = 0.0;

       	for(int shell = 0; shell < SH; shell++){
                for(int i =0; i < NE-1; i++){
                        for(int j =i+1; j < NE; j++){
                                E += (1.0*alloyState.W[i][j][shell]/alloyState.NS[shell]/alloyState.N_3)*alloyState.J[i][j][shell];
			}
		}
	}
	E += alloyState.reglin_intercept;
	return E*alloyState.N_3; 

}

void BondSwap(SamplingMode mode){
	int i, j, k, x, y, z, ii, jj, kk, n, it, shell;
	int ai, aj, cnt, cntE; 
	int Wo[NE][NE][SH];
	double deltaE, dE, E1, E2;

        i = (int) (randd1()*alloyState.N_3);
        n = (int) (randd1()*alloyState.Nneighbors);
        j = alloyState.inputPos[i][n];
    	if(alloyState.Atom[i] != alloyState.Atom[j]){
        	E1 = alloyState.currEtot;
        	ai = alloyState.Atom[i];
        	aj = alloyState.Atom[j];
        	memcpy(Wo, alloyState.W, sizeof(int)*NE*NE*SH);
        	updateW(ai, aj, i, j);
        	alloyState.Atom[i] = aj;
        	alloyState.Atom[j] = ai;
        	E2 = Etot();

        	deltaE = E2 - E1;
        	alloyState.attd++;
		if(mode == metropolis){//metropolis
        		if(Metropolis(alloyState.currEtot, alloyState.currEtot+deltaE) == 1){// accept
        			alloyState.currEtot += deltaE; alloyState.accd++;
        		}else{//reject
                		alloyState.Atom[i] = ai;
                		alloyState.Atom[j] = aj;
                		memcpy(alloyState.W, Wo, sizeof(int)*NE*NE*SH);
        		}
		}
		else{//wanglandau
        		if(WangLandau(alloyState.currEtot, alloyState.currEtot+deltaE, mode) == 1){// accept
                		alloyState.currEtot += deltaE; alloyState.accd++;
        		}else{//reject
                		alloyState.Atom[i] = ai;
                		alloyState.Atom[j] = aj;
                		memcpy(alloyState.W, Wo, sizeof(int)*NE*NE*SH);
        		}

		}
    	}

}

void encode(float* z){
	int i,j,k,t; 
#ifdef TF_BACKEND
	tensorflow::TensorShape data_shape({1, alloyState.VAE_D, alloyState.VAE_D, alloyState.VAE_D, NE});
        tensorflow::Tensor data(tensorflow::DT_FLOAT, data_shape);
        auto data_ = data.flat<float>().data();
  	for(i=0;i<alloyState.N;i++)for(j=0; j<alloyState.N;j++)for(k=0;k<alloyState.N;k++){
		int idx = (j-i+k+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*alloyState.VAE_D*NE + (k-j+i+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*NE + (j+i-k+alloyState.SHIFT+alloyState.PAD)*NE;
		for(t =0; t <NE; t++)
			if(alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k]==t){
				data_[idx+t] = 1;
				break;
			}else{
				data_[idx+t] = 0;
			}
	}
	std::vector<tensorflow::Tensor> outputs;
        alloyState.model[0].Predict(data, outputs, PredictMode::encoder);
	for(i=0;i<3;i++)
		z[i] = outputs[0].flat<float>().data()[i];
#else
	std::string z_key, config_key;
	std::string rankID = std::to_string(mpiState.myrank);
	memset(alloyState.inputConfig,0,sizeof(float)*alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE);
  	for(i=0;i<alloyState.N;i++)for(j=0; j<alloyState.N;j++)for(k=0;k<alloyState.N;k++){
		int idx = (j-i+k+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*alloyState.VAE_D*NE + (k-j+i+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*NE + (j+i-k+alloyState.SHIFT+alloyState.PAD)*NE;
		for(t =0; t <NE; t++)
			if(alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k]==t){
				alloyState.inputConfig[idx+t] = 1;
				break;
			}
	}
	z_key = "output_z_"+rankID;
	config_key = "input_config_"+rankID;
        (*SRclient).put_tensor(config_key, alloyState.inputConfig, {1, alloyState.VAE_D, alloyState.VAE_D, alloyState.VAE_D, NE}, SmartRedis::TensorType::flt,
                                SmartRedis::MemoryLayout::contiguous);
	(*SRclient).run_model(encoder_name+rankID, {config_key}, {z_key});
	(*SRclient).unpack_tensor(z_key, z, {3}, SmartRedis::TensorType::flt,SmartRedis::MemoryLayout::contiguous);
#endif
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

	npos[0] += alloyState.Z_R*v1;
	npos[1] += alloyState.Z_R*v2;
	npos[2] += alloyState.Z_R*v3; 

}

int arg_max(std::vector<float> const& vec){
	return static_cast<int>(std::distance(vec.begin(), std::max_element(vec.begin(), vec.end())));
}

void decode(float* z){
  	int i, j, k, t, tt, id, tid;	
	int sum[NE] = {0};
	int sumo[NE] = {0};
	float output[alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE] = {0};
	std::vector<float> vec(NE); 	
#ifdef TF_BACKEND

	tensorflow::TensorShape data_shape({1, 3});
        tensorflow::Tensor data(tensorflow::DT_FLOAT, data_shape);
        auto data_ = data.flat<float>().data();
	for(i=0;i<3;i++)
		data_[i] = z[i];
	std::vector<tensorflow::Tensor> outputs;
        alloyState.model[1].Predict(data, outputs, PredictMode::decoder);
	memcpy(output, outputs[0].flat<float>().data(), sizeof(float)*alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE);
#else
	std::string z_key, config_key;
	std::string rankID = std::to_string(mpiState.myrank);
	
	z_key = "input_z_"+rankID;
	config_key = "output_config_"+rankID;
        (*SRclient).put_tensor(z_key, z, {1,3}, SmartRedis::TensorType::flt,
                                SmartRedis::MemoryLayout::contiguous);
	(*SRclient).run_model(decoder_name+rankID, {z_key}, {config_key});
	(*SRclient).unpack_tensor(config_key, &output, {alloyState.VAE_D*alloyState.VAE_D*alloyState.VAE_D*NE}, SmartRedis::TensorType::flt,SmartRedis::MemoryLayout::contiguous);
#endif
	for(t=0;t<NE;t++)
		memset(alloyState.elist[t],-1,sizeof(int)*alloyState.N_3);
  	for(i=0;i<alloyState.N;i++)for(j=0; j<alloyState.N;j++)for(k=0;k<alloyState.N;k++){
		int idx = (j-i+k+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*alloyState.VAE_D*NE + (k-j+i+alloyState.SHIFT+alloyState.PAD)*alloyState.VAE_D*NE + (j+i-k+alloyState.SHIFT+alloyState.PAD)*NE;
		for(t=0; t<NE; t++)
			vec[t] = output[idx+t];
		t = arg_max(vec);
        	alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k] = t;
		alloyState.elist[t][sum[t]++] = i*alloyState.N_2+j*alloyState.N+k;
	}
	for(t=0;t<NE;t++)
		sumo[t] = sum[t];
	// fix concentration, better to take prob into account 
	for(t=0;t<NE;t++){
		while(sum[t] > alloyState.NT[t]){
			id = int(randd1()*sumo[t]);
			if(alloyState.elist[t][id] != -1){
				for(tt=0;tt<NE;tt++)
					if(tt!=t && sum[tt] < alloyState.NT[tt]){
						tid = tt;
						break;
					}
				alloyState.Atom[alloyState.elist[t][id]] = tid;
				alloyState.elist[t][id] = -1;
				sum[t]--;
				sum[tid]++;
			} 	
		}
	}
	for(t=0;t<NE;t++)
		assert(sum[t] == alloyState.NT[t]);
	ini_W();

}

void vae_update(SamplingMode mode){
	float z[3] = {0};
	double E1, E2, deltaE;
	int Wo[NE][NE][SH];
        memcpy(alloyState.Atomo, alloyState.Atom, sizeof(short)*alloyState.N_3);
        memcpy(Wo, alloyState.W, sizeof(int)*NE*NE*SH);
	// encode the current config to latent space;
	encode(z);
        /*if(mpiState.myrank == 0 && wlState.TotalSweeps%100 == 0){
		printf("step %d before: z=(%f,%f,%f), E=%f\n", wlState.TotalSweeps, z[0], z[1], z[2], alloyState.currEtot);
		//write_xyz(TotalSweeps);
	}*/
	// random walk in latent space; 
	walk(z);
	// decode the data point to real space; 
        decode(z);
	// measure the energy of new configuration;
	E1 = alloyState.currEtot;
	E2 = Etot(); 
        /*if(mpiState.myrank == 0 && wlState.TotalSweeps%100==0){
		printf("step %d after: z=(%f,%f,%f), E=%f\n", wlState.TotalSweeps, z[0], z[1], z[2], E2);
		//write_xyz(TotalSweeps+1);
	}*/
	// update according to WL 	
	deltaE = E2 - E1;
        alloyState.attd++;
	if(mode == 0){//metropolis
        	if(Metropolis(alloyState.currEtot, alloyState.currEtot+deltaE) == 1){// accept
        		alloyState.currEtot += deltaE; alloyState.accd++;
        	}else{//reject
        		memcpy(alloyState.Atom, alloyState.Atomo, sizeof(short)*alloyState.N_3);
        		memcpy(alloyState.W, Wo, sizeof(int)*NE*NE*SH);
		}
	}else{//wanglandau
        	if(WangLandau(alloyState.currEtot, alloyState.currEtot+deltaE, mode) == 1){// accept
               		alloyState.currEtot += deltaE; alloyState.accd++;
        	}else{//reject
        		memcpy(alloyState.Atom, alloyState.Atomo, sizeof(short)*alloyState.N_3);
        		memcpy(alloyState.W, Wo, sizeof(int)*NE*NE*SH);
        	}
	}
}

void write_pos(){
  char s[512];
  FILE *ofp;
  ofp=fopen("compos.dat","w");
	int i,j,k;
	for(i=0;i<alloyState.N;i++)
		for(j=0;j<alloyState.N;j++)
			for(k=0;k<alloyState.N;k++){
				fprintf(ofp, "%d\t%d\t%d\t%d\t \n",i+k, i+j, j+k,alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k]);
			}
  fclose(ofp);
}

void write_xyz(int frame)
{
  int i,j,k,cnt, shell;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"snap_%d_%d.xyz",frame, mpiState.myrank);
  ofp=fopen(s,"ab+");
  fprintf(ofp, "%d\n", alloyState.N_3); 
  fprintf(ofp,"Eng = %.8lg  MC_step = %d\n", alloyState.currEtot, frame);

  for(shell = 0; shell < SH; shell++)
        for(i =0; i < NE-1; i++)
                for(j =i+1; j < NE; j++){
                        fprintf(ofp, "%.8lg\t",1.0*alloyState.W[i][j][shell]/alloyState.NS[shell]/alloyState.N_3);
                }
  fprintf(ofp, "\n");


  for(i=0;i<alloyState.N;i++)
          for(j=0; j<alloyState.N;j++)
                  for(k=0;k<alloyState.N;k++)
                        fprintf(ofp,"%s %d %d %d \n",alloyState.element[alloyState.Atom[i*alloyState.N_2+j*alloyState.N+k]+1], j-i+k, k-j+i, j+i-k);
  fclose(ofp);

}

void thermoqs()
{
        int i;
        double T,U,Z,C,F,S,M,M2,X;
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
        for(T=alloyState.TTi;T<=alloyState.TTf+alloyState.dTT;T=T+alloyState.dTT)
        {
                U = 0.0;        // Initializing the average energy <E>
                Z = 0.0;        // Initializing the Partition Function
                C = 0.0;        // Initialize the Specific Heat
                F = 0.0;        // Initialize the free energy
                S = 0.0;        // Initialize the entropy
                Bw = 0.0;       //Boltzmann weight
		M = M2 = X = 0.0; // order parameter
                lambda = -1.0e300;      //Normalization shift (max value of DOS considering T)
                centerE = 0.0;  // Taking the center of the energy bin
               //Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
                for(i=0;i<wlState.D1BINS;i++)
                {
                        centerE = ( (i/wlState.invdWLD1+wlState.WLD1min) + 0.5*(wlState.WLD1max - wlState.WLD1min)/(1.0*wlState.D1BINS) )/alloyState.invN;

                        if( (lambda < ((wlState.wllng[i]) - 1.0*(centerE*E_scale)/(T*T_scale))) )// && (wllng[i] > 0.0) )
                                lambda = ((wlState.wllng[i]) - 1.0*(centerE*E_scale)/(T*T_scale));

                };

                //Central Loop for calculating thermodynamic properties from the DOS
                for(i=0;i<wlState.D1BINS;i++)
                {
                        //Taking the center of the bin
                        centerE = ( (i/wlState.invdWLD1+wlState.WLD1min) + 0.5*(wlState.WLD1max - wlState.WLD1min)/(1.0*wlState.D1BINS) )/alloyState.invN;

                        //Boltzmann Factor
                        Bw =  exp( wlState.wllng[i]  - (centerE*E_scale)/(T*T_scale) - lambda );

                        //Partition Function
                        Z = Z + Bw;

                        //Average Energy <E>
                        U = U + (centerE)*Bw;

                        //Average Energy Squared <E^2>
                        C = C + (centerE)*(centerE)*Bw;  

			if(wlState.wlH[i] > 0){
				M += wlState.wllngi[i]/wlState.wlH[i]*Bw;
				M2 += wlState.wllngd[i]/wlState.wlH[i]*Bw;
			}                                                      
                };

                //Internal Energy
                U = U/Z;
                //Normalization of free energy
                if(T == alloyState.TTi)
                {
                        Nfree = -(T*T_scale)/E_scale*( lambda + log(Z) ) - U;
                };
                //Specific Heat
                C = ( (C/Z) - (U*U) )*E_scale*E_scale / (T*T*T_scale*T_scale);
                //Free Energy
                F = -(T*T_scale)/E_scale*( lambda + log(Z) ) - Nfree;
                //Entropy
                S = (U - F)*E_scale/(T*T_scale);

		//VAE order parameter
		M /= Z; M2 /= Z;
		X = (M2 - M*M)*alloyState.N_3/T/T_scale; 

                fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",T,U*alloyState.invN,C*alloyState.invN,F*alloyState.invN,S*alloyState.invN,M,X);
        };

        fclose(therm_op);
}

double L1(){
	float z[3];
	encode(z);
	double M = 0.0;
	for(int i=0;i<3;i++)
		M += fabs(z[i]);

	return M;
}

// VAE order parameter 
void OrderParameter(int idx){
	double M;  
	M = L1();
	assert(idx >=0 && idx < wlState.D1BINS);
        alloyState.op[idx] += M;
        alloyState.op2[idx] += M*M;
}



