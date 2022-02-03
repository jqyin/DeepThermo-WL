#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define  DEFINE_GLOBALS
#include "alloy.hpp"
#include "rand.hpp"
#include "wanglandau.hpp"
//#include "cell.cch"
extern int myrank;

void initialize(const char* filename){
   	int i,j,k;
	FILE* fp;
	
	fp = fopen("source_code/data/dos16.dat", "r");
	if(fp == NULL){
		printf("can't open dos64.dat");
		exit(1);
	}
	gexact = (double*)malloc( 256 * sizeof(double ) );
	for(i = 0; i< 256; i++)
		fscanf(fp, "%lg", &gexact[i]);

	S = (int**) malloc(sizeof(int*) *N);
	for(i = 0; i<N; i++)
		S[i] = (int*) malloc(sizeof(int) *N);

	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			S[i][j] = 1;

	invN=1.0/(N*N);
	currEtot=Etot();
	
}

double Etot(){
	int i,j,k;
	double E;

	E=0.0;
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
				E += 0.5*Esite(i,j);

		return E;

}

void count_nn(int t, FILE* fp){
	double err,norm;
	double Hmin, Hmax, Havg;
	int cnt=0, i;
	err = 0;
#ifdef PROF
	norm = wllng[0]-log(prof[0])-gexact[0];
#else
	norm = wllng[0]-gexact[0];
#endif

	Hmin = 1e100;
	Hmax = -1;
	Havg = 0;
	for(i=0;i<D1BINS;i++){
		if(mask[i] == 1){
#ifdef PROF
			err += fabs(1-(wllng[i]-log(prof[i]) - norm)/gexact[i]);
#else
			wllng[i] -= norm; 
			err += fabs(1-wllng[i]/gexact[i]);
#endif
			if(Hmin > wlH[i])
				Hmin = wlH[i];
			if(Hmax < wlH[i])
				Hmax = wlH[i];
			Havg += wlH[i];

			cnt++;
		}
	}
	Havg /= cnt;
	fprintf(fp, "%d\t%g\t%g\t%g\t%g\t%g\n", t, err/cnt, lnwlf, Havg, (Hmax - Hmin), (Hmax - Hmin)/Havg );
	fflush(fp);
	//Gerr[t/1000] = err / cnt; 
	//Herr[t/1000] = (Hmax - Hmin)/Havg;

}
double Esite(int i, int j){
	double E;
	
	E = - S[i][j]*(S[BC(i+1)][j]+S[i][BC(j+1)]+S[BC(i-1)][j]+S[i][BC(j-1)]);

	return E;
}



void Rot(){
		int i, j;
		double E1, E2, deltaE;

		i = (int) (randd1()*N);
		j = (int) (randd1()*N);

		E1 = Esite(i,j);
		S[i][j] = -S[i][j];
		E2 = Esite(i,j);
		deltaE = E2 - E1;
	
		if(WangLandau(currEtot, currEtot+deltaE) == 1){// accept
			currEtot += deltaE;
			if(lnwlf < ModFactorFinal)
				write_conf();
		}else{//reject
			S[i][j] = -S[i][j];
		}


}

void write_conf(){
  char s[512];
  FILE *ofp;
  sprintf(s,"conf%d.dat",myrank);
  ofp=fopen(s,"a");
	int i,j,k;
	fprintf(ofp, "%g\n",currEtot);
	for(i=0;i<N;i++)
		for(j=0;j<N;j++){
				fprintf(ofp, "%d ",S[i][j]);
			}
	fprintf(ofp,"\n");
//	fclose(fop);
}
void write_mol2(int){

}

//Calculates thermodynamic quantities and writes modification-factor labeled files
void Mag(int iti){
	int i, j, k, n;
//	double delta = sizeHistMag / 2.0; //(maxMag - minMag) ;
	double M, M2, M4;
	double Mx, My, Mz;
	Mx = My = Mz =0.0;
	M = M2 = M4 =0.0;
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
				M += S[i][j];

	M *= invN;
	M2 = M*M;
	M4 = M2*M2;
	Mavg[iti] += M;
	Mavg2[iti] += M2;
	Mavg4[iti] += M4;

/*	n = delta*(M +1); //- minMag);
	++ histMag[ (n >=0 && n < sizeHistMag)? n : (n < 0)? 0 : (sizeHistMag-1)][iti];
	n = delta*2*M2;
	++ histMag2[ (n >=0 && n < sizeHistMag)? n : (n < 0)? 0 : (sizeHistMag-1)][iti];
	n = delta*2*M4;
	++ histMag4[ (n >=0 && n < sizeHistMag)? n : (n < 0)? 0 : (sizeHistMag-1)][iti];

	++ cntMag[iti];
	*/
}


void Tavg(){
	int i,j,k,n;	
	double T,Z, delta, temp;
	double centerE,lambda,Bw;
//	double Mag[sizeHistMag];
//	double Mag2[sizeHistMag];
//	double Mag4[sizeHistMag];
	double M, M2, M4;
	double X, BC;
	char s1[512];
	FILE* fp = fopen("m.dat","w");

	
	//Main Temperature Loop
	for(T=TTi;T<=TTf+dTT;T=T+dTT)
	{	

		M=M2=M4=0.0;
		X=BC=0.0;
/*		for(i=0; i<sizeHistMag; i++){
			Mag[i] = Mag2[i] = Mag4[i] = 0.0;
		}*/

		Z = 0.0;  //Initialize the partition function
		Bw = 0.0;	//Boltzmann weight
		lambda = -1.0e300;	//Normalization shift (max value of DOS considering T)
		centerE = 0.0;	// Taking the center of the energy bin

		//Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
		for(i=0;i<D1BINS;i++)
		{
			centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;
				
			if( (lambda < ((wllng[i]) - 1.0*(centerE)/(T*T_scale)  )) ) 
				lambda = ((wllng[i]) - 1.0*(centerE)/(T*T_scale) );  
		};
		
		//Central Loop for calculating thermodynamic properties from the DOS
		for(i=0;i<D1BINS;i++)
		{			
			//Taking the center of the bin
			centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;
			
			//Boltzmann Factor
			Bw =  exp( wllng[i]  - (centerE)/(T*T_scale)  - lambda );
				
			//Partition Function
			Z = Z + Bw;
			
			if(wlH[i]>0.0)		
			{
				//magnetization;	
//				for(k=0; k<sizeHistMag; k++)
					if(wlH[i] > 0){
						M += Mavg[i]/wlH[i]*Bw;
						M2 += Mavg2[i]/wlH[i]*Bw;
						M4 += Mavg4[i]/wlH[i]*Bw;
						//Mag[k] += 1.0*histMag[k][i]/cntMag[i]*Bw;
						//Mag2[k] +=1.0*histMag2[k][i]/cntMag[i]*Bw;
						//Mag4[k] += 1.0*histMag4[k][i]/cntMag[i]*Bw;
					}

			};
		
		};

/*		for(k=0; k<sizeHistMag; k++){
			M += Mag[k];
			M2 += Mag2[k];
			M4 += Mag4[k];
		}
		for(k=0; k<sizeHistMag; k++){
			Mag[k] /= M;
			Mag2[k] /= M2;
			Mag4[k] /= M4;
		}
//		if(T = 0.5) 
//			write_Mdist(Mag, T);

		delta = 2.0/sizeHistMag;
		M = M2 = M4 = 0.0;
		for(k=0; k<sizeHistMag; k++){
			temp = (k+0.5)*delta -1;
			M += fabs(temp)*Mag[k];
			temp = (k+0.5)*delta/2;
			M2 += temp*Mag2[k];
			M4 += temp*Mag4[k];
		}*/

		M /= Z; M2 /= Z; M4 /= Z;

		X = (M2 - M*M)*N*N*N/T/T_scale;
		BC = 1 - M4 / (M2*M2) / 3.0;

		fprintf(fp,"%g\t%g\t%g\t%g\n",T,M,X,BC);		

	};	
		fclose(fp);
	
}



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
			centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;
				
			if( (lambda < ((wllng[i]) - 1.0*(centerE)/(T*T_scale))) )// && (wllng[i] > 0.0) ) 
				lambda = ((wllng[i]) - 1.0*(centerE)/(T*T_scale));  
			
		};
		
		//Central Loop for calculating thermodynamic properties from the DOS
		for(i=0;i<D1BINS;i++)
		{			
			//Taking the center of the bin
			centerE = ( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) )/invN;
			
			//Boltzmann Factor
			Bw =  exp( wllng[i]  - (centerE)/(T*T_scale) - lambda );
		
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
			Nfree = -(T*T_scale)*( lambda + log(Z) ) - U;
		};
		//Specific Heat
		C = ( (C/Z) - (U*U) ) / (T*T*T_scale*T_scale);	
		//Free Energy
		F = -(T*T_scale)*( lambda + log(Z) ) - Nfree;
		//Entropy
		S = (U - F)/T/T_scale;
			
		fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\n",T,U*invN,C*invN,F*invN,S*invN);	
	};	
		
	fclose(therm_op);
}

void freeS(){
	free(S[0]);
	free(S);
	/*free(Sb[0][0][0]);
	free(Sb[0][0]);
	free(Sb[0]);
	free(Sb);

	free(pos[0][0][0]);
	free(pos[0][0]);
	free(pos[0]);
	free(pos);
	free(posb[0][0][0]);
	free(posb[0][0]);
	free(posb[0]);
	free(posb);
	free(Jr);
//	free(Jr[0][0][0]);
//	free(Jr[0][0]);
//	free(Jr[0]);
//	free(Jr);
	free(vol[0][0][0]);
	free(vol[0][0]);
	free(vol[0]);
	free(vol);
	free(volb[0][0][0]);
	free(volb[0][0]);
	free(volb[0]);
	free(volb);*/
}

/* create int [x][y][r][c] */

double ****create_a_foo ( int max_x, int max_y, int max_r, int max_c ) {
    double ****all_x = (double ****)malloc( max_x * sizeof *all_x );
    double  ***all_y =(double ***) malloc( max_x * max_y * sizeof *all_y );
    double   **all_r = (double **)malloc( max_x * max_y * max_r * sizeof *all_r );
    double    *all_c = (double *)malloc( max_x * max_y * max_r * max_c * sizeof *all_c );
    double ****result = all_x;
    int x, y, r;

    for ( x = 0 ; x < max_x ; x++, all_y += max_y ) {
        result[x] = all_y;
        for ( y = 0 ; y < max_y ; y++, all_r += max_r ) {
            result[x][y] = all_r;
            for ( r = 0 ; r < max_r ; r++, all_c += max_c ) {
                result[x][y][r] = all_c;
            }
        }
    }

    return result;
}
