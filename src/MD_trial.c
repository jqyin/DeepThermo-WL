#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include"MD_trial.h"
#include"water.h"
#include "rand.h"
#include "wanglandau.h"

extern int myrank;

void ini_MD(){
	int i;
	srand(time(NULL)*myrank);
	quat = malloc(sizeof(struct Quat)*N);

	F = malloc(sizeof(struct Vec)*N);
	Tr = malloc(sizeof(struct Vec)*N);
	J = malloc(sizeof(struct Vec)*N);
	V = malloc(sizeof(struct Vec)*N);
	WmolOld = malloc(sizeof(struct Water)*N);
	attemptmd = malloc(sizeof(int)*D1BINS);
	acceptmd = malloc(sizeof(int)*D1BINS);
	actmd = malloc(sizeof(int)*D1BINS);
	attmd = malloc(sizeof(int)*D1BINS);
	for(i=0;i<D1BINS;i++)
		attemptmd[i]= acceptmd[i] = 0;

}

void ini_quat(){
	int i;
	double A11, A12, A13, A21, A22, A23, A31,A32,A33;
	double t1, t2, t3,Nr, Xc, Yc, Zc;
	struct Water test[2];



	for(i=0;i<N;i++){
		Xc = Wmol[i].cx;
		Yc = Wmol[i].cy;
		Zc = Wmol[i].cz;

		t1 = (Xc - Wmol[i].Ox);
		t2 = (Yc - Wmol[i].Oy);
		t3 = (Zc - Wmol[i].Oz);
		Nr = sqrt(t1*t1+t2*t2+t3*t3);
		
		A13 = t1/Nr;
		A23 = t2/Nr;
		A33 = t3/Nr;	

		t1 = (Wmol[i].H1x - Wmol[i].H2x);
		t2 =  (Wmol[i].H1y - Wmol[i].H2y);
		t3 =  (Wmol[i].H1z - Wmol[i].H2z);
		Nr = sqrt(t1*t1+t2*t2+t3*t3);

		A12 = t1/Nr;
		A22 = t2/Nr;
		A32 = t3/Nr;	

		A11 = A22*A33 - A32*A23;
		A21 = A32*A13 - A12*A33;
		A31 = A12*A23 - A22*A13;

		if( A11+A22+A33 > 0){
			Nr = sqrt(1.0+A11+A22+A33)*2;
			quat[i].q0 = 0.25*Nr;
			quat[i].q1 = (A32 - A23)/Nr;
			quat[i].q2 = (A13 - A31)/Nr;
			quat[i].q3 = (A21 - A12)/Nr;
		}else if( (A11>A22) & (A11 > A33)){
			Nr = sqrt(1.0+A11-A22-A33)*2;
			quat[i].q0 = (A32 - A23)/Nr;
			quat[i].q1 = 0.25*Nr;
			quat[i].q2 = (A12 + A21)/Nr;
			quat[i].q3 = (A13 + A31)/Nr;		
		}else if( A22 > A33){
			Nr = sqrt(1.0+A22-A11-A33)*2;
			quat[i].q0 = (A13 - A31)/Nr;
			quat[i].q1 = (A12 + A21)/Nr;
			quat[i].q2 = 0.25*Nr;
			quat[i].q3 = (A23 + A32)/Nr;				
		}else{
			Nr = sqrt(1.0+A33-A11-A22)*2;
			quat[i].q0 = (A21 - A12)/Nr;
			quat[i].q1 = (A13 + A31)/Nr;
			quat[i].q2 = (A23 + A32)/Nr;
			quat[i].q3 = 0.25*Nr;	
		}

	}
//test
/*
	for(i=0;i<N;i++){
		Xc = Rc[i].x;
		Yc = Rc[i].y;
		Zc = Rc[i].z;

	    A11 = quat[i].q0 *quat[i].q0 +quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
        A12 = 2.0 * (quat[i].q1*quat[i].q2+ quat[i].q0 *quat[i].q3);
        A13 = 2.0 * (quat[i].q1*quat[i].q3- quat[i].q0 *quat[i].q2);
        A21 = 2.0 * (quat[i].q1*quat[i].q2- quat[i].q0 *quat[i].q3);
        A22 = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 +quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
        A23 = 2.0 * (quat[i].q2*quat[i].q3+ quat[i].q0 *quat[i].q1);
        A31 = 2.0 * (quat[i].q1*quat[i].q3+ quat[i].q0 *quat[i].q2);
        A32 = 2.0 * (quat[i].q2*quat[i].q3- quat[i].q0 *quat[i].q1);
        A33 = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 +quat[i].q3*quat[i].q3;
			
	   test[i].Ox = (Xc -0.0207858 * A31) ;
	   test[i].Oy = (Yc-0.0207858 * A32) ;
	   test[i].Oz = (Zc-0.0207858 * A33) ;

	   test[i].H1x = (Xc + 0.239997 * A21 +0.1649726  * A31) ;
	   test[i].H1y = (Yc +0.239997 * A22 +0.1649726  * A32) ;
	   test[i].H1z = (Zc +0.239997 * A23 +0.1649726  * A33) ;

	   test[i].H2x = (Xc-0.239997 * A21 +0.1649726  * A31) ;
	   test[i].H2y = (Yc-0.239997 * A22 +0.1649726  * A32) ;
	   test[i].H2z = (Zc-0.239997 * A23 +0.1649726  * A33) ;	
	} */
}

void FT(){
	int i,j,iti, nr, nl;
	double roo6, r12, beta;
	double temp,t, t1, t2, t3;
	double xm[2],ym[2],zm[2],xc[2],yc[2],zc[2];
	double temp1, temp2;
	double** FOx,  **FOy, **FOz;
	double** FMx,  **FMy, **FMz;
	double** FH1x,  **FH1y, **FH1z;
	double** FH2x,  **FH2y, **FH2z;
	double* FaOx, *FaOy,*FaOz,* FaMx, *FaMy,*FaMz,* FaH1x, *FaH1y,*FaH1z,* FaH2x, *FaH2y,*FaH2z;

	FaOx = malloc(sizeof(double)*N);
	FaOy = malloc(sizeof(double)*N);
	FaOz = malloc(sizeof(double)*N);

	FaMx = malloc(sizeof(double)*N);
	FaMy = malloc(sizeof(double)*N);
	FaMz = malloc(sizeof(double)*N);

	FaH1x = malloc(sizeof(double)*N);
	FaH1y = malloc(sizeof(double)*N);
	FaH1z = malloc(sizeof(double)*N);

	FaH2x = malloc(sizeof(double)*N);
	FaH2y = malloc(sizeof(double)*N);
	FaH2z = malloc(sizeof(double)*N);

	FOx = (double**) malloc(sizeof(double*)*N);
	FOy = (double**) malloc(sizeof(double*)*N);
	FOz = (double**) malloc(sizeof(double*)*N);

	FMx = (double**) malloc(sizeof(double*)*N);
	FMy = (double**) malloc(sizeof(double*)*N);
	FMz = (double**) malloc(sizeof(double*)*N);

	FH1x = (double**) malloc(sizeof(double*)*N);
	FH1y = (double**) malloc(sizeof(double*)*N);
	FH1z = (double**) malloc(sizeof(double*)*N);

	FH2x = (double**) malloc(sizeof(double*)*N);
	FH2y = (double**) malloc(sizeof(double*)*N);
	FH2z = (double**) malloc(sizeof(double*)*N);

	for(i=0;i<N;i++){
		FOx[i] = malloc(sizeof(double)*N);
		FOy[i] = malloc(sizeof(double)*N);
		FOz[i] = malloc(sizeof(double)*N);

		FMx[i] = malloc(sizeof(double)*N);
		FMy[i] = malloc(sizeof(double)*N);
		FMz[i] = malloc(sizeof(double)*N);

		FH1x[i] = malloc(sizeof(double)*N);
		FH1y[i] = malloc(sizeof(double)*N);
		FH1z[i] = malloc(sizeof(double)*N);

		FH2x[i] = malloc(sizeof(double)*N);
		FH2y[i] = malloc(sizeof(double)*N);
		FH2z[i] = malloc(sizeof(double)*N);
	}

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			FOx[i][j] =FOy[i][j]= FOz[i][j]=0.0;
			FMx[i][j] =FMy[i][j]= FMz[i][j]=0.0;
			FH1x[i][j] =FH1y[i][j]= FH1z[i][j]=0.0;
			FH2x[i][j] =FH2y[i][j]= FH2z[i][j]=0.0;
		}
		F[i].x = F[i].y = F[i].z = 0.0;
		Tr[i].x = Tr[i].y = Tr[i].z = 0.0;
		FaOx[i] = FaOy[i] = FaOz[i] = 0.0;
		FaMx[i] = FaMy[i] = FaMz[i] = 0.0;
		FaH1x[i] = FaH1y[i] = FaH1z[i] = 0.0;
		FaH2x[i] = FaH2y[i] = FaH2z[i] = 0.0;
	}

	for(i=0 ; i< N-1; i++){
		for(j=i+1; j< N;j++){
			t1 =  (Wmol[i].Ox- Wmol[j].Ox);
			t2 =  (Wmol[i].Oy- Wmol[j].Oy);
			t3 =  (Wmol[i].Oz- Wmol[j].Oz);
			temp = t1*t1+t2*t2+t3*t3;
			roo6 = temp*temp*temp;
			temp1 = 48*(1/roo6 -0.5)/roo6/temp;
			FOx[i][j] += t1*temp1;
			FOy[i][j] += t2*temp1;
			FOz[i][j] += t3*temp1;


#ifdef TIP4P
		   xc[0] = Wmol[i].cx;
		   yc[0] = Wmol[i].cy;
		   zc[0] = Wmol[i].cz;

		   xc[1] = Wmol[j].cx;
		   yc[1] = Wmol[j].cy;
		   zc[1] = Wmol[j].cz;

	//       t = sqrt(  (Wmol[i].Ox - xc[0])*(Wmol[i].Ox - xc[0]) + (Wmol[i].Oy - yc[0])*(Wmol[i].Oy - yc[0]) + (Wmol[i].Oz - zc[0])*(Wmol[i].Oz - zc[0]) );
		   t = Rom/Roc;
		   xm[0] = Wmol[i].Ox - (Wmol[i].Ox - xc[0])*t;  // line equation in 3D;
		   ym[0] = Wmol[i].Oy - (Wmol[i].Oy - yc[0])*t;
		   zm[0] = Wmol[i].Oz - (Wmol[i].Oz - zc[0])*t;

	 //      t = sqrt(  (Wmol[j].Ox - xc[1])*(Wmol[j].Ox - xc[1]) + (Wmol[j].Oy - yc[1])*(Wmol[j].Oy - yc[1]) + (Wmol[j].Oz - zc[1])*(Wmol[j].Oz - zc[1]) );
		   t = Rom/Roc;
		   xm[1] = Wmol[j].Ox - (Wmol[j].Ox - xc[1])*t;  // line equation in 3D;
		   ym[1] = Wmol[j].Oy - (Wmol[j].Oy - yc[1])*t;
		   zm[1] = Wmol[j].Oz - (Wmol[j].Oz - zc[1])*t;

			t1 = (xm[0]- xm[1]);
			t2 = (ym[0]- ym[1]);
			t3 = (zm[0]- zm[1]);
			r12 = t1*t1 + t2*t2 + t3*t3;
			r12 = sqrt(r12)*r12;
			temp2 = 4.0/r12;
			FMx[i][j] += t1*temp2*b;
			FMy[i][j] += t2*temp2*b;
			FMz[i][j] += t3*temp2*b;

			t1 =  (Wmol[i].H1x- xm[1]);
			t2 =  (Wmol[i].H1y- ym[1]);
			t3 =  (Wmol[i].H1z- zm[1]);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= -2.0/r12;
			FH1x[i][j] += t1*temp2*b;
			FH1y[i][j] += t2*temp2*b;
			FH1z[i][j] += t3*temp2*b;

		    t1 =  (Wmol[i].H2x- xm[1]);
			t2 = (Wmol[i].H2y- ym[1]);
			t3 =  (Wmol[i].H2z- zm[1]);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= -2.0/r12;
			FH2x[i][j] += t1*temp2*b;
			FH2y[i][j] += t2*temp2*b;
			FH2z[i][j] += t3*temp2*b;

			t1 = (xm[0]- Wmol[j].H1x);
			t2 = (ym[0]- Wmol[j].H1y);
			t3 = (zm[0]- Wmol[j].H1z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= -2.0/r12;
			FMx[i][j] += t1*temp2*b;
			FMy[i][j] += t2*temp2*b;
			FMz[i][j] += t3*temp2*b;

			 t1 = (Wmol[i].H1x- Wmol[j].H1x);
			 t2 =  (Wmol[i].H1y- Wmol[j].H1y);
			 t3 = (Wmol[i].H1z- Wmol[j].H1z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= 1.0/r12;
			FH1x[i][j] += t1*temp2*b;
			FH1y[i][j] += t2*temp2*b;
			FH1z[i][j] += t3*temp2*b;

			t1 = (Wmol[i].H2x- Wmol[j].H1x);
			t2 =  (Wmol[i].H2y- Wmol[j].H1y);
			t3 =  (Wmol[i].H2z- Wmol[j].H1z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= 1.0/r12;
			FH2x[i][j] += t1*temp2*b;
			FH2y[i][j] += t2*temp2*b;
			FH2z[i][j] += t3*temp2*b;


			t1 =  (xm[0]- Wmol[j].H2x);
			t2 =  (ym[0]- Wmol[j].H2y);
			t3 =  (zm[0]- Wmol[j].H2z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= -2.0/r12;
			FMx[i][j] += t1*temp2*b;
			FMy[i][j] += t2*temp2*b;
			FMz[i][j] += t3*temp2*b;

			t1 =  (Wmol[i].H1x- Wmol[j].H2x);
			t2 =  (Wmol[i].H1y- Wmol[j].H2y);
			t3 =  (Wmol[i].H1z- Wmol[j].H2z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= 1.0/r12;
			FH1x[i][j] += t1*temp2*b;
			FH1y[i][j] += t2*temp2*b;
			FH1z[i][j] += t3*temp2*b;

			t1 =  (Wmol[i].H2x- Wmol[j].H2x);
			t2 =  (Wmol[i].H2y- Wmol[j].H2y);
			t3 =  (Wmol[i].H2z- Wmol[j].H2z);
		    r12 =  t1*t1 + t2*t2 + t3*t3;
			r12 =sqrt(r12)*r12;
			temp2= 1.0/r12;
			FH2x[i][j] += t1*temp2*b;
			FH2y[i][j] += t2*temp2*b;
			FH2z[i][j] += t3*temp2*b;

			FOx[j][i] = -FOx[i][j];
			FOy[j][i] = -FOy[i][j];
			FOz[j][i] = -FOz[i][j];

			FMx[j][i] = -FMx[i][j];
			FMy[j][i] = -FMy[i][j];
			FMz[j][i] = -FMz[i][j];

			FH1x[j][i] = -FH1x[i][j];
			FH1y[j][i] = -FH1y[i][j];
			FH1z[j][i] = -FH1z[i][j];

			FH2x[j][i] = -FH2x[i][j];
			FH2y[j][i] = -FH2y[i][j];
			FH2z[j][i] = -FH2z[i][j];
#endif

		}
	}

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(j!=i){
				FaOx[i] +=  FOx[i][j]  ;
				FaOy[i] +=  FOy[i][j] ;
				FaOz[i] +=  FOz[i][j]  ;

				FaMx[i] +=  FMx[i][j]  ;
				FaMy[i] +=  FMy[i][j] ;
				FaMz[i] +=  FMz[i][j]  ;

				FaH1x[i] +=  FH1x[i][j]  ;
				FaH1y[i] +=  FH1y[i][j] ;
				FaH1z[i] +=  FH1z[i][j]  ;

				FaH2x[i] +=  FH2x[i][j]  ;
				FaH2y[i] +=  FH2y[i][j] ;
				FaH2z[i] +=  FH2z[i][j]  ;
			}
		}
		F[i].x += FaOx[i] + FaMx[i]  + FaH1x[i] +FaH2x[i] ;
		F[i].y += FaOy[i] + FaMy[i]  + FaH1y[i] +FaH2y[i] ;
		F[i].z += FaOz[i] + FaMz[i]  + FaH1z[i] +FaH2z[i] ;
	}

	for(i=0;i<N;i++){
		t1 = Wmol[i].Ox - Wmol[i].cx;
		t2 = Wmol[i].Oy - Wmol[i].cy;
		t3 = Wmol[i].Oz - Wmol[i].cz;
		Tr[i].x += t2*FaOz[i] - t3*FaOy[i];
		Tr[i].y += t3*FaOx[i] - t1*FaOz[i];
		Tr[i].z += t1*FaOy[i] - t2*FaOx[i];


	//       t = sqrt(  (Wmol[i].Ox - Rc[i].x)*(Wmol[i].Ox - Rc[i].x) + (Wmol[i].Oy -  Rc[i].y)*(Wmol[i].Oy - Rc[i].y) + (Wmol[i].Oz - Rc[i].z)*(Wmol[i].Oz -  Rc[i].z) );
	   t = Rom/Roc;
	   xm[0] = Wmol[i].Ox - (Wmol[i].Ox - Wmol[i].cx)*t;  // line equation in 3D;
	   ym[0] = Wmol[i].Oy - (Wmol[i].Oy -	Wmol[i].cy)*t;
	   zm[0] = Wmol[i].Oz - (Wmol[i].Oz - Wmol[i].cz)*t;

		t1 = xm[0] - Wmol[i].cx;
		t2 = ym[0] - Wmol[i].cy;
		t3 = zm[0] - Wmol[i].cz;
		Tr[i].x += t2*FaMz[i] - t3*FaMy[i];
		Tr[i].y += t3*FaMx[i] - t1*FaMz[i];
		Tr[i].z += t1*FaMy[i] - t2*FaMx[i];

		t1 = Wmol[i].H1x - Wmol[i].cx;
		t2 = Wmol[i].H1y - Wmol[i].cy;
		t3 = Wmol[i].H1z - Wmol[i].cz;
		Tr[i].x += t2*FaH1z[i] - t3*FaH1y[i];
		Tr[i].y += t3*FaH1x[i] - t1*FaH1z[i];
		Tr[i].z += t1*FaH1y[i] - t2*FaH1x[i];

		t1 = Wmol[i].H2x - Wmol[i].cx;
		t2 = Wmol[i].H2y - Wmol[i].cy;
		t3 = Wmol[i].H2z - Wmol[i].cz;
		Tr[i].x += t2*FaH2z[i] - t3*FaH2y[i];
		Tr[i].y += t3*FaH2x[i] - t1*FaH2z[i];
		Tr[i].z += t1*FaH2y[i] - t2*FaH2x[i];
	}

	iti=(int) ((currEtot*invN-WLD1min)*invdWLD1);

// using coeff generated from savgol procedure;

/*	case(0): // nl=0; nr=4;   ||  -0.77143  0.18571  0.57143  0.38571  -0.37143;
			sum = -0.77143*wllng[iti] + 0.18571*wllng[iti+1] + 0.57143*wllng[iti+2] + 0.38571*wllng[iti+3] - 0.37143*wllng[iti+4];
			break;
		case(1):// nl=1; nr=3; -0.48571  ||  0.04286  0.28571  0.24286  -0.08571; 
			sum = -0.48571*wllng[iti-1] +  0.04286*wllng[iti] + 0.28571*wllng[iti+1] + 0.24286*wllng[iti+2] - 0.08571*wllng[iti+3];
			break;
		case(2:(D1BINS-3) ):// nl=2; nr=2; -0.2  -0.1  || 0.0  || 0.1  0.2;
			sum = -0.2*wllng[iti-2] - 0.1*wllng[iti-1] + 0.1*wllng[iti+1] + 0.2*wllng[iti+2] ;
			break;
		case(D1BINS-2):// nl=3; nr=1; 0.08571  -0.24286  -0.28571  -0.04286 ||  0.48571;
			sum = 0.08571*wllng[iti-3] - 0.24286*wllng[iti-2] - 0.28571*wllng[iti-1] - 0.04286 *wllng[iti] + 0.48571*wllng[iti+1];
			break;
		case(D1BINS-1)://nl=4; nr=0; 0.37143  -0.38571  -0.57143  -0.18571  0.77143 || ;
			sum = 0.37143*wllng[iti-4] - 0.38571*wllng[iti-3] - 0.57143*wllng[iti-2] - 0.18571 *wllng[iti-1] + 0.77143*wllng[iti];
			break;
		default:
			exit(77);
	} */


// 5 point;
	if(iti >= 2 && iti <= (D1BINS-3) ){
			beta = -0.2*wllng[iti-2] - 0.1*wllng[iti-1] + 0.1*wllng[iti+1] + 0.2*wllng[iti+2] ;
	}else if( iti == 0){
			beta = -0.77143*wllng[iti] + 0.18571*wllng[iti+1] + 0.57143*wllng[iti+2] + 0.38571*wllng[iti+3] - 0.37143*wllng[iti+4];	
	}else if(iti == 1){
			beta = -0.48571*wllng[iti-1] +  0.04286*wllng[iti] + 0.28571*wllng[iti+1] + 0.24286*wllng[iti+2] - 0.08571*wllng[iti+3];	
	}else if(iti == (D1BINS-2) ){
			beta = 0.08571*wllng[iti-3] - 0.24286*wllng[iti-2] - 0.28571*wllng[iti-1] - 0.04286 *wllng[iti] + 0.48571*wllng[iti+1];	
	}else if(iti == (D1BINS-1) ){
			beta = 0.37143*wllng[iti-4] - 0.38571*wllng[iti-3] - 0.57143*wllng[iti-2] - 0.18571 *wllng[iti-1] + 0.77143*wllng[iti];	
	}else{
		exit(77);
	}	

/*
if(lnwlf > 1e-4){
	if(iti >= 2 && iti <= (D1BINS-3) ){
		nl = 2; nr = 2;
		beta = coef[2][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[2][i+2]*wllng[iti - i - 1];
			beta += coef[2][Order-i]*wllng[iti + i +1];
		}
	}else if( iti == 0){
		nl = 0; nr = 4;
		beta = coef[0][1]*wllng[iti]; 
		for(i = 0; i< nr; i++){
			beta += coef[0][Order-i]*wllng[iti + i +1];
		}
	}else if(iti == 1){
		nl = 1; nr = 3;
		beta = coef[1][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[1][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[1][Order-i]*wllng[iti + i +1];
		}
	}else if(iti == (D1BINS-2) ){
		nl = 3; nr = 1;
		beta = coef[3][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[3][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[3][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-1) ){
		nl = 4; nr = 0;
		beta = coef[4][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[4][i+2]*wllng[iti - i - 1];
		}
	}else{
		exit(77);
	}
}else{  // stage II 

	if(iti > 0 && iti < D1BINS-1){
		beta = (wllng[iti+1] - wllng[iti-1])*0.5;
	}else if(iti == 0){
		beta = -1.5*wllng[iti] +2.0*wllng[iti+1]-0.5*wllng[iti+2];
	}else{
		beta = 0.5*wllng[iti-2] -2.0* wllng[iti-1] + 1.5*wllng[iti];
	}
	
}
*/


//11 point;
/*
if(lnwlf > 1e-4){
	if(iti >= 5 && iti <= (D1BINS-6) ){
		nl = 5; nr = 5;
		beta = coef[5][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[5][i+2]*wllng[iti - i - 1];
			beta += coef[5][Order-i]*wllng[iti + i +1];
		}

	}else if( iti == 0){

		nl = 0; nr = 10;
		beta = coef[0][1]*wllng[iti];  
		for(i = 0; i< nr; i++){
			beta += coef[0][Order-i]*wllng[iti + i +1];
		}

	}else if(iti == 1){

		nl = 1; nr = 9;
		beta = coef[1][1]*wllng[iti];  
		for(i = 0; i< nl; i++){
			beta +=  coef[1][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[1][Order-i]*wllng[iti + i +1];
		}

	}else if(iti == 2){

		nl = 2; nr = 8;
		beta = coef[2][1]*wllng[iti];  
		for(i = 0; i< nl; i++){
			beta +=  coef[2][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[2][Order-i]*wllng[iti + i +1];
		}

	}else if(iti == 3){

		nl = 3; nr = 7;
		beta = coef[3][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[3][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[3][Order-i]*wllng[iti + i +1];
		}

	}else if(iti==4){
		nl = 4; nr = 6;
		beta = coef[4][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[4][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[4][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-5) ){
		nl = 6; nr = 4;
		beta = coef[6][1]*wllng[iti];  
		for(i = 0; i< nl; i++){
			beta +=  coef[6][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[6][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-4) ){
		nl = 7; nr = 3;
		beta = coef[7][1]*wllng[iti];  
		for(i = 0; i< nl; i++){
			beta +=  coef[7][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[7][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-3) ){
		nl = 8; nr = 2;
		beta = coef[8][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[8][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[8][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-2) ){
		nl = 9; nr = 1;
		beta = coef[9][1]*wllng[iti]; 
		for(i = 0; i< nl; i++){
			beta +=  coef[9][i+2]*wllng[iti - i - 1];
		}
		for(i = 0; i< nr; i++){
			beta += coef[9][Order-i]*wllng[iti + i +1];
		}	
	}else if(iti == (D1BINS-1) ){
		nl = 10; nr = 0;
		beta = coef[10][1]*wllng[iti];  
		for(i = 0; i< nl; i++){
			beta +=  coef[10][i+2]*wllng[iti - i - 1];
		}
	}else{
		exit(77);
	} 	
}else{  // stage II 

	if(iti > 0 && iti < D1BINS-1){
		beta = (wllng[iti+1] - wllng[iti-1])*0.5;
	}else if(iti == 0){
		beta = wllng[iti+1] - wllng[iti];
	}else{
		beta = wllng[iti] - wllng[iti-1];
	}
	
}*/
    beta /= (dWLD1/4.186/Ue)*Ut; 



/*
	if(iti <  3){
	//	t1 = (wllng[iti] +wllng[i+1] +wllng[i+2]+wllng[i+3]+wllng[i+4] )/5; moving average
	//	t2 = (wllng[i+1] +wllng[i+2]+wllng[i+3]+wllng[i+4]+wllng[i+5] )/5;
		beta =  (wllng[iti+5] - wllng[iti])/(dWLD1/4.184/Ue)/5.0;      //dWLD1 is in unit of kJ/mol ;
	}
	else if(iti > D1BINS-4 )
		beta = (wllng[iti] - wllng[iti-5])/(dWLD1/4.184/Ue)/5.0;
	else{
//		t1 = (wllng[i-3] + wllng[i-2] +wllng[i-1] +wllng[i]+wllng[i+1] )/5;
//		t2 = (wllng[i-1] + wllng[i] +wllng[i+1] +wllng[i+2]+wllng[i+3] )/5;
		beta =  (wllng[iti+2]+wllng[iti+3] - wllng[iti-3] - wllng[iti-2])/5.0/(dWLD1/4.184/Ue)/2.0;
	}
	if(beta < 0.0)
		beta = 0.0;
	else
		beta /= Ut;

//	beta = 1.0;
*/
	for(i=0;i<N;i++){
		F[i].x *= beta;
		F[i].y *= beta;
		F[i].z *= beta;

		Tr[i].x *= beta;
		Tr[i].y *= beta;
		Tr[i].z *= beta;
	}

	free(FaOx);
	free(FaOy);
	free(FaOz);

	free(FaMx);
	free(FaMy);
	free(FaMz);

	free(FaH1x);
	free(FaH1y);
	free(FaH1z);

	free(FaH2x);
	free(FaH2y);
	free(FaH2z);
	for(i=0;i<N;i++){
		free(FOx[i]);
		free(FOy[i]);
		free(FOz[i]);

		free(FMx[i]);
		free(FMy[i]);
		free(FMz[i]);

		free(FH1x[i]);
		free(FH1y[i]);
		free(FH1z[i]);

		free(FH2x[i]);
		free(FH2y[i]);
		free(FH2z[i]);
	}
	free(FOx);
	free(FOy);
	free(FOz);

	free(FMx);
	free(FMy);
	free(FMz);

	free(FH1x);
	free(FH1y);
	free(FH1z);

	free(FH2x);
	free(FH2y);
	free(FH2z);
}


double ini_p(){
	int i;
	double K=0.0;
	double AXX,AXY,AXZ,AYX,AYY,AYZ,AZX,AZY,AZZ;
	double OxI, OyI, OzI;
	double sumx, sumy,sumz;

	sumx=sumy=sumz=0.0;
	for(i=0;i<N;i++){
//		J[i].x = J[i].y = J[i].z = 0.0;
//		V[i].x = V[i].y = V[i].z = 0.0;
		V[i].x = gaussrand();
		V[i].y = gaussrand();
		V[i].z = gaussrand();  

		sumx += V[i].x;
		sumy += V[i].y;
		sumz += V[i].z;

/*		   do{
				eta1 = 1.0 - 2.0*randd1();
				eta2 = 1.0 - 2.0*randd1();
				etasq = eta1*eta1 + eta2*eta2;
		   }while(etasq >1.0);
			//these are the new unit vectors
		J[i].x = 2.0*eta1*sqrt(1.0-etasq);
		J[i].y = 2.0*eta2*sqrt(1.0-etasq);
		J[i].z = 1.0 - 2.0*etasq;
			*/
		
		J[i].x = gaussrand()*sqrt(Ixx);
		J[i].y = gaussrand()*sqrt(Iyy);
		J[i].z = gaussrand()*sqrt(Izz);

	}
	for(i=0;i<N;i++){

		   V[i].x -= sumx/N;
		   V[i].y -= sumy/N;
		   V[i].z -= sumz/N;

			K += 0.5*(V[i].x*V[i].x + V[i].y*V[i].y + V[i].z*V[i].z);
			K += 0.5*(J[i].x*J[i].x/Ixx +J[i].y*J[i].y/Iyy +J[i].z*J[i].z/Izz);

			AXX = quat[i].q0 *quat[i].q0 +quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			AXY = 2.0 * (quat[i].q1*quat[i].q2+ quat[i].q0 *quat[i].q3);
			AXZ = 2.0 * (quat[i].q1*quat[i].q3- quat[i].q0 *quat[i].q2);
			AYX = 2.0 * (quat[i].q1*quat[i].q2- quat[i].q0 *quat[i].q3);
			AYY = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 +quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			AYZ = 2.0 * (quat[i].q2*quat[i].q3+ quat[i].q0 *quat[i].q1);
			AZX = 2.0 * (quat[i].q1*quat[i].q3+ quat[i].q0 *quat[i].q2);
			AZY = 2.0 * (quat[i].q2*quat[i].q3- quat[i].q0 *quat[i].q1);
			AZZ = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 +quat[i].q3*quat[i].q3;

	//       ** CONVERT ANGULAR MOMENTUM TO frame-FIXED **
			   OxI = ( AXX * J[i].x + AYX * J[i].y + AZX * J[i].z );
			   OyI = ( AXY * J[i].x + AYY * J[i].y + AZY * J[i].z );
			   OzI = ( AXZ * J[i].x + AYZ * J[i].y + AZZ * J[i].z );

			   J[i].x = OxI;
			   J[i].y = OyI;
			   J[i].z = OzI;

	}

	return K;
}

void MD(){
	double DT2, SDT, Ei, Ef;
	double AXX,AXY,AXZ,AYX,AYY,AYZ,AZX,AZY,AZZ;
	double JxI, JyI, JzI, OxI, OyI, OzI, VxI, VyI, VzI;
	double QWI, QXI,QYI,QZI,QW1I,QX1I,QY1I,QZ1I;
	double tx,ty,tz, rHH, rhh;
	double Ki, Kf, temp;
	int i, j, t, iti, sign;
	
//	assert(fabs(Etot(Wmol)-currEtot) < 1e-7);
	if(randd1() < 0.5)
		sign = 1;
	else
		sign =-1;

	SDT = sign*DT;
	DT2 = SDT/2.0;
	
	for(i=0;i<N;i++){
	/*	WmolOld[i].Ox=Wmol[i].Ox;
		WmolOld[i].Oy=Wmol[i].Oy;
		WmolOld[i].Oz=Wmol[i].Oz;
		WmolOld[i].H1x=Wmol[i].H1x;
		WmolOld[i].H1y=Wmol[i].H1y;
		WmolOld[i].H1z=Wmol[i].H1z;
		WmolOld[i].H2x=Wmol[i].H2x;
		WmolOld[i].H2y=Wmol[i].H2y;
		WmolOld[i].H2z=Wmol[i].H2z;
		WmolOld[i].cx=Wmol[i].cx;
		WmolOld[i].cy=Wmol[i].cy;
		WmolOld[i].cz=Wmol[i].cz;*/
		WmolOld[i] = Wmol[i];

		Wmol[i].Ox/=Ul;
		Wmol[i].Oy/=Ul;
		Wmol[i].Oz/=Ul;
	    Wmol[i].H1x/=Ul;
		Wmol[i].H1y/=Ul;
		Wmol[i].H1z/=Ul;
		Wmol[i].H2x/=Ul;
		Wmol[i].H2y/=Ul;
		Wmol[i].H2z/=Ul;
		Wmol[i].cx/=Ul;
		Wmol[i].cy/=Ul;
		Wmol[i].cz/=Ul;
	}

	for(t=0;t<Nt; t++){

		ini_quat();
		Ki = ini_p();
		Kf = 0.0;
		FT();

		for(i=0;i<N;i++){
	//      ** AUXILIARY EQUATION MOVES   **
	//       ** ANGULAR MOMENTUM TO TIME T **
			JxI =J[i].x + DT2*Tr[i].x;
			JyI =J[i].y + DT2*Tr[i].y;
			JzI =J[i].z + DT2*Tr[i].z;
	//       ** OBTAIN ROTATION MATRIX AT TIME T **
			AXX = quat[i].q0 *quat[i].q0 +quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			AXY = 2.0 * (quat[i].q1*quat[i].q2+ quat[i].q0 *quat[i].q3);
			AXZ = 2.0 * (quat[i].q1*quat[i].q3- quat[i].q0 *quat[i].q2);
			AYX = 2.0 * (quat[i].q1*quat[i].q2- quat[i].q0 *quat[i].q3);
			AYY = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 +quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			AYZ = 2.0 * (quat[i].q2*quat[i].q3+ quat[i].q0 *quat[i].q1);
			AZX = 2.0 * (quat[i].q1*quat[i].q3+ quat[i].q0 *quat[i].q2);
			AZY = 2.0 * (quat[i].q2*quat[i].q3- quat[i].q0 *quat[i].q1);
			AZZ = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 +quat[i].q3*quat[i].q3;

	//       ** CONVERT ANGULAR MOMENTUM TO BODY-FIXED **
	//       ** FORM AND HENCE TO ANGULAR VELOCITIES   **
			   OxI = ( AXX * JxI + AXY * JyI + AXZ * JzI ) / Ixx;
			   OyI = ( AYX * JxI + AYY * JyI + AYZ * JzI ) / Iyy;
			   OzI = ( AZX * JxI + AZY * JyI + AZZ * JzI ) / Izz;

			   Kf +=  0.5*(Ixx * OxI *OxI + Iyy * OyI *OyI + Izz * OzI *OzI);
	//       ** OBTAIN TIME-DERIVATIVES OF QUATERNIONS **
	//       ** AND ADVANCE TO TIME T+DT/2             **
			   QW1I = ( -quat[i].q1 * OxI -quat[i].q2 * OyI -quat[i].q3 * OzI ) * 0.5;
			   QX1I = (   quat[i].q0 * OxI -quat[i].q3 * OyI +quat[i].q2 * OzI ) * 0.5;
			   QY1I = (   quat[i].q3 * OxI + quat[i].q0 * OyI - quat[i].q1 * OzI ) * 0.5;
			   QZ1I = (  -quat[i].q2 * OxI + quat[i].q1 * OyI + quat[i].q0 * OzI ) * 0.5;
			   QWI = quat[i].q0 + DT2 * QW1I;
			   QXI = quat[i].q1 + DT2 * QX1I;
			   QYI =quat[i].q2 + DT2 * QY1I;
			   QZI =quat[i].q3 + DT2 * QZ1I;


	//       ** OBTAIN ROTATION MATRIX AT TIME T+DT/2 **

			   AXX = QWI *QWI + QXI *QXI - QYI *QYI - QZI *QZI;
			   AXY = 2.0 * ( QXI * QYI + QWI * QZI );
			   AXZ = 2.0 * ( QXI * QZI - QWI * QYI );
			   AYX = 2.0 * ( QXI * QYI - QWI * QZI );
			   AYY = QWI *QWI - QXI *QXI + QYI *QYI - QZI *QZI;
			   AYZ = 2.0 * ( QYI * QZI + QWI * QXI );
			   AZX = 2.0 * ( QXI * QZI + QWI * QYI );
			   AZY = 2.0 * ( QYI * QZI - QWI * QXI );
			   AZZ = QWI *QWI - QXI *QXI - QYI *QYI + QZI *QZI;

	//       ** MOVE THE ANGULAR MOMENTA ALL THE WAY     **
	//       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY     **
	//       ** CONVERT TO BODY-FIXED ANGULAR VELOCITIES **
	//       ** AT TIME T+DT/2                           **

			   J[i].x +=  SDT * Tr[i].x;
			   J[i].y +=  SDT * Tr[i].y;
			   J[i].z +=  SDT * Tr[i].z;

			   OxI = ( AXX * J[i].x + AXY * J[i].y + AXZ * J[i].z ) / Ixx;
			   OyI = ( AYX * J[i].x + AYY * J[i].y + AYZ * J[i].z ) / Iyy;
			   OzI = ( AZX * J[i].x + AZY * J[i].y + AZZ * J[i].z ) / Izz;

	//		   ** OBTAIN TIME-DERIVATIVES OF QUATERNIONS **
	//       ** AND ADVANCE TO T+DT                    **

			   QW1I = ( - QXI * OxI - QYI * OyI - QZI * OzI ) * 0.5;
			   QX1I = (   QWI * OxI - QZI * OyI + QYI * OzI ) * 0.5;
			   QY1I = (   QZI * OxI + QWI * OyI - QXI * OzI ) * 0.5;
			   QZ1I = ( - QYI * OxI + QXI * OyI + QWI * OzI ) * 0.5;
			   quat[i].q0 +=  SDT * QW1I;
			   quat[i].q1 +=  SDT * QX1I;
			   quat[i].q2 +=  SDT * QY1I;
			   quat[i].q3 +=  SDT * QZ1I;


		// renormalize 
			   if( (temp=quat[i].q1*quat[i].q1+quat[i].q2*quat[i].q2+quat[i].q3*quat[i].q3) > 1.0){
						quat[i].q0 = 0.0;
					  temp = sqrt(temp);
					 quat[i].q1 /= temp;
					 quat[i].q2 /= temp;
					 quat[i].q3 /= temp;

			   }else
					quat[i].q0 =(quat[i].q0>0?1:-1)* sqrt(1-temp);
	

		//  assert(fabs( quat[i].q0 *quat[i].q0 + quat[i].q1 *quat[i].q1 +quat[i].q2 *quat[i].q2 +quat[i].q3 *quat[i].q3 -1) < 1e-3);
	//  body fixed frame  rO = (0, 0, -0.0207858),  rM=(0,0,0.02677286), rH =(0,+-0.239997,0.1649726 );
			AXX = quat[i].q0 *quat[i].q0 +quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
		   AXY = 2.0 * (quat[i].q1*quat[i].q2+ quat[i].q0 *quat[i].q3);
		   AXZ = 2.0 * (quat[i].q1*quat[i].q3- quat[i].q0 *quat[i].q2);
			AYX = 2.0 * (quat[i].q1*quat[i].q2- quat[i].q0 *quat[i].q3);
			AYY = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 +quat[i].q2*quat[i].q2 -quat[i].q3*quat[i].q3;
			AYZ = 2.0 * (quat[i].q2*quat[i].q3+ quat[i].q0 *quat[i].q1);
			AZX = 2.0 * (quat[i].q1*quat[i].q3+ quat[i].q0 *quat[i].q2);
			AZY = 2.0 * (quat[i].q2*quat[i].q3- quat[i].q0 *quat[i].q1);
			AZZ = quat[i].q0 *quat[i].q0 -quat[i].q1*quat[i].q1 -quat[i].q2*quat[i].q2 +quat[i].q3*quat[i].q3;

		   Wmol[i].Ox = ( Wmol[i].cx -0.0207858 * AZX) ;
		   Wmol[i].Oy = ( Wmol[i].cy-0.0207858 *  AZY) ;
		   Wmol[i].Oz = ( Wmol[i].cz-0.0207858 *  AZZ) ;

		   Wmol[i].H1x = ( Wmol[i].cx + 0.239997 * AYX +0.1649726  * AZX) ;
		   Wmol[i].H1y = ( Wmol[i].cy +0.239997 * AYY +0.1649726  * AZY) ;
		   Wmol[i].H1z = ( Wmol[i].cz +0.239997 * AYZ +0.1649726  * AZZ) ;

		   Wmol[i].H2x = ( Wmol[i].cx-0.239997 * AYX +0.1649726  * AZX) ;
		   Wmol[i].H2y = ( Wmol[i].cy-0.239997 * AYY +0.1649726  * AZY) ;
		   Wmol[i].H2z = ( Wmol[i].cz-0.239997 * AYZ +0.1649726  * AZZ) ;	


	//      ** MOVE THE LINEAR VELOCITIES ALL THE WAY   **
	//       ** FROM T-DT/2 TO T+DT/2 AND STORE AWAY     **
			  
			  VxI = V[i].x;
			  VyI = V[i].y;
			  VzI = V[i].z;

			  V[i].x = VxI +  SDT * F[i].x; // /(Mo+2*Mh)=1;
			  V[i].y = VyI +  SDT * F[i].y;
			  V[i].z = VzI +  SDT * F[i].z;

			  VxI = 0.5*(VxI + V[i].x);
			  VyI = 0.5*(VyI + V[i].y);
			  VzI = 0.5*(VzI + V[i].z);

			   Kf +=  0.5* ( VxI*VxI + VyI *VyI + VzI *VzI );
	//       ** ADVANCE POSITIONS TO T+DT **
			   tx = SDT * V[i].x;
			   ty = SDT * V[i].y;
			   tz = SDT * V[i].z;

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


/*	#ifndef NDEBUG
					rHH =Roh*Roh*2 - 2*Roh*Roh*cos(Angle*Pi/180.0);
					rhh = (Wmol[i].H1x-Wmol[i].H2x)*(Wmol[i].H1x-Wmol[i].H2x) + (Wmol[i].H1y-Wmol[i].H2y)*(Wmol[i].H1y-Wmol[i].H2y) +(Wmol[i].H1z-Wmol[i].H2z)*(Wmol[i].H1z-Wmol[i].H2z);
					if(fabs(rHH-rhh) > 1e-5){
						exit(123);
					}


	#endif   */
		}
}


	Ei = currEtot;

	iti=(int) ((currEtot*invN-WLD1min)*invdWLD1);
	attemptmd[iti]++;

	if(constrain_MD()==1){

			for(i=0;i<N;i++){
				Wmol[i].Ox*=Ul;
				Wmol[i].Oy*=Ul;
				Wmol[i].Oz*=Ul;
				Wmol[i].H1x*=Ul;
				Wmol[i].H1y*=Ul;
				Wmol[i].H1z*=Ul;
				Wmol[i].H2x*=Ul;
				Wmol[i].H2y*=Ul;
				Wmol[i].H2z*=Ul;
				Wmol[i].cx*=Ul;
				Wmol[i].cy*=Ul;
				Wmol[i].cz*=Ul;
			}
			 Ef=Etot(Wmol);

			if(WangLandau(Ei,Ef, Ki,Kf,0.0,0.0,-1 )==1)  // K in the unit of Ue? or should it be dimensionless since it's gonna be on the shoulder of a exponential.
			{
					acceptmd[iti]++;
					currEtot = Ef;	
					for(i=0;i<N;i++)
						for(j=0;j<N;j++)
							Epold[i][j] = Epold[j][i] = Ep[i][j];
			}
			else
			{
				//reject - return monomer to old position
				for(i=0;i<N;i++){
/*					Wmol[i].Ox=WmolOld[i].Ox;
					Wmol[i].Oy=WmolOld[i].Oy;
					Wmol[i].Oz=WmolOld[i].Oz;

					Wmol[i].H1x=WmolOld[i].H1x;
					Wmol[i].H1y=WmolOld[i].H1y;
					Wmol[i].H1z=WmolOld[i].H1z;

					Wmol[i].H2x=WmolOld[i].H2x;
					Wmol[i].H2y=WmolOld[i].H2y;
					Wmol[i].H2z=WmolOld[i].H2z;

					Wmol[i].cx=WmolOld[i].cx;
					Wmol[i].cy=WmolOld[i].cy;
					Wmol[i].cz=WmolOld[i].cz; */
					Wmol[i] = WmolOld[i];

				}
				currEtot = Ei;
					for(i=0;i<N;i++)
						for(j=0;j<N;j++)
							Ep[i][j] = Ep[j][i] = Epold[i][j];
			};		


	}else{
					//reject - return monomer to old position
               WangLandau(Ei,Ei, 0.0,0.0,0.0,0.0,-1  ); // update old state;
				for(i=0;i<N;i++){
/*					Wmol[i].Ox=WmolOld[i].Ox;
					Wmol[i].Oy=WmolOld[i].Oy;
					Wmol[i].Oz=WmolOld[i].Oz;

					Wmol[i].H1x=WmolOld[i].H1x;
					Wmol[i].H1y=WmolOld[i].H1y;
					Wmol[i].H1z=WmolOld[i].H1z;

					Wmol[i].H2x=WmolOld[i].H2x;
					Wmol[i].H2y=WmolOld[i].H2y;
					Wmol[i].H2z=WmolOld[i].H2z;

					Wmol[i].cx=WmolOld[i].cx;
					Wmol[i].cy=WmolOld[i].cy;
					Wmol[i].cz=WmolOld[i].cz; */
					Wmol[i] = WmolOld[i];

				}
				currEtot = Ei;	
					for(i=0;i<N;i++)
						for(j=0;j<N;j++)
							Ep[i][j] = Ep[j][i] = Epold[i][j];

	}

}


int constrain_MD(){
	double r;
	int i,j;

	r = Lc0*Lc0;
	for(i = 0; i< N-1; i++)
		for(j = i+1; j<N;j++)
			if((Wmol[j].cx-Wmol[i].cx)*(Wmol[j].cx-Wmol[i].cx) >r ||(Wmol[j].cy-Wmol[i].cy)*(Wmol[j].cy-Wmol[i].cy) > r || (Wmol[j].cz-Wmol[i].cz)*(Wmol[j].cz-Wmol[i].cz) >r)
				return 0;
	return 1;

}


double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}
