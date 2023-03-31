#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <omp.h>
int main(){

	int Nx=40,Ny=40,Nt=640; // for a stable system dt=h^2/2*ó, ó=1cm^2/sec, Nt=T/dt
	
	double **U=(double**)malloc(Nx*sizeof(double*));
	double **UN=(double**)malloc(Nx*sizeof(double*));
	double **UNM1=(double**) malloc(Nx*sizeof(double*));
	
	double *X=malloc(Nx*sizeof(double));
	double *Y=malloc(Ny*sizeof(double));
	double *t=malloc(Nt*sizeof(double));
	int i,j,n;
	double c,s=1; //s=ó=1cm^2/sec
	double T,h,t1=151,t2=304,t3=608; //t1=5/dt t=5, t2=10/dt t=10,t3=20/dt t=20 ,for t=0 u(x,y,t)=I(x,y)
	double Lx, Ly;
	double xo,yo,to;
	double dt,dx,dy;
	double Cx2,Cy2,tend,tstart;
    double uth;  // analitic solution
   tstart = omp_get_wtime();
	
	Lx=10.0;
	Ly=10.0;
	c=1.0;
    T=20;	
	h=(Lx-0.0)/(Nx-1);  // h=dx=dy
	dt=h*h/(2*s);
	yo=0.0;
	xo=0.0;
	to=0.0;
	
	
 	 for (i=0; i<Nx; i++){
      	 UNM1[i]=(double*)malloc(Ny*sizeof(double));
      }

  	for (i=0; i<Nx; i++){
      	UN[i]=(double*)malloc(Ny*sizeof(double));
      }

  	
  		for (i=0; i<Nx; i++){
      	U[i]=(double*)malloc(Ny*sizeof(double));
      }
      
      
      
      
   
	  
      
omp_set_dynamic(0);
omp_set_num_threads(1);

#pragma omp parallel private(i,j,n) shared(X,h,Y,t,dt,UN,UNM1,U,Cx2,Cy2,xo,yo,to,Nx,Ny,Nt,Lx,Ly) default(none)  

	
	#pragma omp for 
	
	for(i=0; i<Nx; i++){
		X[i]=xo+i*h;		
	}
	
	#pragma omp for 
	
		for(j=0; j<Ny; j++){
		Y[j]=yo+j*h;
	}
	
	#pragma omp for 
	
		for(n=0; n<=Nt; n++){
		t[n]=n*dt;
	}
	
	
	Cx2=(c*(dt/h))*(c*(dt/h));
	Cy2=(c*(dt/h))*(c*(dt/h));
	
	#pragma omp for 
	
	for (i=0; i<Nx; i++){
		for (j=0;j<Ny; j++){
			UNM1[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j]);      // the time t=0.0           	         
	}
}

	#pragma omp single
      {
			uth=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+0.5*to);		  	
			printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
				n*dt,UNM1[Nx/2][Ny/2],uth);	
			
		}
	
	
	
	#pragma omp for  
	
	for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
			UN[i][j]=0.5* (UNM1[i+1][j]-2*UNM1[i][j]+UNM1[i-1][j])*Cx2  
                      +0.5* (UNM1[i][j+1]-2*UNM1[i][j]+UNM1[i][j-1])*Cy2
                      + UNM1[i][j]+dt*0.5*X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j]);
                      + dt*dt*c*c*(1+1.0/2.0)*(Y[j]*(Ly - Y[j])+ X[i]*(Lx- X[i])); // the time  t=1.0
         }
	    
	}
	

	for (n=2;n<=Nt;n++){
	
#pragma omp for 	 

for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      U[i][j]=(UN[i+1][j]-2*UN[i][j]+UN[i-1][j])*Cx2  
                      +(UN[i][j+1]-2*UN[i][j]+UN[i][j-1])*Cy2
                      + 2*UN[i][j]-UNM1[i][j]
                      + dt*dt*2*c*c*(1+t[n]/2.0)*(Y[j]*(Ly - Y[j])+ X[i]*(Lx- X[i]));
                      
                     
		}
           			
	}
	
	
	
	
	
	//swap 
#pragma omp for 

for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      UNM1[i][j]=UN[i][j];	
  }
}

#pragma omp for 

for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      UN[i][j]=U[i][j];	
  }
}



#pragma omp single
      {
		if (n==t1) {
			uth=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+0.5*t[n]);		  	
			printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
				n*dt,UN[Nx/2][Ny/2],uth);	
			}
		}


#pragma omp single
      {
		if (n==t2) {
			uth=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+0.5*t[n]);		  	
			printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
				n*dt,UN[Nx/2][Ny/2],uth);	
			}
		}


#pragma omp single
      {
		if (n==t3) {
			uth=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+0.5*t[n]);		  	
			printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
				n*dt,UN[Nx/2][Ny/2],uth);	
			}
		}






} // end of time



// free memory

    free(U);
  	free(UN);
  	free(UNM1);
  
  	free(X);
  	free(Y);
  	free(t);

tend = omp_get_wtime();

printf("\n  elapsed time = %.20f, \n", tend - tstart);

return(0);
}

	
	
		


	
