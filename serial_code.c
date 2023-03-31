#include<stdio.h>
#include<stdlib.h>
#include <math.h>

int main(){
	FILE *time0;
	FILE *time5;
	FILE *time10;
	FILE *time20;
	int Nx=40,Ny=40,Nt=640; // for a stable system dt=h^2/2*ó, ó=1cm^2/sec, Nt=T/dt
	
	double **U=(double**)malloc(Nx*sizeof(double*));
	double **UN=(double**)malloc(Nx*sizeof(double*));
	double **UNM1=(double**) malloc(Nx*sizeof(double*));
	
	double *X=malloc(Nx*sizeof(double));
	double *Y=malloc(Ny*sizeof(double));
	double *t=malloc(Nt*sizeof(double));
	int i,j,n;
	double c,s=1; //s=ó=1cm^2/sec
	double T,h,t1=151,t2=304,t3=608;
	double Lx, Ly;
	double xo,yo,to;
	double dt,dx,dy;
	double Cx2,Cy2;

	time0=fopen("time0.txt","w");
	time5=fopen("time5.txt","w");
	time10=fopen("time10.txt","w");
	time20=fopen("time20.txt","w");
	
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
      
      
	
	
	for(i=0; i<Nx; i++){
		X[i]=xo+i*h;		
	}
	
		for(j=0; j<Ny; j++){
		Y[j]=yo+j*h;
	}
	
	
		for(n=0; n<=Nt; n++){
		t[n]=n*dt;
	}
	
	
	Cx2=(c*(dt/h))*(c*(dt/h));
	Cy2=(c*(dt/h))*(c*(dt/h));
	
	
	for (i=0; i<Nx; i++){
		for (j=0;j<Ny; j++){
			UNM1[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j]);      // the time t=0.0           	         
		fprintf(time0, "%f\t",UNM1[i][j]);
		}
	fprintf(time0, "\n");
	}
	
	
	
	for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
			UN[i][j]=0.5* (UNM1[i+1][j]-2.*UNM1[i][j]+UNM1[i-1][j])*Cx2  
                      +0.5* (UNM1[i][j+1]-2.*UNM1[i][j]+UNM1[i][j-1])*Cy2
                      + UNM1[i][j]+dt*0.5*X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j])
                      + dt*dt*c*c*(1.+1.0/2.0)*(Y[j]*(Ly - Y[j])+ X[i]*(Lx- X[i])); // the time  t=1.0 
         }
	    
	}
	
	
	for (n=2;n<=Nt;n++){
	
	
for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      U[i][j]=(UN[i+1][j]-2*UN[i][j]+UN[i-1][j])*Cx2  
                      +(UN[i][j+1]-2*UN[i][j]+UN[i][j-1])*Cy2
                      + 2.*UN[i][j]-UNM1[i][j]
                      + dt*dt*2*c*c*(1.+t[n]/2.0)*(Y[j]*(Ly - Y[j])+ X[i]*(Lx- X[i]));
                      
                     
		}
           			
	}
	
	
	
	//swap 

for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      UNM1[i][j]=UN[i][j];	
  }
}


for (i=1; i<Nx-1; i++){
		for (j=1;j<Ny-1; j++){
		      UN[i][j]=U[i][j];	
  }
}





if (n==t1){
	for (i=0; i<Nx; i++){
		for (j=0;j<Ny; j++){
			
				fprintf(time5, "%f\t",U[i][j]); // time 5/dt
		}
	fprintf(time5, "\n");	
    }  

}




if (n==t2){
	for (i=0; i<Nx; i++){
		for (j=0;j<Ny; j++){
			
				fprintf(time10, "%f\t",U[i][j]);  //time 10/dt
		}
	fprintf(time10, "\n");	
    }  

}



if (n==t3){
	for (i=0; i<Nx; i++){
		for (j=0;j<Ny; j++){
			
				fprintf(time20, "%f\t",U[i][j]);   //time 20/dt
		}
	fprintf(time20, "\n");	
    }  

}



} // end of time


fclose(time0);
fclose(time5);
fclose(time10);
fclose(time20);

// free memory

    free(U);
  	free(UN);
  	free(UNM1);
  
  	free(X);
  	free(Y);
  	free(t);



return(0);
}

	
	
		


	
