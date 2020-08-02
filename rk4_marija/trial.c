#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double t, xgrid, ygrid,gam,Fx,Fy,Fz,pE,k,dx,dy,w,Eo,Bo,delta,kdamp, tfwhm, stable;  
int pri, Ni, Nj;


double Poly(double x)
{double res;
 res=10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
 return res;
}

double Envelope(double x,double t){
	double res,x1,x2; 
	x1=(x-t+2*tfwhm+stable)/tfwhm; 
	x2=(t-x)/tfwhm;
	if (x<t-2*tfwhm-stable) res=0;
	else if (x<t-tfwhm-stable) res=Poly(x1);
	else if (x<t-tfwhm) res = 1;
	else if (x<t) res=Poly(x2);
	else res=0; 
	return res;
}


/*
double Ax(double phi){return 0;}
double Ay(double phi){return delta*Eo*sin(phi);}
double Az(double phi){return sqrt(1-delta*delta)*Eo*cos(phi);}

double Phi(double x, double y, double z, double t){
		return w*t - k*x;
	};

double Der(double(*f)(double), double x){
	double h = 1e-5;
	return (f(x+h/2.)-f(x-h/2.))/h; 
};

double Der2(double(*f)(double), double x){
	double h = 1e-5;
	return (f(x+h)-2*f(x)+f(x-h))/h/h;
}

double Efx(double x, double y, double z){return 0;}
double Efy(double x, double y, double z){return -w*Der(Ay,Phi(x,y,z,t)) *Envelope(x,t);}
double Efz(double x, double y, double z){return -w*Der(Az,Phi(x,y,z,t)) *Envelope(x,t);}

double DerEfx(double x, double y, double z){return 0;}
double DerEfy(double x, double y, double z){return -w*w*Der2(Ay,Phi(x,y,z,t)) *Envelope(x,t);}
double DerEfz(double x, double y, double z){return -w*w*Der2(Az,Phi(x,y,z,t)) *Envelope(x,t);}

double Bfx(double x, double y, double z){return 0;}
double Bfy(double x, double y, double z){return k*Der(Az,Phi(x,y,z,t)) *Envelope(x,t);}
double Bfz(double x, double y, double z){return -k*Der(Ay,Phi(x,y,z,t)) *Envelope(x,t);}

double DerBfx(double x, double y, double z){return 0;}
double DerBfy(double x, double y, double z){return w*k*Der2(Az,Phi(x,y,z,t)) *Envelope(x,t);}
double DerBfz(double x, double y, double z){return -w*k*Der2(Ay,Phi(x,y,z,t)) *Envelope(x,t);}


*/

double Efx(double x, double y, double z)  //Ex interpolation to (x,y,z)
{ double res;
 res=0;
 return res; 
}

double DerEfx(double x, double y, double z) //Ex time derivative at (x,y,z)
{double res;
 res=0;
 return res; 
}

double Efy(double x, double y, double z)  //Ey interpolation to (x,y,z)
{double res;
 res=delta*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerEfy(double x, double y, double z) //Ey time derivative at (x,y,z)
{double res;
 res=delta*w*w*Eo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}

double Efz(double x, double y, double z)  //Ez interpolation to (x,y,z)
{ double res;
 res=-sqrt(1-delta*delta)*w*Eo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerEfz(double x, double y, double z) //Ez time derivative at (x,y,z)
{double res;
 res=sqrt(1-delta*delta)*w*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double Bfx( double x, double y, double z)  //Bx interpolation to (x,y,z)
{ double res;
 res=0;
 return res; 
}

double DerBfx(double x, double y, double z) //Bx time derivative at (x,y,z)
{double res;
 res=0;
 return res; 
}

double Bfy(double x, double y, double z)  //By interpolation to (x,y,z)
{ double res;
 res=sqrt(1-delta*delta)*k*Bo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerBfy(double x, double y, double z) //By time derivative at (x,y,z)
{double res;
 res=-sqrt(1-delta*delta)*w*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double Bfz(double x, double y, double z)  //Bz interpolation to (x,y,z)
{double res;
 res=delta*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerBfz(double x, double y, double z) //Bz time derivative at (x,y,z)
{double res;
 res=delta*w*k*Bo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}


double fun1( double px,double py,double pz,double x,double y, double z)   //  dpx/dt
{double res,A,B,C,D;
A=0;// A=gam*DerEfx(x,y,z)+px*(Efx(x+dx,y)-Efx(x-dx,y))/(2*dx)+py*(Efx(x,y+dy)-Efx(x,y-dy))/(2*dy)+py*DerBfz(x,y,z)-pz*DerBfy(x,y,z);
// printf("A=%f\n",A);
B=0;// B=py*px*(Bfz(x+dx,y)-Bfz(x-dx,y))/(2*dx*gam)+py*py*(Bfz(x,y+dy)-Bfz(x,y-dy))/(2*dy*gam)-pz*px*(Bfy(x+dx,y)-Bfy(x-dx,y))/(2*dx*gam)-pz*py*(Bfy(x,y+dy)-Bfy(x,y-dy))/(2*dy*gam);
// printf("=%f\n",gam);
C=0;// C=Efy(x,y,z)*Bfz(x,y,z)-Efz(x,y,z)*Bfy(x,y,z)+(Bfy(x,y,z)*Bfx(x,y,z)*py-Bfy(x,y,z)*Bfy(x,y,z)*px-Bfz(x,y,z)*Bfz(x,y,z)*px+Bfz(x,y,z)*Bfx(x,y,z)*pz)/gam;
D=0;// D=Efx(x,y,z)*pE/gam+px*pE*pE/gam-gam*px*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fx+kdamp*(A+B+C+D);
 return res;

}

double fun2( double px,double py,double pz,double x,double y, double z)   //  dpy/dt
{double res,A,B,C,D;
 A=0;//A=gam*DerEfy(x,y,z)+px*(Efy(x+dx,y)-Efy(x-dx,y))/(2*dx)+py*(Efy(x,y+dy)-Efy(x,y-dy))/(2*dy)+pz*DerBfx(x,y,z)-px*DerBfz(x,y,z);
B=0;// B=pz*px*(Bfx(x+dx,y)-Bfx(x-dx,y))/(2*dx*gam)+pz*py*(Bfx(x,y+dy)-Bfx(x,y-dy))/(2*dy*gam)-px*px*(Bfz(x+dx,y)-Bfz(x-dx,y))/(2*dx*gam)-px*py*(Bfz(x,y+dy)-Bfz(x,y-dy))/(2*dy*gam);
C=0;// C=Efz(x,y,z)*Bfx(x,y,z)-Efx(x,y,z)*Bfz(x,y,z)+(Bfz(x,y,z)*Bfy(x,y,z)*pz-Bfz(x,y,z)*Bfz(x,y,z)*py-Bfx(x,y,z)*Bfx(x,y,z)*py+Bfx(x,y,z)*Bfy(x,y,z)*px)/gam;
D=0;// D=Efy(x,y,z)*pE/gam+py*pE*pE/gam-gam*py*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fy+kdamp*(A+B+C+D);
 return res;

}

double fun3( double px,double py,double pz,double x,double y, double z)   //  dpz/dt
{double res,A,B,C,D;
A=0;// A=gam*DerEfz(x,y,z)+px*(Efz(x+dx,y)-Efz(x-dx,y))/(2*dx)+px*DerBfy(x,y,z)-py*DerBfx(x,y,z);
B=0;// B=px*px*(Bfy(x+dx,y)-Bfy(x-dx,y))/(2*dx*gam)+px*py*(Bfy(x,y+dy)-Bfy(x,y-dy))/(2*dy*gam)-py*px*(Bfx(x+dx,y)-Bfx(x-dx,y))/(2*dx*gam)-py*py*(Bfx(x,y+dy)-Bfx(x,y-dy))/(2*dy*gam);
C=0;// C=Efx(x,y,z)*Bfy(x,y,z)-Efy(x,y,z)*Bfx(x,y,z)+(Bfx(x,y,z)*Bfz(x,y,z)*px-Bfx(x,y,z)*Bfx(x,y,z)*pz-Bfy(x,y,z)*Bfy(x,y,z)*pz+Bfy(x,y,z)*Bfz(x,y,z)*py)/gam;
D=0;// D=Efz(x,y,z)*pE/gam+pz*pE*pE/gam-gam*pz*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fz+kdamp*(A+B+C+D);
 return res;

}


double fun4(double px,double py, double pz, double x, double y, double z)   // dx/dt
{double res;
 res=px/sqrt(1+px*px+py*py+pz*pz);
 return res;
}

double fun5(double px,double py, double pz, double x, double y, double z)    //  dy/dt
{double res;
 res=py/sqrt(1+px*px+py*py+pz*pz);
 return res;
}

double fun6(double px,double py, double pz, double x, double y, double z)    //  dy/dt
{double res;
 res=pz/sqrt(1+px*px+py*py+pz*pz);
 return res;
}


//procedure RK does runge-kutta integration of the model defined with functions above (fun1,2,3,4,5)

double RK(double T,long long int N, double p01, double p02,double p03, double x01, double x02, double x03) 
{double k11, k12, k13, k14,k15, k16, k21, k22, k23, k24, k25, k26, k31, k32, k33, k34, k35, k36, k41, k42, k43, k44, k45, k46;
 double w1, w2, w3, w4, w5, w6, h;
 double px, py, pz, x, y, z;
 long long int i; 
 FILE *fo;
 
 
 fo=fopen("Out0.txt","w");
 //fprintf(fo,"t		 x		 y	   z  	 px		 py		  pz	 gamma\n");
 h=T/(double)N; 
 
 w1=p01; w2=p02; w3=p03; w4=x01; w5=x02; w6=x03; gam=sqrt(1+p01*p01+p02*p02+p03*p03); t=0;
 

 fprintf(fo, "%.14e    	%.14e 	 %.14e   %.14e	%.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w6, w1, w2, w3, gam); 
 for(i=1; i<N; i++)
 	{
 	 {
 	 t=(double)(i-1) * h;
 	
 	 px=w1; py=w2; pz=w3; x=w4; y=w5; z=w6;
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k11=h*fun1(px,py,pz,x,y,z); k12=h*fun2(px,py,pz,x,y,z); k13=h*fun3(px,py,pz,x,y,z); 
 	 k14=h*fun4(px,py,pz,x,y,z); k15=h*fun5(px,py,pz,x,y,z); k16=h*fun6(px,py,pz,x,y,z);
 	 
 	//------------------------------------------------------------------------
 	 
 	 t=(double)(i-1)*h+h/2.;
 
 	 
 	 px=w1+k11*1/2; py=w2+k12*1/2; pz=w3+k13*1/2; x=w4+k14*1/2; y=w5+k15*1/2; z=w6+k16*1/2;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	  
	 k21=h*fun1(px,py,pz,x,y,z); 
 	 k22=h*fun2(px,py,pz,x,y,z);
 	 k23=h*fun3(px,py,pz,x,y,z);
 	 k24=h*fun4(px,py,pz,x,y,z);
 	 k25=h*fun5(px,py,pz,x,y,z);
 	 k26=h*fun6(px,py,pz,x,y,z);

 	 //-----------------------------------------------------------------------
 	 
 	 px=w1+k21*1/2; py=w2+k22*1/2; pz=w3+k23*1/2; x=w4+k24*1/2; y=w5+k25*1/2; z=w6+k26*1/2;
 	
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k31=h*fun1(px,py,pz,x,y,z); 
 	 k32=h*fun2(px,py,pz,x,y,z);
 	 k33=h*fun3(px,py,pz,x,y,z);
 	 k34=h*fun4(px,py,pz,x,y,z);
 	 k35=h*fun5(px,py,pz,x,y,z);
 	 k36=h*fun6(px,py,pz,x,y,z);

 	 //---------------------------------------------------------------------------
 	  
 	 t=(double)i*h;
 	 
 	 px=w1+k31; py=w2+k32; pz=w3+k33; x=w4+k34; y=w5+k35; z=w6+k36;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k41=h*fun1(px,py,pz,x,y,z);  
 	 k42=h*fun2(px,py,pz,x,y,z);
 	 k43=h*fun3(px,py,pz,x,y,z);
 	 k44=h*fun4(px,py,pz,x,y,z);
 	 k45=h*fun5(px,py,pz,x,y,z);
 	 k46=h*fun6(px,py,pz,x,y,z);

 	 
 	 w1=w1+(k11+2*k21+2*k31+k41)/6; 
 	 w2=w2+(k12+2*k22+2*k32+k42)/6; 
 	 w3=w3+(k13+2*k23+2*k33+k43)/6; 
 	 w4=w4+(k14+2*k24+2*k34+k44)/6; 
 	 w5=w5+(k15+2*k25+2*k35+k45)/6;
 	 w6=w6+(k16+2*k26+2*k36+k46)/6;
 	 
 
 	 //---------------------------------------------------------------------------
 	 
 	 if ((i % pri) == 0)
 	 {gam=sqrt(1+w1*w1+w2*w2+w3*w3); 
 	  fprintf(fo, "%.14e    	%.14e 	 %.14e   %.14e	%.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w6, w1, w2, w3, gam); 
 	  // prints t, x, y, z, px, py, pz, gamma
	  }
	 }
 	}
 	fclose(fo);
}

int main()
{ double p01, p02, p03, x01, x02, x03;
  double  T; 
 long long int N; 
 FILE *foo;
 
 foo=fopen("InputTotBatch.txt","r");
 fscanf(foo,"%lf %lf %lf %lf %lf %lf %lf %lf %lld %i %lf %lf ", &p01, &p02, &p03, &Eo, &delta, &kdamp,&k, &T, &N, &pri, &xgrid, &ygrid);
 tfwhm = 100.;
 stable = 100.;
 
 w=k; Bo=Eo; 
 //p01 = -Eo*Eo/4.;
 fclose(foo);
 x01=0.; x02=0.; x03=100.;
 dx=xgrid/2.;dy=ygrid/2.;
 

 RK(T, N, p01, p02, p03, x01, x02, x03);

 printf("gotovo\n");


 }
 
