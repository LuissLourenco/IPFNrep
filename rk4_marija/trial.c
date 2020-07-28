#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double t, xgrid, ygrid,gam,Fx,Fy,Fz,pE,k,dx,dy,w,Eo,Bo,kdamp, tfwhm;  
int pri, Ni, Nj;


double Poly(double x)
{double res;
 res=10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
 return res;
}

double Envelope(double x,double t)
{double res,x1,x2; 
 x1=(x-t+2*tfwhm)/tfwhm; 
 x2=(t-x)/tfwhm;
 if (x<t-2*tfwhm) res=0;
    else if (x<t-tfwhm) res=Poly(x1);
          else if (x<t) res=Poly(x2);
               else res=0; 
               
 return res;
 }
 

double Efx(double x, double y)  //Ex interpolation to (x,y)
{ double res;
 res=0;
 return res; 
}

double DerEfx(double x, double y) //Ex time derivative at (x,y)
{double res;
 res=0;
 return res; 
}

double Efy(double x, double y)  //Ey interpolation to (x,y)
{double res;
 res=Eo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerEfy(double x, double y) //Ey time derivative at (x,y)
{double res;
 res=-Eo*w*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double Efz(double x, double y)  //Ez interpolation to (x,y)
{ double res;
 res=0;
 return res; 
}

double DerEfz(double x, double y) //Ez time derivative at (x,y)
{double res;
 res=0;
 return res; 
}

double Bfx( double x, double y)  //Bx interpolation to (x,y)
{ double res;
 res=0;
 return res; 
}

double DerBfx(double x, double y) //Bx time derivative at (x,y)
{double res;
 res=0;
 return res; 
}

double Bfy(double x, double y)  //By interpolation to (x,y)
{ double res;
 res=0;
 return res; 
}

double DerBfy(double x, double y) //By time derivative at (x,y)
{double res;
 res=0;
 return res; 
}

double Bfz(double x, double y)  //Bz interpolation to (x,y)
{double res;
 res=Bo*cos(w*t-k*x)*Envelope(x,t);
 return res; 
}

double DerBfz(double x, double y) //Bz time derivative at (x,y)
{double res;
 res=-Eo*w*sin(w*t-k*x)*Envelope(x,t);
 return res; 
}

double fun1( double px,double py,double pz,double x,double y)   //  dpx/dt
{double res,A,B,C,D;
A=0;// A=gam*DerEfx(x,y)+px*(Efx(x+dx,y)-Efx(x-dx,y))/(2*dx)+py*(Efx(x,y+dy)-Efx(x,y-dy))/(2*dy)+py*DerBfz(x,y)-pz*DerBfy(x,y);
// printf("A=%f\n",A);
B=0;// B=py*px*(Bfz(x+dx,y)-Bfz(x-dx,y))/(2*dx*gam)+py*py*(Bfz(x,y+dy)-Bfz(x,y-dy))/(2*dy*gam)-pz*px*(Bfy(x+dx,y)-Bfy(x-dx,y))/(2*dx*gam)-pz*py*(Bfy(x,y+dy)-Bfy(x,y-dy))/(2*dy*gam);
// printf("=%f\n",gam);
C=0;// C=Efy(x,y)*Bfz(x,y)-Efz(x,y)*Bfy(x,y)+(Bfy(x,y)*Bfx(x,y)*py-Bfy(x,y)*Bfy(x,y)*px-Bfz(x,y)*Bfz(x,y)*px+Bfz(x,y)*Bfx(x,y)*pz)/gam;
D=0;// D=Efx(x,y)*pE/gam+px*pE*pE/gam-gam*px*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fx+kdamp*(A+B+C+D);
 return res;

}

double fun2( double px,double py,double pz,double x,double y)   //  dpy/dt
{double res,A,B,C,D;
 A=0;//A=gam*DerEfy(x,y)+px*(Efy(x+dx,y)-Efy(x-dx,y))/(2*dx)+py*(Efy(x,y+dy)-Efy(x,y-dy))/(2*dy)+pz*DerBfx(x,y)-px*DerBfz(x,y);
B=0;// B=pz*px*(Bfx(x+dx,y)-Bfx(x-dx,y))/(2*dx*gam)+pz*py*(Bfx(x,y+dy)-Bfx(x,y-dy))/(2*dy*gam)-px*px*(Bfz(x+dx,y)-Bfz(x-dx,y))/(2*dx*gam)-px*py*(Bfz(x,y+dy)-Bfz(x,y-dy))/(2*dy*gam);
C=0;// C=Efz(x,y)*Bfx(x,y)-Efx(x,y)*Bfz(x,y)+(Bfz(x,y)*Bfy(x,y)*pz-Bfz(x,y)*Bfz(x,y)*py-Bfx(x,y)*Bfx(x,y)*py+Bfx(x,y)*Bfy(x,y)*px)/gam;
D=0;// D=Efy(x,y)*pE/gam+py*pE*pE/gam-gam*py*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fy+kdamp*(A+B+C+D);
 return res;

}

double fun3( double px,double py,double pz,double x,double y)   //  dpz/dt
{double res,A,B,C,D;
A=0;// A=gam*DerEfz(x,y)+px*(Efz(x+dx,y)-Efz(x-dx,y))/(2*dx)+px*DerBfy(x,y)-py*DerBfx(x,y);
B=0;// B=px*px*(Bfy(x+dx,y)-Bfy(x-dx,y))/(2*dx*gam)+px*py*(Bfy(x,y+dy)-Bfy(x,y-dy))/(2*dy*gam)-py*px*(Bfx(x+dx,y)-Bfx(x-dx,y))/(2*dx*gam)-py*py*(Bfx(x,y+dy)-Bfx(x,y-dy))/(2*dy*gam);
C=0;// C=Efx(x,y)*Bfy(x,y)-Efy(x,y)*Bfx(x,y)+(Bfx(x,y)*Bfz(x,y)*px-Bfx(x,y)*Bfx(x,y)*pz-Bfy(x,y)*Bfy(x,y)*pz+Bfy(x,y)*Bfz(x,y)*py)/gam;
D=0;// D=Efz(x,y)*pE/gam+pz*pE*pE/gam-gam*pz*(Fx*Fx+Fy*Fy+Fz*Fz);
 res=Fz+kdamp*(A+B+C+D);
 return res;

}


double fun4(double px,double py, double pz, double x, double y)   // dx/dt
{double res;
 res=px/sqrt(1+px*px+py*py+pz*pz);
 return res;
}

double fun5(double px,double py, double pz, double x, double y)    //  dy/dt
{double res;
 res=py/sqrt(1+px*px+py*py+pz*pz);
 return res;
}


//procedure RK does runge-kutta integration of the model defined with functions above (fun1,2,3,4,5)

double RK(double T,long long int N, double p01, double p02,double p03, double x01, double x02) 
{double k11, k12, k13, k14,k15, k21, k22, k23, k24,k25, k31, k32, k33, k34, k35, k41, k42, k43, k44, k45;
 double w1, w2, w3, w4, w5, h;
 double px, py, pz, x, y;
 long long int i; 
 FILE *fo;
 
 
 fo=fopen("Out0.txt","w");
 //fprintf(fo,"t		 x		 y		 px		 py		  pz	 gamma\n");
 h=T/N; 
 
 w1=p01; w2=p02; w3=p03; w4=x01; w5=x02;  gam=sqrt(1+p01*p01+p02*p02+p03*p03); t=0;
 

 fprintf(fo, "%.14e    	%.14e 	 %.14e  	%.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w1, w2, w3, gam); 
 for(i=1; i<N; i++)
 	{
 	 {
 	 t=(i-1)*h;
 	
 	 px=w1; py=w2; pz=w3; x=w4; y=w5;
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y)+py*Efy(x,y)+pz*Efz(x,y); 
 	 Fx=Efx(x,y)+(py*Bfz(x,y)-pz*Bfy(x,y))/gam; 	 
 	 Fy=Efy(x,y)+(pz*Bfx(x,y)-px*Bfz(x,y))/gam; 
 	 Fz=Efz(x,y)+(px*Bfy(x,y)-py*Bfx(x,y))/gam; 
 	 
 	 k11=h*fun1(px,py,pz,x,y); k12=h*fun2(px,py,pz,x,y); k13=h*fun3(px,py,pz,x,y); 
 	 k14=h*fun4(px,py,pz,x,y); k15=h*fun5(px,py,pz,x,y);
 	 
 	//------------------------------------------------------------------------
 	 
 	 t=(i-1)*h+h/2;
 
 	 
 	 px=w1+k11*1/2; py=w2+k12*1/2; pz=w3+k13*1/2; x=w4+k14*1/2; y=w5+k15*1/2;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y)+py*Efy(x,y)+pz*Efz(x,y); 
 	 Fx=Efx(x,y)+(py*Bfz(x,y)-pz*Bfy(x,y))/gam; 	 
 	 Fy=Efy(x,y)+(pz*Bfx(x,y)-px*Bfz(x,y))/gam; 
 	 Fz=Efz(x,y)+(px*Bfy(x,y)-py*Bfx(x,y))/gam; 
 	  
	 k21=h*fun1(px,py,pz,x,y); 
 	 k22=h*fun2(px,py,pz,x,y);
 	 k23=h*fun3(px,py,pz,x,y);
 	 k24=h*fun4(px,py,pz,x,y);
 	 k25=h*fun5(px,py,pz,x,y);

 	 //-----------------------------------------------------------------------
 	 
 	 px=w1+k21*1/2; py=w2+k22*1/2; pz=w3+k23*1/2; x=w4+k24*1/2; y=w5+k25*1/2;
 	
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y)+py*Efy(x,y)+pz*Efz(x,y); 
 	 Fx=Efx(x,y)+(py*Bfz(x,y)-pz*Bfy(x,y))/gam; 	 
 	 Fy=Efy(x,y)+(pz*Bfx(x,y)-px*Bfz(x,y))/gam; 
 	 Fz=Efz(x,y)+(px*Bfy(x,y)-py*Bfx(x,y))/gam; 
 	 
 	 k31=h*fun1(px,py,pz,x,y); 
 	 k32=h*fun2(px,py,pz,x,y);
 	 k33=h*fun3(px,py,pz,x,y);
 	 k34=h*fun4(px,py,pz,x,y);
 	 k35=h*fun5(px,py,pz,x,y);

 	 //---------------------------------------------------------------------------
 	  
 	 t=i*h;
 	 
 	 px=w1+k31; py=w2+k32; pz=w3+k33; x=w4+k34; y=w5+k35;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y)+py*Efy(x,y)+pz*Efz(x,y); 
 	 Fx=Efx(x,y)+(py*Bfz(x,y)-pz*Bfy(x,y))/gam; 	 
 	 Fy=Efy(x,y)+(pz*Bfx(x,y)-px*Bfz(x,y))/gam; 
 	 Fz=Efz(x,y)+(px*Bfy(x,y)-py*Bfx(x,y))/gam; 
 	 
 	 k41=h*fun1(px,py,pz,x,y);  
 	 k42=h*fun2(px,py,pz,x,y);
 	 k43=h*fun3(px,py,pz,x,y);
 	 k44=h*fun4(px,py,pz,x,y);
 	 k45=h*fun5(px,py,pz,x,y);

 	 
 	 w1=w1+(k11+2*k21+2*k31+k41)/6; 
 	 w2=w2+(k12+2*k22+2*k32+k42)/6; 
 	 w3=w3+(k13+2*k23+2*k33+k43)/6; 
 	 w4=w4+(k14+2*k24+2*k34+k44)/6; 
 	 w5=w5+(k15+2*k25+2*k35+k45)/6;
 	 
 
 	 //---------------------------------------------------------------------------
 	 
 	 if ((i % pri) == 0)
 	 {gam=sqrt(1+w1*w1+w2*w2+w3*w3); 
 	  fprintf(fo, "%.14e    	%.14e 	 %.14e  	%.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w1, w2, w3, gam); 
 	  // prints t, x, y, px, py, pz, gamma
	  }
	 }
 	}
 	fclose(fo);
}

int main()
{ double p01, p02, p03, x01, x02;
  double  T; 
 long long int N; 
 FILE *foo;
 
 foo=fopen("InputTotBatch.txt","r");
 fscanf(foo,"%lf %lf %lf %lf %lf %lf %lld %i %lf %lf ", &p01, &p02, &Eo,&kdamp,&k, &T, &N, &pri, &xgrid, &ygrid);
 tfwhm=40.0;
 
 w=k; Bo=Eo; 
 fclose(foo);
 p03=0.0; x01=0.0; x02=0.0; 
  dx=xgrid/2.0;dy=ygrid/2.0;
 

 RK(T, N, p01, p02, p03, x01, x02);

 printf("gotovo\n");


 }
 
 /*
 
 gcc -O2 trial.c -o trial
./trial  
*/