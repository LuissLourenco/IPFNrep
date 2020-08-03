
double tfwhm=10.0;
double stable=10.0;

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

double w0 = 10;
double lambda = 1;
double n = 1;
double E0 = 1;
double omega = 1;

double zr = M_PI * w0 * w0 * n / lambda;
double k = 2 * M_PI * n / lambda;

double w(double x){ return w0*sqrt(1 + x*x / (zr*zr)); }
double R(double x){ return x * (1 + zr*zr / (x*x)); }
double psi(double x){ return atan(x/zr); }

double myfunc1(double* x, double* par){
	double z = x[0];
	double r = x[1];
	return E0 * (w0/w(z)) * exp(-r*r / (w(z)*w(z))) * cos(omega * par[0] - k*z - k*r*r/(2*R(z)) + psi(z)) * Envelope(z, par[0]);
}

void gaussianbeam(){

	TF2* f1 = new TF2("f1", myfunc1, -5, 5, -20, 20, 1); 
	f1->SetTitle(";z;r;Gaussian Beam");
	///gStyle->SetPalette(kMint);
	f1->SetNpx(100);
	f1->SetNpy(100);

	f1->SetParameter(0, 15); //SET TIME

	TCanvas* c1 = new TCanvas("c", "", 1200, 600);

	c1->cd();

	f1->Draw("SURF2"); 
	c1->SaveAs("GaussianBeam.png");

}