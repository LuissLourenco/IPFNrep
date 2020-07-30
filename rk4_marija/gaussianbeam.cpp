
double w0 = 1;
double lambda = 1;
double n = 1;
double E0 = 1;

double zr = M_PI * w0 * w0 * n / lambda;
double k = 2 * M_PI * n / lambda;

double w(double x){ return w0*sqrt(1 + x*x / (zr*zr)); }
double R(double x){ return x * (1 + zr*zr / (x*x)); }
double psi(double x){ return atan(x/zr); }

double myfunc1(double* x, double* par){
	double z = x[0];
	double r = x[1];
	return E0 * (w0/w(z)) * exp(-r*r / (w(z)*w(z))) * cos(k*z + k*r*r/(2*R(z)) - psi(z));
}

void gaussianbeam(){

	TF2* f1 = new TF2("f1", myfunc1, -10, 10, -10, 10, 0);
	f1->SetTitle(";x;t;Gaussian Beam");
	gStyle->SetPalette(kMint);

	TCanvas* c1 = new TCanvas("c", "", 1200, 600);

	c1->cd();
	f1->Draw("COLZ");

	c1->SaveAs("GaussianBeam.png");

}