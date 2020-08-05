
#include "Read.cpp"

double tfwhm=50.0;
double stable=40.0;

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
double E0 = 2;
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


void trajetoria_bonita(){

	// t x y z px py pz gama
	vector<vector<double>> sol = ReadFile2Vec("Out0.txt");
	int n = sol.size();

	double xmin=0;
	double xmax=0;
	double ymax=0;

	for(int i=0; i<n; i++){
		if(sol[i][1]<xmin) xmin = sol[i][1];
		if(sol[i][1]>xmax) xmax = sol[i][1];
		if(abs(sol[i][2])>ymax) ymax = abs(sol[i][2]);
	}

	TF2* f1 = new TF2("f1", myfunc1, xmin, xmax, -ymax, ymax, 1); 
	f1->SetTitle("Beam;x;y;");
	gStyle->SetPalette(kBlueYellow);
	f1->SetNpx(300);
	f1->SetNpy(300);

	f1->SetMinimum(-E0);
	f1->SetMaximum(E0);

	TGraph* trajetoria = new TGraph(n);

	trajetoria->SetTitle(";x;y;z");
	trajetoria->SetMarkerStyle(8);
	trajetoria->SetMarkerSize(0.2);
	trajetoria->SetMarkerColor(kRed);


	TCanvas* c1 = new TCanvas("c", "", 1200, 600);

	int frame = 1;
	for(int i=0; i<n; i+=200){

		cout<<" "<<sol[i][0]<<endl;

		trajetoria->SetPoint(i, sol[i][1], sol[i][2]);

		f1->SetParameter(0, sol[i][0]); //SET TIME

		f1->Draw("colz");
		trajetoria->Draw("SAME P");

		//c1->SaveAs((string("saved/frame")+to_string(frame)+string(".png")).c_str());
		//frame++;

		c1->Update();
		c1->Modified();

	}

}