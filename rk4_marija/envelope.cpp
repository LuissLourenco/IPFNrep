
double tfwhm=20.0;
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

double myfunc1(double* x, double* par){
	return Envelope(x[0], x[1]);
}

double myfunc2(double* x, double* par){
	return Envelope(x[0], par[0]);
}



void envelope(){


	TF2* f1 = new TF2("f1", myfunc1, -100, 100, -100, 100, 0);
	f1->SetTitle(";x;t;Envelope");
	gStyle->SetPalette(kMint);

	TF1* f2 = new TF1("f2", myfunc2, -100, 100, 1);
	f2->SetParameter(0,0);
	f2->SetTitle(";x;Envelope");

	TArrow *arr1 = new TArrow(-10, 0.9, 30, 0.9, 0.02, ">");
	arr1->SetLineWidth(2);
	arr1->SetLineColor(kRed);


	TCanvas* c1 = new TCanvas("c", "", 1200, 600);
	c1->Divide(2,1);

	c1->cd(1);
	f1->Draw("SURF2");
	c1->cd(2);
	f2->Draw();
	arr1->Draw();

	c1->SaveAs("Envelope.png");

}