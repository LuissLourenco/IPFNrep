#include "TVectorField.h"

TVectorField::TVectorField(int nr, double* x1r, double* y1r, double* x2r, double* y2r){
	n = nr; 
	x1 = x1r; 
	y1 = y1r; 
	x2 = x2r; 
	y2 = y2r;

	x_min = x1[0];
	x_max = x1[0];
	y_min = y1[0];
	y_max = y1[0];

	for(int i = 0; i < n; i++){
		if(x1[i] < x_min) x_min = x1[i];
		if(x1[i] > x_max) x_max = x1[i];
		if(y1[i] < y_min) y_min = y1[i];
		if(y1[i] > y_max) y_max = y1[i];
	}

	CreateArrows();

}

TVectorField::TVectorField(DataSet x1r, DataSet y1r, DataSet x2r, DataSet y2r){
	if(x1r.size() != x2r.size() || x1r.size() != y1r.size() || x1r.size() != y2r.size()){
		cout << "Tamanhos de Vetores Incompatíveis!" << endl;
		return;
	}
	n = x1r.size();

	x_min = x1r[0].val();
	x_max = x1r[0].val();
	y_min = y1r[0].val();
	y_max = y1r[0].val();

	x1 = new double[n];
	y1 = new double[n];
	x2 = new double[n];
	y2 = new double[n];

	for(int i = 0; i < n; i++){
		x1[i] = x1r[i].val();
		y1[i] = y1r[i].val();
		x2[i] = x2r[i].val();
		y2[i] = y2r[i].val();
		if(x1[i] < x_min) x_min = x1[i];
		if(x1[i] > x_max) x_max = x1[i];
		if(y1[i] < y_min) y_min = y1[i];
		if(y1[i] > y_max) y_max = y1[i];
	}

	CreateArrows();

}


TVectorField::~TVectorField(){
	delete[] x1; 
	delete[] y1; 
	delete[] x2; 
	delete[] y2; 
	delete f;
	for(int i = 0; i < n; i++) delete arr[i];
	delete[] arr;
}

void TVectorField::Draw(){
	CreateArrows();
	f = new TH1D("", "", 100, x_min, x_max);
	f->GetYaxis()->SetRangeUser(y_min, y_max);
	f->SetStats(0);
	f->Draw("AXIS");
	for(int i = 0; i < n; i++) arr[i]->Draw("SAME >");
}

void TVectorField::SetLimits(double v1, double v2, double v3, double v4){
	x_min = v1;
	y_min = v2;
	x_max = v3;
	y_max = v4;
}


void TVectorField::CreateArrows(){
	arr = new TArrow*[n];
	double dx, dy, r;
	int color;
	for(int i = 0; i < n; i++){
		dx = x2[i] - x1[i];
		dy = y2[i] - y1[i];
		r = atan2(dy, dx);
		arr[i] = new TArrow(x1[i], y1[i], x1[i]+cos(r)*scale, y1[i]+sin(r)*scale, 0.005, ">");
		color = (2/(1+exp(-sqrt(dx*dx+dy*dy)*0.5)) - 1) * TColor::GetNumberOfColors();
		arr[i]->SetLineColor(TColor::GetColorPalette(color));
		if(sqrt(dx*dx+dy*dy) <= 1E-10) arr[i]->SetLineWidth(0);
	}
}