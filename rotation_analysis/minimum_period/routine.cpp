#include "../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"
#include "TApplication.h"

#include "../simulate/simulate_helper.cpp"
#include "single_analise.cpp"

using namespace std;

#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TFrame.h"
#include "TPad.h"
#include "TAxis3D.h"
#include "TView3D.h"

using namespace std;

class Axis3D{

public:
	Axis3D(double a, double b, double c, double d, double e, double f){
		
		view = (TView3D*) TView::CreateView(1);

		axis = new TAxis3D();
		axis->SetAxisColor(1);
		axis->SetLabelColor(1);

		SetRange(a, b, c, d, e, f);

	};
	~Axis3D(){};

	//==============//

	void SetRange(double a, double b, double c, double d, double e, double f){
		x_min = a;y_min = b;z_min = c;x_max = d;y_max = e;z_max = f;
		view->SetRange(a, b, c, d, e, f);
	};
	void RotateView(double psi, double theta){
		view->RotateView(psi, theta);
	};
	void SetGrid(double xlow_a, double dx_a, double ylow_a, double dy_a, double zlow_a, double dz_a){
		xlow = xlow_a; dx = dx_a; ylow = ylow_a; dy = dy_a;	zlow = zlow_a; dz = dz_a;
		draw_grid = true;
	};
	TAxis3D* GetAxis(){
		return axis;
	};
	void Draw(){
		axis->Draw();
		outter_cube();
		if(draw_grid) Grid();
	};

private:
	double x_min, y_min, z_min, x_max, y_max, z_max;
	double xlow=0, dx=0, ylow=0, dy=0, zlow=0, dz=0;
	TView3D *view;
	TAxis3D* axis;
	bool draw_grid = false;

	void outter_cube(){
		double** points = new double*[3];
		points[0] = new double[2]; points[1] = new double[2]; points[2] = new double[2];

		double A[3] = {x_min, y_min, z_min};
		double B[3] = {x_min, y_max, z_min};
		double C[3] = {x_max, y_max, z_min};
		double D[3] = {x_max, y_min, z_min};
		double E[3] = {x_min, y_min, z_max};
		double F[3] = {x_min, y_max, z_max};
		double G[3] = {x_max, y_max, z_max};
		double H[3] = {x_max, y_min, z_max};

		TPolyLine3D** outline = new TPolyLine3D*[12];

		auto get_line = [] (double* X1, double* X2){
			double** res = new double*[3];
			res[0] = new double[2]; res[1] = new double[2]; res[2] = new double[2];
			res[0][0] = X1[0];
			res[0][1] = X2[0];
			res[1][0] = X1[1];
			res[1][1] = X2[1];
			res[2][0] = X1[2];
			res[2][1] = X2[2];
			return res;
		};

		points = get_line(A, B); outline[0] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(B, C); outline[1] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(C, D); outline[2] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(D, A); outline[3] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(E, F); outline[4] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(F, G); outline[5] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(G, H); outline[6] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(H, E); outline[7] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(A, E); outline[8] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(B, F); outline[9] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(C, G); outline[10] = new TPolyLine3D(2, points[0], points[1], points[2]);
		points = get_line(D, H); outline[11] = new TPolyLine3D(2, points[0], points[1], points[2]);
	
		
		for(int i = 0; i < 12; i++)
			outline[i]->Draw();

	};
	void Grid(){
		vector<TPolyLine3D*> grid;
		auto get_line = [] (double* X1, double* X2){
			double** res = new double*[3];
			res[0] = new double[2]; res[1] = new double[2]; res[2] = new double[2];
			res[0][0] = X1[0];
			res[0][1] = X2[0];
			res[1][0] = X1[1];
			res[1][1] = X2[1];
			res[2][0] = X1[2];
			res[2][1] = X2[2];
			return res;
		};
		double aux;
		aux = xlow;
		double **P;
		while(aux < x_max){
			double P1[3] = {aux, y_min, z_min};
			double P2[3] = {aux, y_max, z_min};
			P = get_line(P1, P2);
			grid.push_back(new TPolyLine3D(2, P[0], P[1], P[2]));
			aux += dx;
		}

		aux = ylow;
		while(aux < y_max){
			double P1[3] = {x_min, aux, z_min};
			double P2[3] = {x_max, aux, z_min};
			P = get_line(P1, P2);
			grid.push_back(new TPolyLine3D(2, P[0], P[1], P[2]));
			aux += dy;
		}

		aux = zlow;
		while(aux < z_max){
			double P1[3] = {x_min, y_min, aux};
			double P2[3] = {x_min, y_max, aux};
			P = get_line(P1, P2);
			grid.push_back(new TPolyLine3D(2, P[0], P[1], P[2]));
			double P3[3] = {x_min, y_max, aux};
			double P4[3] = {x_max, y_max, aux};
			P = get_line(P3, P4);
			grid.push_back(new TPolyLine3D(2, P[0], P[1], P[2]));
			aux += dz;
		}

		for(int i = 0; i < grid.size(); i++)
			grid[i]->Draw();

	};

};




int main(int argc, char** argv){

	simulation_data sim1;
	
	sim1.kdamp = 0;
	sim1.pri = 10;
	sim1.wave_type = 3;
	sim1.stable = 1000000;//2E300;
	sim1.delta = 0;
	sim1.w0 = 5;
	sim1.lambda = 1;
	sim1.l = 1;
	sim1.p = 0;

	sim1.rmin = 2;
	sim1.rmax = 2.1;
	sim1.dr = 1;
	sim1.phimin = 0;
	sim1.phimax = 1;
	sim1.dphi = 30;

	sim1.pri = 10;

	sim1.Eo = 5;
	sim1.px0 = -1;
	sim1.tfwhm = 500;

	sim1.T = 800;
	sim1.N = sim1.T / 0.005;
	sim1.directory = "";

	simulation_data sim2 = sim1;
	sim2.tfwhm = 5;
	
	//sim1.print();

	//run_simulations(sim1, 6); // 500
	run_simulations(sim2, 6); // 5

	//analise("Out00000.txt", "Plot00000.png", "Log00000.txt");
	
	TApplication* theApp = new TApplication("App", &argc, argv);

	int n_points, n_cols;
	//double** values = ReadFile("Out00000.txt", &n_cols, &n_points, true);
	double** values = ReadFile("Out_a0.1_p0.txt", &n_cols, &n_points, true);

	cout << "N_POINTS = " << n_points << endl;
	int start = 0;
	int end = n_points-1;

	DataSet T = DataSet(n_points, values[0]).subDataSet(start, end);
	DataSet X = DataSet(n_points, values[1]).subDataSet(start, end);
	DataSet Y = DataSet(n_points, values[2]).subDataSet(start, end);
	DataSet Z = DataSet(n_points, values[3]).subDataSet(start, end);

	TGraph* g = GetTGraph(T, Y);

	TPolyLine3D* traj = new TPolyLine3D(X.size(), X.array(), Y.array(), Z.array());
	traj->SetLineColor(4);
	traj->SetLineWidth(2);


	TCanvas* c1 = new TCanvas("c", "", 1500, 1500);
	//c1->SetRightMargin(0.0);
	//c1->SetLeftMargin(0.0);
	//c1->SetBottomMargin(0.0);
	//c1->SetTopMargin(0.0);

	c1->SetGrid();
	c1->cd();

   	Axis3D* axis = new Axis3D(X.getMin().val(), Y.getMin().val(), Z.getMin().val(), X.getMax().val(), Y.getMax().val(), Z.getMax().val());
	axis->RotateView(-30, 80);
	//axis->SetGrid(0, 1, -10, 2, -10, 2);

	axis->GetAxis()->GetXaxis()->SetNdivisions(5);
	axis->GetAxis()->SetXTitle("x");
	axis->GetAxis()->SetYTitle("y");
	axis->GetAxis()->SetZTitle("z");
	axis->GetAxis()->SetTitleOffset(1.5);

	c1->Divide(2, 1);
	c1->cd(1); axis->Draw(); traj->Draw();
	c1->cd(2); g->Draw("APL");


	theApp->Run();


	return 0;
}