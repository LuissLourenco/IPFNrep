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

	sim1.Eo = 8;
	sim1.px0 = -1;
	sim1.tfwhm = 50;

	sim1.T = 300;
	sim1.N = sim1.T / 0.005;
	sim1.directory = "";
	
	//sim1.print();

	run_simulations(sim1, 6); // 500
	//run_simulations(sim2, 6); // 5

	sim1.print();
/*
	int rise_min = 10;
	int rise_max = 1500;
	int drise = 10;
	char file[256];
	for(int rise = rise_min; rise <= rise_max; rise+=drise){
*/
		/*
		cout << "SIMULATING TFWHM = " << rise << endl;

		sim1.tfwhm = rise;
		run_simulations(sim1, 6);
		analise("Out00000.txt", "Plot00000.png", "Log00000.txt", (double)rise/sim1.T);
		cout << "---------------" << (double)rise/sim1.T <<endl;


		system("sleep 0.001;");
		sprintf(file, "mv Out00000.txt tfwhm/Out_tfwhm%03d.txt", rise);
		system(file);
		sprintf(file, "mv Plot00000.png tfwhm/Plot_tfwhm%03d.png", rise);
		system(file);
		sprintf(file, "mv Log00000.txt tfwhm/Log_tfwhm%03d.txt", rise);
		system(file);
		system("sleep 0.001;");
		*/
	/*	

		
		char f1[256], f2[256], f3[256];
		sprintf(f1, "tfwhm/Out_tfwhm%03d.txt", rise);
		sprintf(f2, "tfwhm/Plot_tfwhm%03d.png", rise);
		sprintf(f3, "tfwhm/Log_tfwhm%03d.txt", rise);
		analise(f1, f2, f3);

	}

	char aux[512];
	vector<string> files = { 

"tfwhm/Log_tfwhm010.txt",
"tfwhm/Log_tfwhm020.txt",
"tfwhm/Log_tfwhm030.txt",
"tfwhm/Log_tfwhm040.txt",
"tfwhm/Log_tfwhm050.txt",
"tfwhm/Log_tfwhm060.txt",
"tfwhm/Log_tfwhm070.txt",
"tfwhm/Log_tfwhm080.txt",
"tfwhm/Log_tfwhm090.txt",
"tfwhm/Log_tfwhm100.txt",
"tfwhm/Log_tfwhm110.txt",
"tfwhm/Log_tfwhm120.txt",
"tfwhm/Log_tfwhm130.txt",
"tfwhm/Log_tfwhm140.txt",
"tfwhm/Log_tfwhm150.txt",
"tfwhm/Log_tfwhm160.txt",
"tfwhm/Log_tfwhm170.txt",
"tfwhm/Log_tfwhm180.txt",
"tfwhm/Log_tfwhm190.txt",
"tfwhm/Log_tfwhm200.txt",
"tfwhm/Log_tfwhm210.txt",
"tfwhm/Log_tfwhm220.txt",
"tfwhm/Log_tfwhm230.txt",
"tfwhm/Log_tfwhm240.txt",
"tfwhm/Log_tfwhm250.txt",
"tfwhm/Log_tfwhm260.txt",
"tfwhm/Log_tfwhm270.txt",
"tfwhm/Log_tfwhm280.txt",
"tfwhm/Log_tfwhm290.txt",
"tfwhm/Log_tfwhm300.txt",
"tfwhm/Log_tfwhm310.txt",
"tfwhm/Log_tfwhm320.txt",
"tfwhm/Log_tfwhm330.txt",
"tfwhm/Log_tfwhm340.txt",
"tfwhm/Log_tfwhm350.txt",
"tfwhm/Log_tfwhm360.txt",
"tfwhm/Log_tfwhm370.txt",
"tfwhm/Log_tfwhm380.txt",
"tfwhm/Log_tfwhm390.txt",
"tfwhm/Log_tfwhm400.txt",
"tfwhm/Log_tfwhm410.txt",
"tfwhm/Log_tfwhm420.txt",
"tfwhm/Log_tfwhm430.txt",
"tfwhm/Log_tfwhm440.txt",
"tfwhm/Log_tfwhm450.txt",
"tfwhm/Log_tfwhm460.txt",
"tfwhm/Log_tfwhm470.txt",
"tfwhm/Log_tfwhm480.txt",
"tfwhm/Log_tfwhm490.txt",
"tfwhm/Log_tfwhm500.txt",
"tfwhm/Log_tfwhm510.txt",
"tfwhm/Log_tfwhm520.txt",
"tfwhm/Log_tfwhm530.txt",
"tfwhm/Log_tfwhm540.txt",
"tfwhm/Log_tfwhm550.txt",
"tfwhm/Log_tfwhm560.txt",
"tfwhm/Log_tfwhm570.txt",
"tfwhm/Log_tfwhm580.txt",
"tfwhm/Log_tfwhm590.txt",
"tfwhm/Log_tfwhm600.txt",
"tfwhm/Log_tfwhm610.txt",
"tfwhm/Log_tfwhm620.txt",
"tfwhm/Log_tfwhm630.txt",
"tfwhm/Log_tfwhm640.txt",
"tfwhm/Log_tfwhm650.txt",
"tfwhm/Log_tfwhm660.txt",
"tfwhm/Log_tfwhm670.txt",
"tfwhm/Log_tfwhm680.txt",
"tfwhm/Log_tfwhm690.txt",
"tfwhm/Log_tfwhm700.txt",
"tfwhm/Log_tfwhm710.txt",
"tfwhm/Log_tfwhm720.txt",
"tfwhm/Log_tfwhm730.txt",
"tfwhm/Log_tfwhm740.txt",
"tfwhm/Log_tfwhm750.txt",
"tfwhm/Log_tfwhm760.txt",
"tfwhm/Log_tfwhm770.txt",
"tfwhm/Log_tfwhm780.txt",
"tfwhm/Log_tfwhm790.txt",
"tfwhm/Log_tfwhm800.txt",
"tfwhm/Log_tfwhm810.txt",
"tfwhm/Log_tfwhm820.txt",
"tfwhm/Log_tfwhm830.txt",
"tfwhm/Log_tfwhm840.txt",
"tfwhm/Log_tfwhm850.txt",
"tfwhm/Log_tfwhm860.txt",
"tfwhm/Log_tfwhm870.txt",
"tfwhm/Log_tfwhm880.txt",
"tfwhm/Log_tfwhm890.txt",
"tfwhm/Log_tfwhm900.txt",
"tfwhm/Log_tfwhm910.txt",
"tfwhm/Log_tfwhm920.txt",
"tfwhm/Log_tfwhm930.txt",
"tfwhm/Log_tfwhm940.txt",
"tfwhm/Log_tfwhm950.txt",
"tfwhm/Log_tfwhm960.txt",
"tfwhm/Log_tfwhm970.txt",
"tfwhm/Log_tfwhm980.txt",
"tfwhm/Log_tfwhm990.txt",
"tfwhm/Log_tfwhm1000.txt",
"tfwhm/Log_tfwhm1010.txt",
"tfwhm/Log_tfwhm1020.txt",
"tfwhm/Log_tfwhm1030.txt",
"tfwhm/Log_tfwhm1040.txt",
"tfwhm/Log_tfwhm1050.txt",
"tfwhm/Log_tfwhm1060.txt",
"tfwhm/Log_tfwhm1070.txt",
"tfwhm/Log_tfwhm1080.txt",
"tfwhm/Log_tfwhm1090.txt",
"tfwhm/Log_tfwhm1100.txt",
"tfwhm/Log_tfwhm1110.txt",
"tfwhm/Log_tfwhm1120.txt",
"tfwhm/Log_tfwhm1130.txt",
"tfwhm/Log_tfwhm1140.txt",
"tfwhm/Log_tfwhm1150.txt",
"tfwhm/Log_tfwhm1160.txt",
"tfwhm/Log_tfwhm1170.txt",
"tfwhm/Log_tfwhm1180.txt",
"tfwhm/Log_tfwhm1190.txt",
"tfwhm/Log_tfwhm1200.txt",
"tfwhm/Log_tfwhm1210.txt",
"tfwhm/Log_tfwhm1220.txt",
"tfwhm/Log_tfwhm1230.txt",
"tfwhm/Log_tfwhm1240.txt",
"tfwhm/Log_tfwhm1250.txt",
"tfwhm/Log_tfwhm1260.txt",
"tfwhm/Log_tfwhm1270.txt",
"tfwhm/Log_tfwhm1280.txt",
"tfwhm/Log_tfwhm1290.txt",
"tfwhm/Log_tfwhm1300.txt",
"tfwhm/Log_tfwhm1310.txt",
"tfwhm/Log_tfwhm1320.txt",
"tfwhm/Log_tfwhm1330.txt",
"tfwhm/Log_tfwhm1340.txt",
"tfwhm/Log_tfwhm1350.txt",
"tfwhm/Log_tfwhm1360.txt",
"tfwhm/Log_tfwhm1370.txt",
"tfwhm/Log_tfwhm1380.txt",
"tfwhm/Log_tfwhm1390.txt",
"tfwhm/Log_tfwhm1400.txt",
"tfwhm/Log_tfwhm1410.txt",
"tfwhm/Log_tfwhm1420.txt",
"tfwhm/Log_tfwhm1430.txt",
"tfwhm/Log_tfwhm1440.txt",
"tfwhm/Log_tfwhm1450.txt",
"tfwhm/Log_tfwhm1460.txt",
"tfwhm/Log_tfwhm1470.txt",
"tfwhm/Log_tfwhm1480.txt",
"tfwhm/Log_tfwhm1490.txt",
"tfwhm/Log_tfwhm1500.txt"

	};

	FILE* fout = fopen("Data_tfwhm.txt", "w");
	fprintf(fout, "tfwhm\tphi\tr\tperiod\traio_max\tpx_mean\tx_f\teta");

	FILE* fin;
	int rise_print;
	for(int i = 0; i < files.size(); i++){
		fin = fopen(files[i].c_str(), "r");
		fgets(aux, 512, fin);
		aux[0] = 48;
		sscanf(files[i].c_str(), "tfwhm/Log_tfwhm%d", &rise_print); 
		fprintf(fout, "\n%d\t%s", rise_print, aux);
		fclose(fin);
	}
	fclose(fout);
*/



	//analise("Out00000.txt", "Plot00000.png", "Log00000.txt");
	
	TApplication* theApp = new TApplication("App", &argc, argv);

	int n_points, n_cols;
	double** values = ReadFile("Out00000.txt", &n_cols, &n_points, true);
	//double** values = ReadFile("Sim_a0.1_p0.txt", &n_cols, &n_points, true);

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