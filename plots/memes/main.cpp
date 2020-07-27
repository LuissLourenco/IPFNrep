#include <iostream>
#include <chrono> 
using namespace std::chrono; 
#include <vector>
#include <string>
#include <cmath>
#include <ODEpoint.h>

#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;



vector<ODEpoint> solucao(double a0, double delta, double T, double dt, string frame){

	vector<ODEpoint> sol;
	
	double* p = new double[6];
	for(int i=0; i<6; i++) p[i] = 0;

	double t = 0; 
	double phi = 0;

	if(frame == "lab"){
		auto start = high_resolution_clock::now(); 
		while(t<T){


			//f(phi) = 0 -> Método de Newton
			phi = t-p[0];
			for(int i=0; i<T; i++){
				double funcao = 0.25*a0*a0*(phi+(2*delta*delta-1)/2*sin(2*phi))+phi-t;
				double derivada = 0.25*a0*a0*(1+(2*delta*delta-1)*cos(2*phi))+1;
				double phi1 = phi - funcao/derivada;

				if(abs(phi1-phi)<1e-10){
					phi = phi1;
					break;
				}

				if(i==9){
					cout << " *** método de Newton não convergiu *** " <<endl;  
					cout << " *** abs(phi1-phi) = " << abs(phi1-phi) << " *** " <<endl;
				}

				phi = phi1;
			}

			p[0] = 0.25*a0*a0*(phi+(2*delta*delta-1)/2*sin(2*phi));
			p[1] = delta * a0 * sin(phi);
			p[2] = -sqrt(1-delta*delta) * a0 * cos(phi);
			p[3] = 0.25*a0*a0*(1+(2*delta*delta-1)*cos(2*phi));
			p[4] = delta*a0*cos(phi);
			p[5] = sqrt(1-delta*delta)*a0*sin(phi);


			sol.push_back(ODEpoint(t, p, 6));
			t += dt;	
		}
		auto stop = high_resolution_clock::now(); 
		auto duration = duration_cast<milliseconds>(stop - start); 
		cout << " *** " << duration.count() << " ms ***" << endl; 
	}

	else if(frame == "rest"){

		double gama0 = sqrt(1+a0*a0/2);
		double q = a0/2/gama0;

		auto start = high_resolution_clock::now(); 
		while(t<T){


			//f(phi) = 0 -> Método de Newton
			phi = t-p[0];
			for(int i=0; i<10; i++){
				double funcao = (2*delta*delta-1)/2*q*q*sin(2*phi)+phi-t;
				double derivada = (2*delta*delta-1)*q*q*cos(2*phi)+1;
				double phi1 = phi - funcao/derivada;

				if(abs(phi1-phi)<1e-10){
					phi = phi1;
					break;
				}

				if(i==9){
					cout << " *** método de Newton não convergiu *** " <<endl;  
					cout << " *** abs(phi1-phi) = " << abs(phi1-phi) << " *** " <<endl;
				}

				phi = phi1;
			}

			p[0] = (2*delta*delta-1)/2*q*q*sin(2*phi);
			p[1] = 2*delta*q*sin(phi);
			p[2] = -2*sqrt(1-delta*delta)*q*cos(phi);
			p[3] = (2*delta*delta-1)*a0*a0/4/gama0*cos(2*phi);
			p[4] = delta*a0*cos(phi);
			p[5] = sqrt(1-delta*delta)*a0*sin(phi);


			sol.push_back(ODEpoint(t, p, 6));
			t += dt;	
		}
		auto stop = high_resolution_clock::now(); 
		auto duration = duration_cast<milliseconds>(stop - start); 
		cout << " *** " << duration.count() << " ms ***" << endl; 
	}


	delete[] p;

	return sol;
};



int main(int argc, char** argv){


	auto sol1 = solucao(1, 1, 20, 0.0001, "lab");
	

	int n1 = sol1.size();
	TGraph** g1 = new TGraph*[8];
	for(int i=0; i<8; i++)
		g1[i] = new TGraph(n1);

	for(int i=0; i<n1; i++){
		g1[0]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[0]);
		g1[1]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[1]);
		g1[2]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[2]);
		g1[3]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[3]);
		g1[4]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[4]);
		g1[5]->SetPoint(i, sol1[i].Get()[6], sol1[i].Get()[5]);
		g1[6]->SetPoint(i, sol1[i].Get()[0], sqrt(sol1[i].Get()[1]*sol1[i].Get()[1]+sol1[i].Get()[2]*sol1[i].Get()[2]));
		g1[7]->SetPoint(i, sol1[i].Get()[3], sqrt(sol1[i].Get()[4]*sol1[i].Get()[4]+sol1[i].Get()[5]*sol1[i].Get()[5]));
	}

	TMultiGraph* mg11 = new TMultiGraph();
	TMultiGraph* mg12 = new TMultiGraph();
	mg11->Add(g1[1]);
	mg11->Add(g1[2]);
	mg12->Add(g1[4]);
	mg12->Add(g1[5]);

	g1[0]->SetLineColor(kBlue);
	g1[1]->SetLineColor(kBlue);
	g1[2]->SetLineColor(kViolet);
	g1[3]->SetLineColor(kBlue);
	g1[4]->SetLineColor(kBlue);
	g1[5]->SetLineColor(kViolet);
	g1[6]->SetLineColor(kBlue);
	g1[7]->SetLineColor(kBlue);

	g1[0]->SetTitle(";t;x");
	g1[1]->SetTitle("y");
	g1[2]->SetTitle("z");
	g1[3]->SetTitle(";t;px");
	g1[4]->SetTitle("py");
	g1[5]->SetTitle("pz");
	g1[6]->SetTitle(";r_{long};r_{perp}");
	g1[7]->SetTitle(";p_{long};p_{perp}");
	mg11->SetTitle(";t");
	mg12->SetTitle(";t");

	TCanvas* c1 = new TCanvas("c1", "", 1000, 650);
	c1->Divide(3, 2);

	c1->cd(1);
	g1[0]->Draw("AL");
	c1->cd(4);
	mg11->Draw("AL");
	c1->cd(2);
	g1[3]->Draw("AL");
	c1->cd(5);
	mg12->Draw("AL");
	c1->cd(3);
	g1[6]->Draw("AL");
	c1->cd(6);
	g1[7]->Draw("AL");

	TLegend* leg11 = (c1->cd(4))->BuildLegend(0.7, 0.15, 0.9, 0.3);
	leg11->SetBorderSize(0);
	leg11->SetFillStyle(0);

	TLegend* leg12 = (c1->cd(5))->BuildLegend(0.7, 0.15, 0.9, 0.3);
	leg12->SetBorderSize(0);
	leg12->SetFillStyle(0);

	c1->SaveAs("lab_PolLin_1.png");

	









	for(int i=0; i<8; i++)
		delete g1[i];
	delete[] g1;

	delete c1;
	delete mg11;
	delete mg12;

	return 0;
}
