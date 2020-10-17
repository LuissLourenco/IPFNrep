#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void analise_folder(string directory, int n_terminals){

	system("g++ -o2 single_analise.cpp -lm -o single_analise `root-config --cflags --glibs`");

	int n_files;
	string *files;
	char plot[64], log[64], cmd[516];

	sprintf(cmd, "rm %splots/*", directory.c_str());
	system(cmd);
	sprintf(cmd, "rm %slogs/*", directory.c_str());
	system(cmd);

	files = list_dir(directory.c_str(), &n_files);

	sprintf(cmd, "mkdir %splots", directory.c_str());
	system(cmd);
	sprintf(cmd, "mkdir %slogs", directory.c_str());
	system(cmd);

	for(int i = 0; i < n_files; i++){
		sprintf(plot, "%splots/plot%s.png", directory.c_str(), files[i].substr(string(directory).size(), 8).c_str());
		sprintf(log , "%slogs/log%s.txt"  , directory.c_str(), files[i].substr(string(directory).size(), 8).c_str());
		if(i % n_terminals == 0){
			sprintf(cmd, "./single_analise %s %s %s; sleep 1", files[i].c_str(), plot, log);
			system(cmd);
		}else{
			sprintf(cmd, "gnome-terminal --tab -- bash -ic './single_analise %s %s %s'", files[i].c_str(), plot, log);
			system(cmd);
		}
	}

}

int main(int argc, char** argv){

	/*
	CHANGE THE DIRECTORY OF THE SIMULATIONS HERE
	n_terminals IS THE NUMBER OS TABS THAT WILL BE RAN AT THE SAME TIME
	AFTER FINISHED RUNNING, CHECK EVERY FILE AT directory/plots/ AND DELETE
		THE ONES THAT ARE TRASH
	*/

	string directory = "../outputs/Data_Teste_dividir_gamma0_v3/"; // INCLUDE LAST SLASH
	int n_terminals = 6;

	directory = "../outputs/PicardLindelof02/";
	analise_folder(directory, n_terminals);
	/*directory = "../outputs/Data_Stopped_01.5/";
	analise_folder(directory, n_terminals);
	directory = "../outputs/Data_Stopped_02.5/";
	analise_folder(directory, n_terminals);
	directory = "../outputs/Data_Stopped_03.5/";
	analise_folder(directory, n_terminals);
	directory = "../outputs/Data_Stopped_04.5/";
	analise_folder(directory, n_terminals);

*/
	//directory = "../outputs/Data_Luis/";
	//analise_folder(directory, n_terminals);



	/*

	for(int i=1; i<=6; i++){
		for(int j=1; j<=3; j++){

			double Eo = 5.*i;
			double px0 = -15. - j*5.;
			char aux[128];
			sprintf(aux, "../outputs/PicardLindelof_a0_%02.0lf_px0_%02.0lf/", Eo, px0);
			string directory = aux;
			analise_folder(directory, 6);
		}
	}
	


	*/
	return 0;

}