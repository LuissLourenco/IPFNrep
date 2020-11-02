#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void analise_folder(string directory, int n_terminals){

	system("g++ -o2 single_analise2.cpp -lm -o single_analise2 `root-config --cflags --glibs`");

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
			sprintf(cmd, "./single_analise2 %s %s %s; sleep 1", files[i].c_str(), plot, log);
			system(cmd);
		}else{
			sprintf(cmd, "gnome-terminal --tab -- bash -ic './single_analise2 %s %s %s'", files[i].c_str(), plot, log);
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

	string directory = "../outputs/Folder/"; // INCLUDE LAST SLASH
	int n_terminals = 4;
	//analise_folder(directory, n_terminals);

	/*
	for(int i=1; i<=2; i++){
		for(int j=3; j<=3; j++){
			double Eo = (double)(5*i);
			double px0 = (double)(15+j*5);
			char aux[128];
			sprintf(aux, "../outputs/Data01_a0_%02.lf_p0_%02.lf/", Eo, px0);
			directory = aux;
			analise_folder(directory,n_terminals);
		}
	}*/

	directory = "../outputs/Teste_cool/";
	//analise_folder(directory,4);
	directory = "../outputs/Teste_pl02/";
	//analise_folder(directory,4);
	directory = "../outputs/Teste_pl0-1/";

	directory = "../outputs/teste1/";
	//analise_folder(directory,1);
	directory = "../outputs/teste2/";
	//analise_folder(directory,1);
	directory = "../outputs/teste3/";
	//analise_folder(directory,1);

	directory = "../outputs/Data01_a0_20_p0_20/";
	analise_folder(directory,1);
	directory = "../outputs/Data0-1_a0_20_p0_20/";
	analise_folder(directory,1);
	
	return 0;

}