#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void get_final_txt(string directory, string file_out){

	char aux[512];
	int n_files;
	sprintf(aux, "%splots/", directory.c_str());
	string* files = list_dir(aux, &n_files);

	sprintf(aux, "%s%s", directory.c_str(), file_out.c_str());
	FILE* fout = fopen(aux, "w");
	fprintf(fout, "y\tz\tphi\tr\tperiod\tthetamin\tthetamax\traio_max\tpx_mean");

	printf("WRITING TO FILE %s\n", aux);

	FILE* fin;
	int file_to_read;
	for(int i = 0; i < n_files; i++){
		sprintf(aux, "%splots/plotOut%%05d.txt", directory.c_str());
		sscanf(files[i].c_str(), aux, &file_to_read);
		sprintf(aux, "%slogs/logOut%05d.txt", directory.c_str(), file_to_read);
		fin = fopen(aux, "r");
		fgets(aux, 512, fin);
		fprintf(fout, "\n%s", aux);
		fclose(fin);
	}
	fclose(fout);

	printf("FINISHED\n");

}

int main(int argc, char **argv){
	
	/*
	CHANGE THE MAIN DIRECTORY OF THE PLOTS AND THE FILE NAME
	*/

	string directory = "../outputs/Data_Teste_dividir_gamma0_v2/"; // INCLUDE LAST SLASH
	string file_out = "Data_Out.txt";


	directory = "../outputs/Data_Stopped_00.5/";
	get_final_txt(directory, file_out);	
	directory = "../outputs/Data_Stopped_01.5/";
	get_final_txt(directory, file_out);	
	directory = "../outputs/Data_Stopped_02.5/";
	get_final_txt(directory, file_out);	
	directory = "../outputs/Data_Stopped_03.5/";
	get_final_txt(directory, file_out);	
	directory = "../outputs/Data_Stopped_04.5/";
	get_final_txt(directory, file_out);	


	directory = "../outputs/Data_Stopped_04.0/";
	get_final_txt(directory, file_out);	
	directory = "../outputs/Data_Stopped_05.0/";
	get_final_txt(directory, file_out);	
	


	//directory = "../outputs/Data_Luis/";
	//get_final_txt(directory, file_out);	

	return 0;

}