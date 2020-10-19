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
	fprintf(fout, "phi\tr\tperiod\traio_max\tpx_mean");

	printf("WRITING TO FILE %s\n", aux);

	FILE* fin;
	int file_to_read;
	for(int i = 0; i < n_files; i++){
		sprintf(aux, "%splots/plotOut%%05d.txt", directory.c_str());
		sscanf(files[i].c_str(), aux, &file_to_read);
		sprintf(aux, "%slogs/logOut%05d.txt", directory.c_str(), file_to_read);
		fin = fopen(aux, "r");
		fgets(aux, 512, fin);
		aux[0] = 48;
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



	directory = "../outputs/Testephi0_a0_5_p0_10/";

	int px_arr[18] = {5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 15, 15, 15, 15, 15, 15};
	int Eo_arr[18] = {5, 10, 15, 20, 25, 30, 5, 10, 15, 20, 25, 30, 5, 10, 15, 20, 25, 30};

	char aux1[128];
	char aux2[128];
	for(int i = 0; i < 18; i++){
		sprintf(aux1, "../outputs/Data_px%02d_a%02d_v2/", px_arr[i], Eo_arr[i]);
		sprintf(aux2, "Data_px%02d_a%02d.txt", px_arr[i], Eo_arr[i]);
	
		directory = aux1;
		file_out = aux2;

		get_final_txt(directory, file_out);	

	}


	return 0;

}