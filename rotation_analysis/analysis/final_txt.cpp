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

	string directory; // INCLUDE LAST SLASH
	string file_out;



	for(int i=1; i<=6; i++){
		for(int j=1; j<=3; j++){
			double Eo = (double)(5*i);
			double px0 = -(double)(15+j*5);
			char aux1[128];
			char aux2[128];
			sprintf(aux1, "../outputs/Data01_a0_%02.lf_p0_%02.lf/", Eo, -px0);
			sprintf(aux2, "Data_px%02.lf_a%02.lf.txt", -px0, Eo);
			directory = string(aux1);
			file_out = string(aux2);
			cout << directory << endl;
			get_final_txt(directory, file_out);
		}
	}

	return 0;

}