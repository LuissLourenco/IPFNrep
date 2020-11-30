#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void get_final_txt(string directory, string file_out, int n_files){

	char aux[512];
	
	sprintf(aux, "%s%s", directory.c_str(), file_out.c_str());
	FILE* fout = fopen(aux, "w");
	fprintf(fout, "file\tphi\tr\tperiod\traio_max\tpx_mean\tx_final\teta");

	printf("WRITING TO FILE %s\n", aux);

	for(int i = 0; i < n_files; i++){
		sprintf(aux, "%slogs/logOut%05d.txt", directory.c_str(), i);
		if(FILE *fin = fopen(aux, "r")){
			fgets(aux, 512, fin);
			aux[0] = 48;
			fprintf(fout, "\n%i\t%s", i, aux);
			fclose(fin);
		}else{
			fprintf(fout, "\n%i\t-1\t-1\t-1\t-1\t-1\t-1\t-1", i);
		}		
	}
	fclose(fout);

	printf("FINISHED\n");

}

int main(int argc, char **argv){
	
	/*
	CHANGE THE MAIN DIRECTORY OF THE PLOTS AND THE FILE NAME
	*/

	string directory = "../outputs/Data_Finner/"; // INCLUDE LAST SLASH
	string file_out = "Data_Out.txt";
	get_final_txt(directory, file_out, 3124);


	return 0;

}