#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

int main(int argc, char **argv){
	
	/*
	CHANGE THE MAIN DIRECTORY OF THE PLOTS AND THE FILE NAME
	*/

	string directory = "../outputs/Data/"; // INCLUDE LAST SLASH
	string file_out = "Data_Out.txt";

	char aux[256];
	int n_files;
	sprintf(aux, "%splots/", directory.c_str());
	string* files = list_dir(aux, &n_files);

	sprintf(aux, "%s%s", directory.c_str(), file_out.c_str());
	FILE* fout = fopen(aux, "w");
	fprintf(fout, "y\tz\tphi\tr\tperiod\tthetamin\tthetamax");

	FILE* fin;
	int file_to_read;
	for(int i = 0; i < n_files; i++){
		sprintf(aux, "%splots/plotOut%%05d.txt", directory.c_str());
		sscanf(files[i].c_str(), aux, &file_to_read);
		sprintf(aux, "%slogs/logOut%05d.txt", directory.c_str(), file_to_read);
		fin = fopen(aux, "r");
		fgets(aux, 256, fin);
		fprintf(fout, "\n%s", aux);
		printf("%s\n", aux);
		fclose(fin);
	}
	fclose(fout);

}