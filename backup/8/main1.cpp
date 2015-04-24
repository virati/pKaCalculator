#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

//Define Structures
struct atom {
	int type;
	int phys_loc[3];
	int chain_loc;
	int charge;
	double radius;
	
	atom *prev;
	atom *next;
};

struct residue {
	char type;
	float alpha_c_loc[3];
	int chain_address;
	
	residue *prev;
	residue *next;	
}

struct pdb_line {
	string token;
	int type;	//see int_type list for pdbFields	
}

class pdbreg {
	public:
		pdbreg();
		void printLine(int x);
		void addLine();
		
	private:
		vector<pdb_line> token;
}
	
atom mol_init;

void initialize_tree() {
			
	mol_init.type=0;
	//mol_init.phys_loc = {0,0,0};
	mol_init.chain_loc = 0;
	mol_init.charge=0;
	mol_init.radius=0;
	
	mol_init.prev = NULL;
	mol_init.next = NULL;
	
}

int pdb_tag_id(string *tag) {
		
	if(*tag == "HEADER") {
		return 1;	
	}
	else if(*tag == "HETATM") {
		return 4;	
	}
	else if(*tag == "COMPND") {
		return 5;
	}
	else if(*tag == "CONECT") {
		return 6;
	}
	else if(*tag == "SOURCE") {
		return 2;
	}
	else if(*tag == "END") {
		return 3;
	}
	return 0;
}

int pdb_grab_line() {
	
}



int main(int argc, char *argv[]) {
	//Main string vector for pdb data; each pdb line becomes new element
	vector<string> pdb_data;
	
	//Setup file input stream and temporary string for input
	ifstream pdb_origin(argv[1]);
	string temp;
	
	//Action variable to determine how to build our protein struct
	int action = 0;
	
	//Initialize our molecule
	initialize_tree();
	
	//Traverse our molecule
	atom *root = &mol_init;
	atom *writer = &mol_init;
	int traversed = 0;
	
	//Read in the pdb_data lines
	while(getline(pdb_origin, temp)) {
		pdb_data.push_back(temp);
	}
	
	//Going through the lines
	int i;
	int j;

	//Parsing the Actual lines
	string buf;
	istringstream ss;
	vector<string> tokens;
	int line_items[pdb_data.size()];
	int line_address[pdb_data.size()];
	line_address[0] = 0;
	
	
	//Loop that goes through each line, pushes each part (seperated by whitespace)
	//and then keeps track of how many components are in each line
	for(i=0; i < pdb_data.size(); i++) {
		line_items[i] = 0;
		
		ss.clear();
		ss.str(pdb_data[i]);
		
		while(ss >> buf) {
			tokens.push_back(buf);
			line_items[i]++;
		}			
		
		line_address[i+1] = line_address[i] + line_items[i];		
	}
	
	//cout << line_address[2] << " " << line_address[3] << " " << line_address[4] << " " << line_address[5] << " " << endl;
	
	
	
	
	//Output any line from the PDB
	int call_line;
	cout << "PDB data has " << pdb_data.size() << " lines" << endl;
	cout << "Which line would you like to see?: ";
	cin >> call_line;
	call_line--;
	
	for(j=line_address[call_line]; j < line_address[call_line + 1]; j++) {
			cout << tokens[j] << " ";
	}
	
	
	
	cout << endl;
	
	pdb_origin.close();
		
	return 0;
}
