#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include "prot_gl.h"

//#include "proter.h"
using namespace std;

//ProPka Variables
float Chb = 0;
float d1 = 0;
float d2 = 0;
float Clocal=0.001;
float Cglobal=0.001;
float Cchgchg=0;

int tempinput;

//Define Structures/classes
class atom {
	public:
	string type;
	string name;
	float locat[3];
	int pdb_address;
};

class residue {
	public:
		residue(int position);

		string type;
		float alpha_c_loc[3];
		int chain_address;
		float pka_m;
		float del_pka;
		float tpka;
		string important_atom;
		
		residue *prev;
		residue *next;	
};

residue::residue(int position) {
	type='0';
	chain_address=position;
	prev=NULL;	
}

struct pdb_line {
	string line;
	vector<string> tokens;
	int type;	//see int_type list for pdbFields	
};

class pdbreg {
	public:
		pdbreg();
		pdbreg(string name);
		string pLine(int x);
		string pLine(int x, int y);
		string pAtom(int x);
		void addpdbLine(pdb_line *new_item);
		string returnName();
		void splitLines();
		void printFull(int x);
		void generateMolinfo();
		void generateResinfo();
		int getsize(int x);
		string returnResInfo(int x,int y);
		int fillAtoms();
		atom* atomreg;
		
		
	private:
		vector<pdb_line> docu;
		vector<pdb_line> moleculeInfo;
		vector<pdb_line> residueInfo;
		string moleculename;
		
		
};

int pdbreg::fillAtoms() {
	
	atomreg = new atom[moleculeInfo.size()];
	stringstream ss;
	
	float tempLoc;
	int tempInt;
	int x=0;
	
	for(x=0; x<moleculeInfo.size(); x++) {
		atomreg[x].type = moleculeInfo[x].tokens[11];
		atomreg[x].name = moleculeInfo[x].tokens[3];
		for(int y=0; y<3; y++) {
			ss.clear();
			ss.str("");
			
			ss << moleculeInfo[x].tokens[6+y];
			ss >> tempLoc;
			atomreg[x].locat[y] = tempLoc;
		}
		
		ss.clear();
		ss.str("");
		ss << moleculeInfo[x].tokens[1];
		ss >> tempInt;
		atomreg[x].pdb_address = tempInt;
		
	}
	return x;
}

//Calculation functions

string pdbreg::returnResInfo(int x, int y) {
	return residueInfo[x].tokens[y];
}

int pdbreg::getsize(int x) {
	if(x == 2) {
		return moleculeInfo.size();
	}
	else if(x == 3) {
		return residueInfo.size();
	}
}

void pdbreg::generateMolinfo() {
	for(int x=0; x < docu.size(); x++) {
		if(docu[x].tokens[0] == "ATOM" | docu[x].tokens[0] == "HETATM") {
			moleculeInfo.push_back(docu[x]);
		}
	}
}

void pdbreg::generateResinfo() {
	int currentNum = 0;
	int tierNum = 0;
	for(int x=0; x < docu.size(); x++) {
		if(docu[x].tokens[0] == "ATOM" | docu[x].tokens[0] == "HETATM") {
			stringstream ss(docu[x].tokens[5]);
			ss >> currentNum;
			
			if(!(currentNum == tierNum) && docu[x].tokens[2] == "CA") {
				residueInfo.push_back(docu[x]);
				tierNum = currentNum;
			} 				
		}
	}
}


string pdbreg::pAtom(int x) {
	return moleculeInfo[x].line;
}

void pdbreg::printFull(int x) {
	if(x==1) {
		for(int x=0; x < docu.size();x++) {
			cout << docu[x].line << endl;
		}
	}
	else if(x==2) {
		for(int x=0; x < moleculeInfo.size();x++) {
			cout << moleculeInfo[x].line << endl;
		}
	}
	else if(x==3) {
		for(int x=0; x < residueInfo.size();x++) {
			cout << residueInfo[x].line << endl;
		}
	}
	
}

void pdbreg::splitLines() {
	string buf;
	istringstream ss;
	
	for(int x=0;x < docu.size();x++) {
		ss.clear();
		ss.str(docu[x].line);
		
		while(ss >> buf) {
			docu[x].tokens.push_back(buf);
		}			
	}
}

void pdbreg::addpdbLine(pdb_line *new_item) {
	docu.push_back(*new_item);
}

string pdbreg::pLine(int x) {
	return docu[x].line;
}

string pdbreg::pLine(int x, int y) {
	
	return docu[x].tokens[y];
}

pdbreg::pdbreg(string name) {
	moleculename = name;
}

pdbreg::pdbreg() {
	moleculename="transient";
}

string pdbreg::returnName() {
	return moleculename;
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

//Begin Program stuff
residue seed(0);
pdbreg m_register;

residue *moveDown(residue *root, int move) {
	residue *parser;
	parser=root->next;
	for(int x=0; x < move; x++) {
		parser = parser->next;
	}
	
	return parser;
}

//Menu Code
void menu(residue *root, int length) {
	int templine;
	int temptoken = 0;	
	bool fin = false;
	int option;
	int inpres=0;
	residue *rparser = root;

	
	while(!fin) {	
		cout << endl << "[1] - To See lines" << endl << "[2] - To See Residues" << endl << "Option: ";
		option = 0;
		cin >> option;
		rparser = root;	
		
		switch(option) {
			case 0:
				fin = true;
				break;
				
			case 1:
				cout << "Which line would you like to see?: ";
				cin >> templine;
				cout << endl << "Which element? (0 to see whole line) ";
				cin >> temptoken;
				
				templine--;
				if(temptoken == 0) {
					cout << endl << "This is the first line: " << endl << m_register.pLine(templine) << endl;
				}
				else if(temptoken == -1) {
					m_register.printFull(1);
					cout << endl;
				}
				else {
					cout << endl << "---------------------------------";
					cout << endl << "This is the first line: " << endl << m_register.pLine(templine) << endl << endl;
					cout << "This is the " << temptoken << " token: " << endl << m_register.pLine(templine,temptoken-1) << endl << "-----------------------" << endl << endl;
				}	
				break;
				
			case 2:
				
				cout << "which residue would you like to see?: ";
				cin >> inpres;
				//cout << "0";
				inpres--;
				
				//cout << "1";
				
				/*for(int z=0; z < inpres+1; z++) {
					rparser = rparser->next;
				}
				*/
				rparser = moveDown(rparser,inpres+1);
				
				cout << "2";
				
				cout << endl << endl << "Residue is:" << endl;
				cout << rparser->type << " at " << rparser->chain_address << endl << "Coords: ";
				for(int z=0; z<3; z++) {
					cout << rparser->alpha_c_loc[z] << " ";
				}
				cout << endl;
			break;
				
			case 3:
				rparser = rparser->next;				
				for(int x=0; x < length; x++) {
					cout << rparser->type << ": pka of " << rparser->pka_m << endl;
					rparser = rparser->next;
				}
			break;
			
		}
	}
}

float return_m_pka(string r_type) {
	if(r_type == "ARG") {
		return 12.5;
	}
	else if(r_type == "ASP") {
		return 3.9;
	}
	else if(r_type == "CYS") {
		return 8.3;
	}
	else if(r_type == "GLU") {
		return 4.3;
	}
	else if(r_type == "HIS") {
		return 6.0;
	}
	else if(r_type == "LYS") {
		return 10.5;
	}
	else if(r_type == "SER") {
		return 13;
	}
	else if(r_type == "THR") {
		return 13;
	}
	else if(r_type == "TYR") {
		return 10.1;
	}	
	else {
		return 0;
	}
}

string return_main_residue(string r_type) {
	if(r_type == "ARG") {
		return "NH1";
	}
	else if(r_type == "ASP") {
		return "OE1";
	}
	else if(r_type == "CYS") {
		return "S1";
	}
	else if(r_type == "GLU") {
		return "OE1";
	}
	else if(r_type == "HIS") {
		return "NE2";
	}
	else if(r_type == "LYS") {
		return "NZ";
	}
	else if(r_type == "SER") {
		return "OE1";
	}
	else if(r_type == "THR") {
		return "OE1";
	}
	else if(r_type == "TYR") {
		return "OE1";
	}	
	else {
		return 0;
	}
}

float return_distance(float td1[3], float td2[3]) {
	float dist;

	dist = sqrt(pow(td1[1]-td2[1],2) + pow(td1[2]-td2[2],2) + pow(td1[3] - td2[3],2));
	return dist;
}

float return_angle(float td1[3], float td2[3]) {
	float angle;
	
	float magt1 = sqrt(pow(td1[1],2) + pow(td1[2],2) + pow(td1[3],2));
	float magt2 = sqrt(pow(td2[1],2) + pow(td2[2],2) + pow(td2[3],2));
	
	angle = (td1[1]*td2[1] + td1[2]*td2[2] + td1[3]*td2[3])/(magt1 * magt2);
}

void apply_m_pka(residue *root, int length) {
	residue *parser;
	
	parser = root->next;
	
	for(int x=0; x < length; x++) {
		parser->pka_m = return_m_pka(parser->type);
		//parser->important_atom = return_main_residue(parser->type);
		parser = parser->next;
	}

}

float return_globaldes_pka(residue *target, int atomreg_size) {	
	int N15a = 0;
	float rposition[3];
	
	float del_dist;
	float atompos[3];
	
	//residue *target = rparser;
	//int temp1;
	
	for(int y=0; y<2; y++) {
		rposition[y] = target->alpha_c_loc[y];
	}

	for(int x=0; x<atomreg_size; x++) {
		for(int y=0; y<2; y++) {
			//cout << endl << m_register.atomreg[x].locat[0];
			atompos[y]=m_register.atomreg[x].locat[y];
		}
		del_dist=return_distance(rposition,atompos);
		if(del_dist < 15.5 && m_register.atomreg[x].type != "H") {
			N15a++;
		}
	}
	cout << endl << "The number of Atoms:" << N15a << endl;
	float chngPka = (N15a - 400) * Cglobal;
	return chngPka;
}

float return_localdes_pka(residue *target, int atomreg_size, float Rlocal) {	
	int Na = 0;
	float rposition[3];
	
	float del_dist;
	float atompos[3];
	
	//residue *target = rparser;
	//int temp1;
	
	for(int y=0; y<2; y++) {
		rposition[y] = target->alpha_c_loc[y];
	}

	for(int x=0; x<atomreg_size; x++) {
		for(int y=0; y<2; y++) {
	//		cout << endl << m_register.atomreg[x].locat[0];
			atompos[y]=m_register.atomreg[x].locat[y];
		}
		del_dist=return_distance(rposition,atompos);
		if(del_dist < Rlocal && m_register.atomreg[x].type != "H") {
			Na++;
		}
	}
	cout << endl << "The number of Atoms:" << Na << endl;
	float chngPka = (Na) * Clocal;
	return chngPka;
}

float return_hb_sdc(residue *target, int atomreg_size) {
	int central_position = target->chain_address;
	
	
}


//Main Function
int main(int argc, char *argv[]) {
	//Setup file input stream and temporary string for input
	ifstream pdb_origin(argv[1]);
	string temp;
	
	
	pdb_line *next_pdb;
	
	while(getline(pdb_origin, temp)) {
		next_pdb = new pdb_line();
		next_pdb->line = temp;
		
		m_register.addpdbLine(next_pdb);
	}
	
	m_register.splitLines();
	
	//Choices here:
	//Genmolinfo generates just the atomic/hetatm information
	//genresinfo generate the residues of the molecule
	
	
	cout << "here is the residue info" << endl << endl;
	m_register.generateMolinfo();
	m_register.generateResinfo();
	
	//Print the residue info (x=3)
	m_register.printFull(3);
	
	//Construct our molecule
	residue *rparser;
	residue *root;
	residue *tail;
	
	root = &seed;
	rparser = root;
	
	float temp_pos;
	int mol_size = 0;
	
	stringstream ss("");
	
	for(int x=0; x < m_register.getsize(3);x++) {
		rparser->next = new residue(x+1);
		rparser->next->type = m_register.returnResInfo(x,3);
		
		for(int y=0; y < 3; y++) {
			ss.clear();
			ss.str("");
			ss << m_register.returnResInfo(x,6+y);
			ss >> temp_pos;
			//cout << temp_pos << endl;
			//cout << endl;
			////cout << m_register.returnResInfo(x,6+y);
			
			rparser->next->alpha_c_loc[y] = temp_pos;
		}
		ss.str("");
		ss.clear();
		ss << m_register.returnResInfo(x,5);
				
		ss >> rparser->next->chain_address;
		rparser->next->prev = rparser;
		rparser = rparser->next;
		mol_size++;
	}
	
	
	tail = rparser;
	rparser = root;
		
//	menu(root,mol_size);

	cout << endl;
	

	apply_m_pka(root,mol_size);
	//cout << "hello!";
	//menu(root, mol_size);	
	
		
	int atomreg_size = m_register.fillAtoms();
	//cout << atomreg_size << endl;
	
	rparser = root;
	cout << "test";
	float chngPka[mol_size];// = 0;
	
	for(int x=0; x<mol_size; x++) {
		chngPka[x] = 0.0;
		if(rparser->pka_m != 0.0) {
			chngPka[x] += return_globaldes_pka(rparser,atomreg_size);
			chngPka[x] += return_localdes_pka(rparser,atomreg_size,1);
		}
	//cout << chngPka << endl;
		rparser->tpka = rparser->pka_m + chngPka[x];
		cout << endl << rparser->tpka;
		rparser = rparser->next;
	}
	
	rparser = root->next;
	cout << "This is the added pka: " << endl;
	for(int x=0; x<mol_size; x++) {
		cout << rparser->type << rparser->chain_address << " is now " << rparser->tpka << endl;
		rparser = rparser->next;
	}
	
	
	
	/*
	 *
	//OpenGL Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(500,500);
	glutInitWindowPosition(100,100);
	glutCreateWindow("PDB Viewer");
	glutDisplayFunc(display);
	//glutIdleFunc(display);
	glutMainLoop();
		
	//Main run loop	
	/*
	 * m_register.printFull(1);
	 * 
	cout << endl;

	menu();	
	*/
}


