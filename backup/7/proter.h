#ifndef ADD_H
#define ADD_H


//Define Structures/classes
class residue {
	public:
		residue(int position);

		string type;
		float alpha_c_loc[3];
		int chain_address;
		float pka_m;
		
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
		
	private:
		vector<pdb_line> docu;
		vector<pdb_line> moleculeInfo;
		vector<pdb_line> residueInfo;
		string moleculename;
};

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

//Calculation functions
float model_pka(string r_type) {
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

void menu() {
	int templine;
	int temptoken = 0;	
	bool fin = false;
	
	while(!fin) {	
		cout << "Which line would you like to see? (-2 to exit): ";
		cin >> templine;
		if(templine == -2) {
			fin = true;
			break;			
		}
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
	}
}

#endif /* ADD_H */
