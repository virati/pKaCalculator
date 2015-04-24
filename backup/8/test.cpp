#include <string>
#include <sstream>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
	stringstream ss;
	int x;
	
	string a;
	string b;
	
	a = " 343";
	b = "569";
	
	ss << a;
	
	//ss.str("23");
	ss >> x;
	
	cout << "test " << x << endl;
	
	//ss.str(string());
	//ss.str("");
	ss.clear();
	ss << b;
	int y = 0;
	ss >> y;
	
	cout << "tests " << y;
	
	return 0;
}

//stringstream ss(m_register.returnResInfo(x,6+y));
