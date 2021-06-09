#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <set>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>
#include <gemmi/unitcell.hpp>

#define PI 3.1415926535897932

using namespace std;

struct Properties {
  double area;
  double energy;
};

double molecular_mass(vector<gemmi::SmallStructure::Site> const &sites) {
  double mass = 0;
  for (auto site: sites)
    mass += site.element.weight() * site.occ;
  return mass;
}

gemmi::Vec3 randomSphereVector() {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> normal_distrib (0.0,1.0);
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}

void PrintStringVector(vector<string> v) {
	for(int i=0;i<v.size();++i)
		cout << v[i] << endl;
  // cout << endl;
}

void PrintV3(gemmi::Vec3 v){
	cout << v.x << "\t";
	cout << v.y << "\t";
	cout << v.z << endl;
}

string ReplaceString(string subject, const string& search, const string& replace) {
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
}

string strip(string sep, string inpt)
{
  int start_it = 0;
  int end_it = inpt.size()-sep.size();
  while (inpt.substr(start_it,sep.size())==sep)
    start_it = start_it + sep.size();
  while (inpt.substr(end_it,sep.size())==sep)
    end_it = end_it - sep.size();
  return inpt.substr(start_it, end_it-start_it+1);
}

vector<string> SplitString(string s, string sep){
  vector<string> v;
	string temp = "";
	for(int i=0;i<s.length();++i){
		if(s.substr(i,sep.size())==sep){
      if(temp!=""){
			  v.push_back(temp);
      }
			temp = "";
      i=i+sep.size()-1;
		}
		else{
			temp.push_back(s[i]);
		}
	}
  if(temp!=""){
    v.push_back(temp);
  }
  return v;
}

// Put it in a private / public class system
map<string, vector<string> >  ReadFF(string forcefield) {
  map<string, vector<string> > forcefield_dict;
  vector<string> L;
  string raspa_dir = getenv("RASPA_DIR");
  string path = raspa_dir + "/share/raspa/forcefield/" + forcefield + "/force_field_mixing_rules.def";
  cout << path << endl;
  ifstream MyFile(path);
  if (!MyFile) {
      cerr << "Couldn't open forcefield file.\n";
  }
  string myText;
  while (getline (MyFile, myText)) {
    cout << myText << endl;
    L.push_back(myText);
  }
  vector<string> forcefieldDefInfo(L.begin() + 7, L.end() - 2);
  PrintStringVector(forcefieldDefInfo);
  cout << L[6].substr(2) << endl;
  vector <string> columns_values = SplitString(L[6].substr(2), ", ");
  PrintStringVector(columns_values);

  for (size_t i = 0; i < forcefieldDefInfo.size(); ++i ) {
    vector<string> split_row_temp = SplitString(ReplaceString(forcefieldDefInfo[i],"\t"," "), " ");
    int k = 0;
    for (size_t j = 0; j < columns_values.size(); ++j ) {
      forcefield_dict[columns_values[j]].push_back(split_row_temp[k]);
      k++;
    }
  }
  // PrintStringVector(forcefield_dict["sigma"]);
  return forcefield_dict;
  // print into a json
  // print head(n)
}

vector<string> get_epsilon_sigma(string element, map<string, vector<string> > forcefield_dict) {
  vector<string> epsilon_sigma = {};
  int i=0; 
  // add "_" for framework atoms
  while (i<forcefield_dict["atom type"].size() && forcefield_dict["atom type"][i]!=element)
    i++;
  epsilon_sigma.push_back(strip("\t",forcefield_dict["epsilon"][i]));
  epsilon_sigma.push_back(strip("\t",forcefield_dict["sigma"][i]));
  return epsilon_sigma;
}

int main(int argc, char* argv[])
{
  auto block = gemmi::cif::read_file(argv[1]).sole_block();
  auto st = gemmi::make_small_structure_from_block(block);
  cout << st.spacegroup_hm << endl;
  cout << "Number of atoms: " << st.sites.size() << endl;

  string forcefield = "UFF";
  map<string, vector<string> > forcefield_dict = ReadFF(forcefield);

  int num_steps = 200;
  double cutoff = 12.0;
  string element_ads = "Xe";
  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);
  for (auto site: st.sites) {
    string element_host = site.type_symbol;
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;
    // cout << "element: " << element_host << endl;
    // cout << "sigma host: " << sigma_host << endl;
    // cout << "radius: " << radius << endl;
    
  }
  gemmi::Vec3 U = randomSphereVector();
  gemmi::Vec3 V = randomSphereVector();
  PrintV3(U);
  PrintV3(V);
  double d = gemmi::Position(U).dist(gemmi::Position(V));
  cout << d << endl;
  double radius = 2.0;
  U *= 2.0;
  PrintV3(U);
  cout << U.length() << endl;
}
// Understand how the PBC works
// TODO calculate the energies on surface
// LJ function