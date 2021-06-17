#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <math.h>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>
#include <gemmi/unitcell.hpp>

// #define PI 3.1415926535897932
#define R 8.31446261815324e-3 // kJ/K/mol

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

string strip(string sep, string inpt) {
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
map<string, vector<string> >  ReadFF(string forcefield_path) {
  map<string, vector<string> > forcefield_dict;
  vector<string> L;
  // cout << forcefield_path << endl;
  ifstream MyFile(forcefield_path);
  if (!MyFile) {
      cerr << "Couldn't open forcefield file.\n";
  }
  string myText;
  while (getline (MyFile, myText)) {
    L.push_back(myText);
  }
  vector<string> forcefieldDefInfo(L.begin() + 7, L.end() - 2);
  // PrintStringVector(forcefieldDefInfo);
  vector <string> columns_values = SplitString(L[6].substr(2), ", ");
  // PrintStringVector(columns_values);

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

double LJEnergy_shifted(gemmi::Vec3 ads_position, vector<tuple<double, double, gemmi::Position> > neighbors, double cutoff) {
  double Energy = 0;
  for(auto neigh: neighbors) {
    gemmi::Position pos_neigh = get<2>(neigh);
    double distance = ads_position.dist(pos_neigh);
    if (distance<cutoff) {
      double epsilon = get<0>(neigh);
      double sigma = get<1>(neigh);
      Energy += 4 * epsilon * ( pow(sigma / distance,12) - pow(sigma / cutoff,12) - pow(sigma / distance,6) + pow(sigma / cutoff,6) );
    }
  }
  return R * Energy;
}

// void CreateFracToCartMatrix(double a, double b, double c, double alpha, double beta, double gamma, vector< vector<double> > &frac_to_cart) {
//   const double deg_to_rad = PI/180;
//   alpha = deg_to_rad * alpha;
//   beta = deg_to_rad * beta;
//   gamma = deg_to_rad * gamma;
//   double n = ( cos(alpha) - cos(gamma)*cos(beta) ) / sin(gamma);
//   frac_to_cart[0][0] = a;              // lx
//   frac_to_cart[0][1] = b * cos(gamma); // xy
//   frac_to_cart[0][2] = c * cos(beta);  // xz
//   frac_to_cart[1][0] = 0;              // 
//   frac_to_cart[1][1] = b * sin(gamma); // ly
//   frac_to_cart[1][2] = c * n;          // yz
//   frac_to_cart[2][0] = 0;              //
//   frac_to_cart[2][1] = 0;              //
//   frac_to_cart[2][2] = c * sqrt( pow(sin(beta),2) - pow(n,2) ); // lz
// }

int main(int argc, char* argv[])
{
  auto t_start = std::chrono::high_resolution_clock::now();
  // Try catch error (when file format wrong)
  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  int num_steps = stoi(argv[5]);
  string element_ads = argv[6];

  auto block = gemmi::cif::read_file(structure_file).sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);

  map<string, vector<string> > forcefield_dict = ReadFF(forcefield_path);

  // cout << structure.spacegroup_hm << endl;
  // cout << "Number of atoms: " << structure.sites.size() << endl;
  // cout << "a: " << structure.cell.a << endl;
  // cout << "b: " << structure.cell.b << endl;
  // cout << "c: " << structure.cell.c << endl;
  // cout << "alpha: " << structure.cell.alpha << endl;
  // cout << "beta: " << structure.cell.beta << endl;
  // cout << "gamma: " << structure.cell.gamma << endl;

  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  // // loop to create a map of potential neighbors
  int n_max = (int)round(cutoff/structure.cell.a)+1;
  int m_max = (int)round(cutoff/structure.cell.b)+1;
  int l_max = (int)round(cutoff/structure.cell.c)+1;

  // cout << n_max << " " << m_max << " " << l_max << endl;

  vector<tuple<double, double, gemmi::Position> > neighbors;
  for (auto site: structure.sites) {  
    string element_host = site.type_symbol;  
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    // Lorentz-Berthelot
    double epsilon = sqrt(stod(epsilon_sigma_temp[0])*epsilon_ads);
    double sigma = (stod(epsilon_sigma_temp[1])+sigma_ads)/2;
    for (int n = -n_max; (n<n_max+1); ++n){
      for (int m = -m_max; (m<m_max+1); ++m) {
        for (int l = -l_max; (l<l_max+1); ++l) {
          gemmi::Fractional coord = site.fract;
          coord.x = coord.x + n;
          coord.y = coord.y + m;
          coord.z = coord.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord));
          neighbors.push_back(make_tuple(epsilon, sigma, pos));
        }
      }
    }
  }// can be optimized by putting conditions on the positions (if distance from closest <cutoff)
  // cout << neighbors.size() << endl;
  // loop over the sites to calculate LJ_energies
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  for (auto site: structure.sites) {
    string element_host = site.type_symbol;
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;

    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));

    vector<double> list_energy_lj;
    for (int i = 0; i < num_steps; i++) {
      gemmi::Vec3 V = randomSphereVector();
      V *= radius;
      double energy_lj = LJEnergy_shifted(V+Vsite, neighbors, cutoff);
      boltzmann_energy_lj += exp(-energy_lj/(R*temperature)) * energy_lj;
      sum_exp_energy += exp(-energy_lj/(R*temperature));
      // cout << energy_lj << endl;
    }
  }
  auto t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  cout << structure_file << "," << boltzmann_energy_lj/sum_exp_energy - R*temperature << "," << elapsed_time_ms/1000 << endl;
}