#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <set>
#include <limits>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

using namespace std;
#define R 8.31446261815324e-3 // kJ/K/mol

string strip(string sep, string inpt) {
  int start_it = 0;
  int end_it = inpt.size()-sep.size();
  while (inpt.substr(start_it,sep.size())==sep)
    start_it = start_it + sep.size();
  while (inpt.substr(end_it,sep.size())==sep)
    end_it = end_it - sep.size();
  return inpt.substr(start_it, end_it-start_it+1);
}

string ReplaceString(string subject, const string& search, const string& replace) {
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
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

map<string, vector<string> >  ReadFF(string forcefield_path) {
  map<string, vector<string> > forcefield_dict;
  vector<string> L;
  ifstream MyFile(forcefield_path);
  if (!MyFile) {
      cerr << "Couldn't open forcefield file.\n";
  }
  string myText;
  while (getline (MyFile, myText)) {
    L.push_back(myText);
  }
  vector<string> forcefieldDefInfo(L.begin() + 7, L.end() - 2);
  vector <string> columns_values = SplitString(L[6].substr(2), ", ");

  for (size_t i = 0; i < forcefieldDefInfo.size(); ++i ) {
    vector<string> split_row_temp = SplitString(ReplaceString(forcefieldDefInfo[i],"\t"," "), " ");
    int k = 0;
    for (size_t j = 0; j < columns_values.size(); ++j ) {
      forcefield_dict[columns_values[j]].push_back(split_row_temp[k]);
      k++;
    }
  }
  return forcefield_dict;
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

gemmi::Vec3 randomSphereVector() {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> normal_distrib (0.0,1.0);
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}

bool check_free(gemmi::Vec3 ads_position, vector<tuple<double, double, gemmi::Position> > neighbors, double acc_threshold) {
  bool free = true;
  for (auto neigh: neighbors) {
    gemmi::Position pos_neigh = get<2>(neigh);
    double distance = ads_position.dist(pos_neigh);
    double sigma = get<1>(neigh);
    if (distance<sigma*acc_threshold) {
      free = false;
      break;
    }
  }
  return free;
}

double LJEnergy_shifted(gemmi::Vec3 ads_position, vector<tuple<double, double, gemmi::Position> > neighbors, double cutoff, double acc_threshold) {
  double Energy = 0;
  for(auto neigh: neighbors) {
    gemmi::Position pos_neigh = get<2>(neigh);
    double distance = ads_position.dist(pos_neigh);
    double sigma = get<1>(neigh);
    if (distance<sigma*acc_threshold) {
      Energy = std::numeric_limits<double>::infinity();
      break;
    }
    else if (distance<cutoff) {
      double epsilon = get<0>(neigh);
      double sigma = get<1>(neigh);
      Energy += 4 * epsilon * ( pow(sigma / distance,12) - pow(sigma / cutoff,12) - pow(sigma / distance,6) + pow(sigma / cutoff,6) );
    }
  }
  return R * Energy;
}

int main(int argc, char* argv[]) {
  auto t_start = std::chrono::high_resolution_clock::now();
  // Parameters
  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  int num_steps = stoi(argv[5]);
  string element_ads = argv[6];
  double acc_threshold = stod(argv[7]); 
//   double acc_threshold = 0.8;

  auto block = gemmi::cif::read_file(structure_file).sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);

  map<string, vector<string> > forcefield_dict = ReadFF(forcefield_path);

  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  // loop to create a map of potential neighbors
  int n_max = (int)floor(cutoff/structure.cell.a)+1;
  int m_max = (int)floor(cutoff/structure.cell.b)+1;
  int l_max = (int)floor(cutoff/structure.cell.c)+1;

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

  // loop over the sites to calculate LJ_energies
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  for (auto site: structure.sites) {
    string element_host = site.type_symbol;
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;

    int surf = 0;
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));

    vector<double> list_energy_lj;
    for (int i = 0; i < num_steps; i++) {
      gemmi::Vec3 V = randomSphereVector();
      V *= radius;
        surf++;
        double energy_lj = LJEnergy_shifted(V+Vsite, neighbors, cutoff, acc_threshold);
	if (energy_lj == std::numeric_limits<double>::infinity()) {
	}
	else {
	  double exp_en = exp(-energy_lj/(R*temperature));
          boltzmann_energy_lj += exp_en * energy_lj;
          sum_exp_energy += exp_en;
	}
    }
  }
  auto t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  cout << structure_file << "," << boltzmann_energy_lj/sum_exp_energy - R*temperature << "," << elapsed_time_ms/1000 << endl;
}
