#include <iostream>
#include <fstream>

#include <math.h>
#include <chrono>
#include <random>

#include <gemmi/cif.hpp>       // file -> cif::Document
#include <gemmi/smcif.hpp>     // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>  // Space Group manipulation
#include <gemmi/neighbor.hpp>  // Neighbor Search
#include <gemmi/elem.hpp>  // Neighbor Search


#define R 8.31446261815324e-3 // kJ/K/mol

using namespace std;
namespace cif = gemmi::cif;

gemmi::Vec3 randomSphereVector() {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> normal_distrib (0.0,1.0);
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}

vector<gemmi::Vec3> generateSphereRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> normal_distrib (0.0,1.0);
  for (int i = 0; i < num_steps; i++) {
    auto coord = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
    v.push_back(coord.normalized());
  }
  return v;
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
  ifstream MyFile(forcefield_path);
  if (!MyFile) {
      cerr << "Couldn't open forcefield file.\n";
  }
  string myText;
  while (getline (MyFile, myText)) {
    L.push_back(myText);
  }
  vector <string> columns_values = SplitString(L[6].substr(2), ", ");
  vector<string> forcefieldDefInfo(L.begin() + 7, L.end() - 2);

  for (size_t j = 1; j < columns_values.size(); ++j ) {
    forcefield_dict[columns_values[0]].push_back(columns_values[j]);
  }
  for (size_t i = 0; i < forcefieldDefInfo.size(); ++i ) {
    vector<string> split_row_temp = SplitString(ReplaceString(forcefieldDefInfo[i],"\t"," "), " ");
    for (size_t j = 1; j < split_row_temp.size(); ++j ) {
      forcefield_dict[split_row_temp[0]].push_back(split_row_temp[j]);
    }
  }
  return forcefield_dict;
}

vector<string> get_epsilon_sigma(string element, map<string, vector<string> > forcefield_dict) {
  vector<string> epsilon_sigma = {};
  try {
    vector<string> value = forcefield_dict.at(element);
    epsilon_sigma.push_back(value[1]);
    epsilon_sigma.push_back(value[2]);
  }
  catch (const std::out_of_range&) {
    cout << "In get_epsilon_sigma, Key \"" << element.c_str() << "\" not found" << endl;
  }
  return epsilon_sigma;
}

double LJEnergy_shifted(gemmi::Vec3 ads_position, vector<tuple<double, double, gemmi::Position> > neighbors, double cutoff_sq) {
  double Energy = 0;
  for(auto neigh: neighbors) {
    gemmi::Position pos_neigh = get<2>(neigh);
    double distance_sq = ads_position.dist_sq(pos_neigh);
    if (distance_sq<cutoff_sq) {
      double epsilon = get<0>(neigh);
      double sigma_sq = pow(get<1>(neigh),2);
      Energy += 4 * epsilon * ( pow(sigma_sq/distance_sq,6) - pow(sigma_sq/cutoff_sq,6) - pow(sigma_sq/distance_sq,3) + pow(sigma_sq/cutoff_sq,3) );
    }
  }
  return R * Energy;
}

int main(int argc, char* argv[])
{
  auto t_start = std::chrono::high_resolution_clock::now();
  // Try catch error (when file format wrong)
  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = pow(cutoff,2);
  int num_steps = stoi(argv[5]);
  string element_ads = argv[6];

  map<string, vector<string> > forcefield_dict = ReadFF(forcefield_path);

  auto doc = cif::read_file(structure_file);
  auto block = doc.sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);
  // Setup spacegroup using number and reset the images properly (don't use hm notations)
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const std::string* val = block.find_value(tag)) {
      int spacegroup_number = (int) cif::as_number(*val);
      const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(spacegroup_number);
      structure.cell.set_cell_images_from_spacegroup(sg);
      break;
    }

  // cout << "Number of unique atoms: " << structure.sites.size() << endl;
  // cout << "Number of all unitcell atoms: " << structure.get_all_unit_cell_sites().size() << endl;

  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  // // Defines the size of the supracell
  int n_max = (int)floor(1.3*cutoff/structure.cell.a)+1;
  int m_max = (int)floor(1.3*cutoff/structure.cell.b)+1;
  int l_max = (int)floor(1.3*cutoff/structure.cell.c)+1;
  // cout << n_max << " " << m_max << " " << l_max << endl;

  // Creates a list of sites within the cutoff
  vector<tuple<double, double, gemmi::Position> > supracell_sites;
  for (auto site: structure.get_all_unit_cell_sites()) {
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
          supracell_sites.push_back(make_tuple(epsilon, sigma, pos));
        }
      }
    }
  }
  // TODO Tester plusieurs techniques
  vector<gemmi::Vec3> sphere_distr_vector = generateSphereRandom(num_steps);

  double mass = 0;
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  for (auto site: structure.sites) {
    // Get LJ parameters
    string element_host = site.type_symbol;
    gemmi::Element el(element_host.c_str());
    mass += el.weight();
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;
    // Get site coordinates
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));

    vector<double> list_energy_lj;
    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      double energy_lj = LJEnergy_shifted(V+Vsite,supracell_sites, cutoff_sq);
      boltzmann_energy_lj += exp(-energy_lj/(R*temperature)) * energy_lj;
      sum_exp_energy += exp(-energy_lj/(R*temperature));
    }
  }
  double Na = 6.02214076e23;
  auto Framework_density = (structure.cell.images.size()+1)*1e-3*mass/(Na*structure.cell.volume*1e-30); // kg/m3
  auto t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Time (s)
  cout << structure_file << "," << boltzmann_energy_lj/sum_exp_energy - R*temperature << "," << 1e-3*sum_exp_energy/(R*temperature)/(structure.sites.size()*num_steps)/Framework_density << "," << elapsed_time_ms/1000 << endl;
}