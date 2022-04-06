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

double LJEnergy_shifted(gemmi::Vec3 ads_position, double epsilon_ads, double sigma_ads, map<string, vector<string> > forcefield_dict, vector<gemmi::NeighborSearch::Mark*> neighbor_marks, double cutoff) {
  double Energy = 0;
  for(auto mark: neighbor_marks) {
    gemmi::Position pos_neigh = mark->pos();
    double distance = ads_position.dist(pos_neigh);
    if (distance<cutoff) {
      string element_host = mark->element.name();
      vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
      // Lorentz-Berthelot
      double epsilon = sqrt(stod(epsilon_sigma_temp[0])*epsilon_ads);
      double sigma = (stod(epsilon_sigma_temp[1])+sigma_ads)/2;
      Energy += 4 * epsilon * ( pow(sigma / distance,12) - pow(sigma / cutoff,12) - pow(sigma / distance,6) + pow(sigma / cutoff,6) );
    }
  }
  return R * Energy;
}

// TODO use dist_sq
// TODO change type for
double LJEnergy_shifted_old(gemmi::Vec3 ads_position, vector<tuple<double, double, gemmi::Position> > neighbors, double cutoff) {
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

  // cout << "check: " << structure.cell.is_crystal() << endl;

  // cout << "Spacegroup number: " << sg->number << endl;
  // cout << "Space Group HM: " << structure.spacegroup_hm << endl;
  // cout << "Spacegroup hm: " << sg->hm << endl;

  // cout << "Number of images: " << structure.cell.images.size() << endl;
  cout << "Number of unique atoms: " << structure.sites.size() << endl;
  cout << "Number of all unitcell atoms: " << structure.get_all_unit_cell_sites().size() << endl;

  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  // // loop to create a map of potential neighbors
  int n_max = (int)floor(1.5*cutoff/structure.cell.a)+1;
  int m_max = (int)floor(1.5*cutoff/structure.cell.b)+1;
  int l_max = (int)floor(1.5*cutoff/structure.cell.c)+1;

  cout << n_max << " " << m_max << " " << l_max << endl;

  // TODO Trouver une meilleure façon de trouver les autres sites
  gemmi::SmallStructure suprastructure;
  suprastructure.cell = structure.cell;
  suprastructure.spacegroup_hm = structure.spacegroup_hm;
  suprastructure.wavelength = structure.wavelength;
  vector<tuple<double, double, gemmi::Position> > supracell_sites;
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
          gemmi::SmallStructure::Site site;
          site.fract.x = coord.x + n;
          site.fract.y = coord.y + m;
          site.fract.z = coord.z + l;
          site.type_symbol = element_host;
          const char *element_name = element_host.c_str();
          site.element = gemmi::find_element(element_name);
          suprastructure.sites.push_back(site);
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(site.fract));
          supracell_sites.push_back(make_tuple(epsilon, sigma, pos));
        }
      }
    }
  }// can be optimized by putting conditions on the positions (if distance from closest <cutoff)
  cout << suprastructure.sites.size() << endl;
  cout << suprastructure.get_all_unit_cell_sites().size() << endl;

  // TODO Don't take atoms outside
  auto NS = gemmi::NeighborSearch(structure, 1.5*cutoff);
  NS.populate(false);

  // TODO Tester plusieurs techniques
  vector<gemmi::Vec3> sphere_distr_vector = generateSphereRandom(num_steps);

  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  for (auto site: structure.sites) {
    // Get LJ parameters
    string element_host = site.type_symbol;
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;
    // Get site coordinates
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));
    // TODO Determine neighbors
    auto neighbor_marks = NS.find_atoms(gemmi::Position(Vsite), '\0', cutoff + radius);
    // cout << neighbor_marks.size() << endl;
    // for (auto mark: neighbor_marks) {
    //   cout << mark->x << " " << mark->y << " "  << mark->z << " "  << endl;
    // }

    vector<double> list_energy_lj;
    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      // TODO change how LJEnergy calculation works (use of distances)
      double energy_lj = LJEnergy_shifted(V+Vsite, epsilon_ads, sigma_ads, forcefield_dict, neighbor_marks, cutoff);
      // double energy_lj = LJEnergy_shifted_old(V+Vsite,supracell_sites, cutoff);
      boltzmann_energy_lj += exp(-energy_lj/(R*temperature)) * energy_lj;
      sum_exp_energy += exp(-energy_lj/(R*temperature));
      // TODO add Henry coefficient
    }
  }
  auto t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  cout << structure_file << "," << boltzmann_energy_lj/sum_exp_energy - R*temperature << "," << elapsed_time_ms/1000 << endl;
}
// Les images générées par symmétrie ne marche pas (trouver un autre moyen)
// get_epsilon_sigma trop lent