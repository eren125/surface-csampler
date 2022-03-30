#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <set>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

using namespace std;
#define R 8.31446261815324e-3 // kJ/K/mol

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


double molecular_mass(std::vector<gemmi::SmallStructure::Site> const &sites)
{
  double mass = 0;
  for (auto site: sites)
    mass += site.element.weight() * site.occ;
  return mass;
}
struct Properties {
  double area;
  double energy;
};

gemmi::Vec3 randomSphereVector() {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> normal_distrib (0.0,1.0);
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}

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

int main(int argc, char* argv[])
{
  auto t_start = std::chrono::high_resolution_clock::now();
  string structure_file = argv[1];
  int num_steps = 200;
  auto block = gemmi::cif::read_file(structure_file).sole_block();
  auto st = gemmi::make_small_structure_from_block(block);
  // std::cout << st.spacegroup_hm << std::endl;
  // std::cout << "Number of atoms: " << st.sites.size() << std::endl;
  auto allsites = st.get_all_unit_cell_sites();
  // std::cout << "Number of atoms: " << allsites.size() << std::endl;
  // std::set<int> elements;
  // for (auto site: st.sites)
  //   elements.insert(site.element.atomic_number());
  // std::cout << "Atom types present:";
  // for (auto i: elements)
  //   std::cout << " " << gemmi::Element(i).name();
  // std::cout << std::endl;
  
  std::map<std::string, Properties> prop;

  string element_ads = "Xe";
  double T = 298.0;
  string forcefield_path = "force_field_mixing_rules.def";
  map<string, vector<string> > forcefield_dict = ReadFF(forcefield_path);
  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  double boltzmann_energy_proxy = 0;
  double sum_exp_energy = 0;
  double total_area = 0;
  for (auto site: allsites)
  {
    if (prop.count(site.label) != 0)
      continue;
    string element_host = site.type_symbol;
    vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_host + "_", forcefield_dict);
    double epsilon_host = stod(epsilon_sigma_temp[0]);
    double sigma_host = stod(epsilon_sigma_temp[1]);
    double radius = pow(2,1/6) * (sigma_ads+sigma_host)/2;
        // Build a list of neighbors
    std::vector<std::pair<gemmi::SmallStructure::Site, gemmi::Position> > neighbors;
    for (auto n: allsites)
    {
      if (n.label == site.label)
	      continue;
      auto vec_fract = (n.fract - site.fract).wrap_to_zero();
      auto vec_cart = st.cell.orthogonalize_difference(vec_fract);
      auto dist2 = vec_cart.length_sq();
      if (dist2 < 4*radius*radius)
        neighbors.push_back(std::make_pair(n, vec_cart));
    }
    // std::cout << neighbors.size() << std::endl;
    Properties m;
    int surf = 0;
    for (int i = 0; i < num_steps; i++)
    {
      bool free = true;
      auto v = randomSphereVector();
      v *= radius;
      for (auto neigh: neighbors)
      {
	      auto d = gemmi::Vec3(neigh.second) - v;
	      if (d.length_sq() < radius * radius)
	      {
	        free = false;
	        break;
	      }
      }
      if (free)
	      surf++;
    }
    m.area = double(surf) / num_steps;
    // std::cout << site.label << " " << m.area << " " << -4 * epsilon_host << std::endl;
    prop[site.label] = m;
    if (m.area != 0) {
      boltzmann_energy_proxy += -4 * m.area * R * epsilon_host * exp(4*epsilon_host/T);
      sum_exp_energy += m.area * exp(4*epsilon_host/T);
    }
  }
  // std::cout << "Molecular weigth: " << molecular_mass(allsites) << std::endl;
  // std::cout << "Unit cell volume: " << st.cell.volume << std::endl;
  // std::cout << "Density: " << molecular_mass(allsites) / st.cell.volume / 6.02e23 * 1e30 * 1e-6 << std::endl;
  auto t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  std::cout << structure_file << ',' << boltzmann_energy_proxy/sum_exp_energy << ',' << elapsed_time_ms/1000 << std::endl;
}