#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

#define PI 3.1415926535897932

using namespace std;
default_random_engine generator;
normal_distribution<double> normal_distrib;
struct Properties {
  double area;
  double energy;
};

double molecular_mass(vector<gemmi::SmallStructure::Site> const &sites)
{
  double mass = 0;
  for (auto site: sites)
    mass += site.element.weight() * site.occ;
  return mass;
}

gemmi::Vec3 randomSphereVector()
{
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}

void PrintStringVector(vector<string> v){
	for(int i=0;i<v.size();++i)
		cout << v[i] << endl;
  // cout << endl;
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
    vector<string> split_row_temp = SplitString(forcefieldDefInfo[i], " ", split_row_temp);
    int k = 0;
    for (size_t j = 0; j < columns_values.size(); ++j ) {
      forcefield_dict[columns_values[j]].push_back(split_row_temp[k]);
      if (columns_values[j]=="atom type")
        forcefield_dict["element"].push_back(strip("_",strip("\t",split_row_temp[k])));
      k++;
    }
  }
  PrintStringVector(forcefield_dict["element"]);
  return forcefield_dict;
  // print into a json
  // print head(n)
}

int main(int argc, char* argv[])
{
  auto block = gemmi::cif::read_file(argv[1]).sole_block();
  auto st = gemmi::make_small_structure_from_block(block);
  cout << st.spacegroup_hm << endl;
  cout << "Number of atoms: " << st.sites.size() << endl;
  auto allsites = st.get_all_unit_cell_sites();
  cout << "Number of atoms: " << allsites.size() << endl;
  set<int> elements;
  for (auto site: st.sites)
    elements.insert(site.element.atomic_number());
  cout << "Atom types present:";
  for (auto i: elements)
    cout << " " << gemmi::Element(i).name();
  cout << endl;
  map<string, Properties> prop;

  string forcefield = "UFF";
  map<string, vector<string> > forcefield_dict = ReadFF(forcefield);

  int num_steps = 200;
  double cutoff = 12.0;
  // string adsorbent_element = "Xe";
  // for (auto site: allsites)
  // {
  //   if (prop.count(site.label) != 0)
  //     continue;
  //   double sigma;
  //   double probe_radius = sigma * pow(2.0,1/6);
  //   // Build a list of neighbors
  //   vector<pair<gemmi::SmallStructure::Site, gemmi::Position> > neighbors;
  //   for (auto n: allsites)
  //   {
  //     if (n.label == site.label)
	//       continue;
  //     auto vec_fract = (n.fract - site.fract).wrap_to_zero();
  //     auto vec_cart = st.cell.orthogonalize_difference(vec_fract);
  //     auto dist2 = vec_cart.length_sq();
  //     if (dist2 < pow(cutoff,2) )
  //       neighbors.push_back(make_pair(n, vec_cart));
  //   }
  //   Properties m;
  //   int surf = 0;
  //   for (int i = 0; i < num_steps; i++)
  //   {
  //     bool free = true;
  //     auto v = randomSphereVector();
  //     for (auto neigh: neighbors)
  //     {
  //       auto d = gemmi::Vec3(neigh.second) - probe_radius * v;
  //       if (d.length_sq() < pow(probe_radius, 2))
  //       {
  //       free = false;
  //       break;
  //       }
  //     }
  //     if (free) {
  //       surf++;
  //     }
  //   }
  //   m.area = double(surf) / num_steps;
  //   cout << site.label << " " << m.area << endl;
  //   prop[site.label] = m;
  // }
  // cout << "Molecular weigth: " << molecular_mass(allsites) << endl;
  // cout << "Unit cell volume: " << st.cell.volume << endl;
  // cout << "Density: " << molecular_mass(allsites) / st.cell.volume / 6.02e23 * 1e30 * 1e-6 << endl;
}
// Understand how the PBC works
// TODO calculate the energies on surface
// LJ function