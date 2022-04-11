#include <local/readff.h>      // read force_mixing_rules.def (raspa type file) and import key headers
#include <functional>          // std::function, std::negate

#include <cmath>
#include <limits>
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

#include <cstdint>

double invsqrtQuake( double number ) {
  double y = number;
  double x2 = y * 0.5;
  std::int64_t i = *(std::int64_t *) &y;
  // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
  i = 0x5fe6eb50c7b537a9 - (i >> 1);
  y = *(double *) &i;
  y = y * (1.5 - (x2 * y * y));   // 1st iteration
  return y;
  }

// Direction independent
vector<gemmi::Vec3> generateSphereNormalRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator_x (seed);
  default_random_engine generator_y (seed+500);
  default_random_engine generator_z (seed+1000);
  normal_distribution<double> normal_distrib (0.0,1.0);
  gemmi::Vec3 coord = gemmi::Vec3(0,0,0);
  double norm_sq = 0;
  for (int i = 0; i < num_steps; i++) {
    coord = gemmi::Vec3(normal_distrib(generator_x), normal_distrib(generator_y), normal_distrib(generator_z));
    norm_sq = coord.length_sq();
    v.push_back(coord*invsqrtQuake(norm_sq));
  }
  return v;
}

// Projection of a cylindre on a sphere (cos_phi being the height) Archimedes Theorem
vector<gemmi::Vec3> generateSphereAngleRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> uniform01(0.0, 1.0);
  double cos_theta; double sin_theta;
  double cos_phi; double sin_phi;
  gemmi::Vec3 coord;
  for (int i = 0; i < num_steps; i++) {
    cos_theta = cos(2 * M_PI * uniform01(generator));
    cos_phi = 1 - 2 * uniform01(generator);
    sin_theta = sqrt(1 - cos_theta * cos_theta);
    sin_phi = sqrt(1 - cos_phi * cos_phi);
    coord = gemmi::Vec3( sin_phi * cos_theta, sin_phi * sin_theta, cos_phi );
    v.push_back(coord);
  }
  return v;
}

vector<gemmi::Vec3> generateSphereCubeRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> uniform01(0.0, 1.0);
  double x; double y; double z;
  double norm_sq;
  gemmi::Vec3 coord;
  for (int i = 0; i < num_steps; i++) {
    x = uniform01(generator); y = uniform01(generator); z = uniform01(generator);
    norm_sq = x*x + y*y + z*z;
    while ( norm_sq > 1 ) {
      x = uniform01(generator); y = uniform01(generator); z = uniform01(generator);
      norm_sq = x*x + y*y + z*z;
    }
    coord = gemmi::Vec3(x,y,z);
    v.push_back(coord.normalized());
  }
  return v;
}

// Golden number spirals
vector<gemmi::Vec3> generateSphereSpirals(int num_steps) {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> uniform01(0.0, 1.0);
  vector<gemmi::Vec3> v;
  gemmi::Vec3 coord;
  double d_theta = M_PI * (3-sqrt(5)) ;
  double theta = 2 * M_PI * uniform01(generator);
  double d_cos_phi = 2.0/num_steps;
  double cos_phi = 1 + d_cos_phi/2;
  double sin_theta; double cos_theta;
  double sin_phi;
  for (int i = 0; i < num_steps; i++) {
    theta += d_theta;
    cos_theta = cos(theta);
    sin_theta = sqrt(1 - cos_theta * cos_theta);
    cos_phi -= d_cos_phi;
    sin_phi = sqrt(1 - cos_phi * cos_phi);
    coord = gemmi::Vec3( sin_phi * cos_theta, sin_phi * sin_theta, cos_phi );
    v.push_back(coord);
  }
  return v;
}

vector<gemmi::Vec3> generateSphereGeodesic(int num_steps) {
  vector<gemmi::Vec3> v;
  return v;
}

vector<gemmi::Vec3> generateSphereThomson(int num_steps) {
  vector<gemmi::Vec3> v;
  return v;
}

int main(int argc, char* argv[])
{
  // Set up Input Variables
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  int num_steps = stoi(argv[5]);
  string element_ads = argv[6];
  double access_coeff_sq = 0;
  if (argv[7]) {access_coeff_sq = stod(argv[7])*stod(argv[7]);}

  // Set key values
  double sqrt_8 = sqrt(8);
  double Na = 6.02214076e23;

  // Inialize key variables
  double distance_sq;
  double epsilon;
  double sigma_sq;
  double sigma_6;
  double a;
  double exp_energy;
  double exp_energy_energy;

  // function<double(void)> ljEnergy;
  // if (access_coeff_sq == 0) {ljEnergy = bind(ljEnergy_shifted, inv_cutoff_6, inv_cutoff_12, cref(distance_sq), cutoff_sq, cref(epsilon), cref(sigma_sq), cref(sigma_6) );}
  // else {ljEnergy  = bind(ljEnergy_shifted_acc, access_coeff_sq, inv_cutoff_6, inv_cutoff_12, cref(distance_sq), cutoff_sq, cref(epsilon), cref(sigma_sq), cref(sigma_6) ); }
  
  // Error catch
  if ( num_steps < 0 || temperature < 0 ) {throw invalid_argument( "Received negative value for the Number of Steps or the Temperature or the Accessibility Coefficient" );}
  if ( access_coeff_sq > 1 ) {throw invalid_argument( "Accessibility Coefficient above 1 (Read the purpose of this coeff)" );}

  // Read Forcefield Infos
  map<string, vector<string> > forcefield_dict = ReadFF(forcefield_path);

  // Read Structure cif files
  cif::Document doc = cif::read_file(structure_file);
  cif::Block block = doc.sole_block();
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

  vector <string> epsilon_sigma_temp = get_epsilon_sigma(element_ads, forcefield_dict);
  double epsilon_ads = stod(epsilon_sigma_temp[0]);
  double sigma_ads = stod(epsilon_sigma_temp[1]);

  // Defines the size of the supracell
  int n_max = (int)floor(1.3*cutoff/structure.cell.a)+1;
  int m_max = (int)floor(1.3*cutoff/structure.cell.b)+1;
  int l_max = (int)floor(1.3*cutoff/structure.cell.c)+1;

  // Creates a list of sites within the cutoff
  int N_sites = structure.get_all_unit_cell_sites().size();
  double supracell_sites [N_sites*(2*n_max+1)*(2*m_max+1)*(2*l_max+1)][6];
  int k = 0;
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
          supracell_sites[k][0] = pos.x;
          supracell_sites[k][1] = pos.y;
          supracell_sites[k][2] = pos.z;
          supracell_sites[k][3] = epsilon;
          supracell_sites[k][4] = sigma*sigma;
          supracell_sites[k][5] = supracell_sites[k][4]*supracell_sites[k][4]*supracell_sites[k][4];
          k++;
        }
      }
    }
  }
  // TODO Tester plusieurs techniques
  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereNormalRandom(num_steps);
  vector<gemmi::Vec3> sphere_distr_vector = generateSphereAngleRandom(num_steps);
  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereCubeRandom(num_steps);
  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereSpirals(num_steps);

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
    double radius = sqrt_8 * (sigma_ads+sigma_host)/2;
    // Get site coordinates
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));

    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      // Lennard Jones interaction energies
      double energy_lj = 0;
      gemmi::Vec3 pos_neigh;
      for(auto pos_epsilon_sigma: supracell_sites) {
        pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
        distance_sq = (V+Vsite).dist_sq(pos_neigh);
        sigma_sq = pos_epsilon_sigma[4];
        epsilon = pos_epsilon_sigma[3];
        sigma_6 = pos_epsilon_sigma[5];
        if (distance_sq < sigma_sq*access_coeff_sq) {
          energy_lj = numeric_limits<double>::infinity();
          break;
        }
        else if (distance_sq < cutoff_sq) {
          double inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
          double inv_distance_12 = inv_distance_6 * inv_distance_6;
          energy_lj += epsilon * sigma_6 * ( sigma_6 * inv_distance_12 - inv_distance_6 - sigma_6 * inv_cutoff_12 + inv_cutoff_6 );
        }
      }
      energy_lj *= 4*R;
      a = energy_lj/(R*temperature);
      exp_energy = 0;
      exp_energy_energy = 0;
      if (a < 30) {exp_energy = exp(-a); exp_energy_energy = exp_energy*energy_lj;}
      sum_exp_energy += exp_energy;
      boltzmann_energy_lj += exp_energy_energy;
    }
  }
  
  double Framework_density = (structure.cell.images.size()+1)*1e-3*mass/(Na*structure.cell.volume*1e-30); // kg/m3
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Time (s)
  cout << structure_file << "," << boltzmann_energy_lj/sum_exp_energy - R*temperature << "," << 1e-3*sum_exp_energy/(R*temperature)/(structure.sites.size()*num_steps)/Framework_density << "," << elapsed_time_ms/1000 << endl;
}