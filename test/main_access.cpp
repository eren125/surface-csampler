#include <local/forcefield.hpp> // define forcefield parameters
#include <local/sphere.hpp>     // Define sampling methodology

#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/smcif.hpp>      // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>   // Space Group manipulation
#include <gemmi/unitcell.hpp>

// Set key constants
#define R 8.31446261815324e-3 // kJ/mol/K
#define sqrt_2 1.414213562373095
#define N_A 6.02214076e23    // part/mol

using namespace std;
namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
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
  string element_guest_str = argv[6];
  double access_coeff_sq = 0;
  if (argv[7]) {access_coeff_sq = stod(argv[7])*stod(argv[7]);}
  // Error catch
  if ( num_steps < 0 || temperature < 0 ) {throw invalid_argument( "Received negative value for the Number of Steps or the Temperature or the Accessibility Coefficient" );}
  if ( access_coeff_sq > 1 ) {throw invalid_argument( "Accessibility Coefficient above 1 (Read the purpose of this coeff)" );}

  // Read Forcefield Infos
  LennardJones::Parameters ff_params;
  if (forcefield_path != "DEFAULT") {
    ff_params.read_lj_from_raspa(forcefield_path);
  }

  // Inialize key variables
  string element_host_str;
  double dist = 0;
  double distance_sq = 0;
  double epsilon = 0;
  double sigma = 0;
  double sigma_sq = 0;
  double sigma_6 = 0;
  double exp_energy = 0;

  // Adsorption infos from forcefield dictionary
  pair<double,double> epsilon_sigma = ff_params.get_epsilon_sigma(element_guest_str, true);
  double epsilon_guest = epsilon_sigma.first;
  double sigma_guest = epsilon_sigma.second;

  // Read Structure cif files
  cif::Document doc = cif::read_file(structure_file);
  cif::Block block = doc.sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);
  // Setup spacegroup using number and reset the images properly (don't use hm notations)
  int spacegroup_number = 1;
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const string* val = block.find_value(tag)) {
      spacegroup_number = (int) cif::as_number(*val);
      break;
    }
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(spacegroup_number);
  structure.cell.set_cell_images_from_spacegroup(sg);
  vector<gemmi::SmallStructure::Site> unique_sites = structure.sites;
  vector<gemmi::SmallStructure::Site> all_sites = structure.get_all_unit_cell_sites();

  // Cell parameters
  double a = structure.cell.a; double b = structure.cell.b; double c = structure.cell.c; 
  double deg_rad = M_PI/180;
  double alpha = deg_rad*structure.cell.alpha; double beta = deg_rad*structure.cell.beta; double gamma = deg_rad*structure.cell.gamma;
  // Definition of the enlarged cutoff
  double cell_diag_halved = 0.5 * sqrt( a*a + b*b + c*c + 2*(b*c*cos(alpha) + c*a*cos(beta) + a*b*cos(gamma)) );
  double large_cutoff = cutoff + sigma_guest + cell_diag_halved;
  // Cell vector definition
  double a_x = structure.cell.orth.mat[0][0]; double b_x = structure.cell.orth.mat[0][1]; double c_x = structure.cell.orth.mat[0][2];
  double a_y = structure.cell.orth.mat[1][0]; double b_y = structure.cell.orth.mat[1][1]; double c_y = structure.cell.orth.mat[1][2];
  double a_z = structure.cell.orth.mat[2][0]; double b_z = structure.cell.orth.mat[2][1]; double c_z = structure.cell.orth.mat[2][2]; 
  // Minimal box setting for a triclinic cell to contain a sphere of radius `large_cutoff`
  int n_max = abs(int(large_cutoff * sqrt((b_y*c_x-b_x*c_y)*(b_y*c_x-b_x*c_y) 
                      + (b_z*c_x-b_x*c_z)*(b_z*c_x-b_x*c_z) 
                      + (b_z*c_y-b_y*c_z)*(b_z*c_y-b_y*c_z)) 
                      / (a_z*b_y*c_x - a_y*b_z*c_x - a_z*b_x*c_y 
                      + a_x*b_z*c_y + a_y*b_x*c_z - a_x*b_y*c_z))) + 1;
  int m_max = abs(int(large_cutoff * sqrt((c_y*a_x-c_x*a_y)*(c_y*a_x-c_x*a_y) 
                      + (c_z*a_x-c_x*a_z)*(c_z*a_x-c_x*a_z) 
                      + (c_z*a_y-c_y*a_z)*(c_z*a_y-c_y*a_z)) 
                      / (b_z*c_y*a_x - b_y*c_z*a_x - b_z*c_x*a_y 
                      + b_x*c_z*a_y + b_y*c_x*a_z - b_x*c_y*a_z))) + 1;
  int l_max = abs(int(large_cutoff * sqrt((a_y*b_x-a_x*b_y)*(a_y*b_x-a_x*b_y) 
                      + (a_z*b_x-a_x*b_z)*(a_z*b_x-a_x*b_z) 
                      + (a_z*b_y-a_y*b_z)*(a_z*b_y-a_y*b_z)) 
                      / (c_z*a_y*b_x - c_y*a_z*b_x - c_z*a_x*b_y 
                      + c_x*a_z*b_y + c_y*a_x*b_z - c_x*a_y*b_z))) + 1;

  // Creates a list of sites within the cutoff
  vector<array<double,6>> supracell_sites;
  gemmi::Fractional coord_temp;
  string element_host_str_temp = "X";

  for (auto site: all_sites) {
    element_host_str = site.type_symbol;
    if (element_host_str != element_host_str_temp) {
      epsilon_sigma = ff_params.get_epsilon_sigma(element_host_str, false);
      // Lorentz-Berthelot
      epsilon = sqrt( epsilon_sigma.first * epsilon_guest );
      sigma = 0.5 * ( epsilon_sigma.second+sigma_guest );
    }
    element_host_str_temp = element_host_str;
    gemmi::Fractional coord = site.fract;
    for (int n = -n_max; (n<n_max+1); ++n){
      for (int m = -m_max; (m<m_max+1); ++m) {
        for (int l = -l_max; (l<l_max+1); ++l) {
          // calculate a distance from centre box
          array<double,6> pos_epsilon_sigma;
          coord_temp.x = coord.x + n;
          coord_temp.y = coord.y + m;
          coord_temp.z = coord.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord_temp));
          pos_epsilon_sigma[0] = pos.x;
          pos_epsilon_sigma[1] = pos.y;
          pos_epsilon_sigma[2] = pos.z;
          pos_epsilon_sigma[3] = epsilon;
          pos_epsilon_sigma[4] = sigma * sigma;
          pos_epsilon_sigma[5] = pos_epsilon_sigma[4] * pos_epsilon_sigma[4] * pos_epsilon_sigma[4];
          supracell_sites.push_back(pos_epsilon_sigma);
        }
      }
    }
  }

  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereNormalRandom(num_steps);
  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereAngleRandom(num_steps);
  // vector<gemmi::Vec3> sphere_distr_vector = generateSphereCubeRandom(num_steps);
  vector<gemmi::Vec3> sphere_distr_vector = generateSphereSpirals(num_steps);

  double mass = 0;
  double boltzmann_energy_lj = 0;

  double sum_exp_energy = 0;
  double inv_distance_6;
  double inv_distance_12;
  double sigma_host;
  double radius;
  double energy_lj;
  vector<array<double,6>> neighbor_sites;
  int count_acc = 0;
  double surface = 0;

  for ( gemmi::SmallStructure::Site site: unique_sites ) {
    // Get LJ parameters
    element_host_str = site.type_symbol;
    sigma_host = ff_params.get_sigma(element_host_str, false);
    radius = sqrt_2 * (sigma_guest+sigma_host);
    surface += 4 * M_PI * radius * radius;

    gemmi::Element el(element_host_str.c_str());
    mass += el.weight();
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));
    // Cell list pruning to have only the sites that are within (cutoff + radius) of the unique site
    neighbor_sites = {};
    for (array<double,6> pos_epsilon_sigma : supracell_sites) {
      dist = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]).dist(Vsite);
      if (dist < cutoff + radius) {neighbor_sites.push_back(pos_epsilon_sigma);}
    }
    // Loop around the sphere surface of the unique site
    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      // Lennard Jones interaction energies
      energy_lj = 0;
      gemmi::Vec3 pos_neigh;
      for(array<double,6> pos_epsilon_sigma : neighbor_sites) {
        pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
        distance_sq = (V+Vsite).dist_sq(pos_neigh);
        sigma_sq = pos_epsilon_sigma[4];
        if (distance_sq < sigma_sq*access_coeff_sq) {
          // epsilon = pos_epsilon_sigma[3];
          // sigma_6 = pos_epsilon_sigma[5];
          // inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
          // inv_distance_12 = inv_distance_6 * inv_distance_6;
          // energy_lj = epsilon * sigma_6 * ( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
          energy_lj = 1e10;
          // count_acc -= 1;
          break;
        }
        else if (distance_sq < cutoff_sq) {
          epsilon = pos_epsilon_sigma[3];
          sigma_6 = pos_epsilon_sigma[5];
          inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
          inv_distance_12 = inv_distance_6 * inv_distance_6;
          energy_lj += epsilon * sigma_6 * ( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
        }
      }
      energy_lj *= 4*R;
      if ( energy_lj <= -R*temperature ){
        count_acc++;
      }
      exp_energy = exp(-energy_lj/(R*temperature)); 
      sum_exp_energy += exp_energy;
      boltzmann_energy_lj += exp_energy*energy_lj;
      // count_acc++;
    }
  }
  double sym_images = structure.cell.images.size()+1;
  double Framework_density = sym_images * 1e-3 * mass/(N_A*structure.cell.volume*1e-30); // kg/m3
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry_surface = 1e-3*sum_exp_energy/(R*temperature)/(unique_sites.size()*num_steps)/Framework_density;    // mol/kg/Pa
  double area_accessible = 1e7 * sym_images * surface * count_acc / (num_steps * unique_sites.size() * structure.cell.volume * Framework_density); // m2/g

  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  cout << structure_file << "," << enthalpy_surface << "," << henry_surface << "," << area_accessible << "," << elapsed_time_ms*0.001 << endl;
}
