#include <local/forcefield.hpp> // define forcefield parameters
#include <local/sphere.hpp>     // Define sampling methodology

#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/smcif.hpp>      // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>   // Space Group manipulation
#include <gemmi/unitcell.hpp>
#include <unordered_map>

#define R 8.31446261815324e-3 // kJ/mol/K
#define sqrt_2 1.414213562373095
#define min_factor 1.122462048309373  // 2^(1/6)
#define N_A 6.02214076e23    // part/mol

using namespace std;
namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
  // Set up Input Variables
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  auto structure_file = argv[1];
  
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
  int count=0;
  for ( gemmi::SmallStructure::Site site: unique_sites ) {
    count++;
  }
  cout << count << endl;
}
