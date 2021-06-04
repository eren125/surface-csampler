#include <iostream>
#include <random>
#include <set>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

#define PI 3.1415926535897932

using namespace std;
namespace cif = gemmi::cif;

int main() {
  cif::Document doc = cif::read_file("KAXQIL_clean.cif");
  for (cif::Block& block : doc.blocks)
    for (auto cc : block.find("_chem_comp.", {"id", "formula_weight"}))
      cout << cc[0] << " weights " << cc[1] << endl;
}

