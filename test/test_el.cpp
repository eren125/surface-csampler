#include <iostream>            // std::cout
#include <string>
#include <cmath>
#include <gemmi/elem.hpp>      // Atomic Element manipulation
#include <local/forcefield.hpp>      // Atomic Element manipulation

using namespace std;

// Print epsilons and sigmas into format then copy paste 4 by 4

int main() {
LennardJones::Parameters ff_params;

string host_el_str = "C";
// double epsilon = ff_params.get_epsilon(host_el_str, false);
// double sigma = ff_params.get_sigma(host_el_str, false);
auto epsilon_sigma = ff_params.get_epsilon_sigma(host_el_str, false);
double epsilon = epsilon_sigma.first;
double sigma = epsilon_sigma.second;
cout << epsilon << endl;
cout << sigma << endl;

// string forcefield_path = "force_field_mixing_rules.def";
// ff_params.read_lj_from_raspa(forcefield_path);

// cout << ff_params.host_epsilons[host_el.ordinal()] << endl;
// cout << ff_params.host_sigmas[host_el.ordinal()] << endl;

// string guest_el_str = "Xe";
// gemmi::Element guest_el(guest_el_str.c_str());
// cout << ff_params.guest_epsilons[ff_params.guest_noble_index[guest_el.ordinal()]] << endl;
// cout << ff_params.guest_sigmas[ff_params.guest_noble_index[guest_el.ordinal()]] << endl;

// typedef const char elname_t[3];
// static constexpr elname_t names[119] = {
//     "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O", "F", "Ne",
//     "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
//     "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
//     "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
//     "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
//     "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
//     "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
//     "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
//     "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
//     "Tl", "Pb", "Bi", "Po", "At", "Rn",
//     "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
//     "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
//     "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
//     "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
//   };

// cout << "EPSILON UFF" << endl;
// for (int i=0; i<119; i++) {
//   gemmi::Element host_el(names[i]);
//   cout << "/*" << names[i] << "*/ " << ff_params.host_epsilons[host_el.ordinal()] << ", ";
//   if (i % 4 == 0) {
//     cout << endl;
//   }
// }
// cout << endl;
// cout << "SIGMA UFF" << endl;
// for (int i=0; i<119; i++) {
//   gemmi::Element host_el(names[i]);
//   cout << "/*" << names[i] << "*/ " << ff_params.host_sigmas[host_el.ordinal()] << ", ";
//   if (i % 4 == 0) {
//     cout << endl;
//   }
// }
// cout << endl;
}
