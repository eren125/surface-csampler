#include <iostream>
#include <random>
#include <set>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

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
std::default_random_engine generator;
std::normal_distribution<double> normal_distrib;
gemmi::Vec3 randomSphereVector()
{
  auto v = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
  return v.normalized();
}
int main(int argc, char* argv[])
{
  auto block = gemmi::cif::read_file(argv[1]).sole_block();
  auto st = gemmi::make_small_structure_from_block(block);
  std::cout << st.spacegroup_hm << std::endl;
  std::cout << "Number of atoms: " << st.sites.size() << std::endl;
  auto allsites = st.get_all_unit_cell_sites();
  std::cout << "Number of atoms: " << allsites.size() << std::endl;
  std::set<int> elements;
  for (auto site: st.sites)
    elements.insert(site.element.atomic_number());
  std::cout << "Atom types present:";
  for (auto i: elements)
    std::cout << " " << gemmi::Element(i).name();
  std::cout << std::endl;
  std::map<std::string, Properties> prop;
  for (auto site: allsites)
  {
    if (prop.count(site.label) != 0)
      continue;
    // Build a list of neighbors
    std::vector<std::pair<gemmi::SmallStructure::Site, gemmi::Position> > neighbors;
    for (auto n: allsites)
    {
      if (n.label == site.label)
	continue;
      auto vec_fract = (n.fract - site.fract).wrap_to_zero();
      auto vec_cart = st.cell.orthogonalize_difference(vec_fract);
      auto dist2 = vec_cart.length_sq();
      if (dist2 < 12*12)
        neighbors.push_back(std::make_pair(n, vec_cart));
    }
    Properties m;
    int surf = 0;
    for (int i = 0; i < 200; i++)
    {
      bool free = true;
      auto v = randomSphereVector();
      for (auto neigh: neighbors)
      {
	auto d = gemmi::Vec3(neigh.second) - 1.4 /*probe radius*/ * v;
	if (d.length_sq() < 1.4 * 1.4)
	{
	  free = false;
	  break;
	}
      }
      if (free)
	surf++;
    }
    m.area = double(surf) / 200;
    std::cout << site.label << " " << m.area << std::endl;
    prop[site.label] = m;
  }
  std::cout << "Molecular weigth: " << molecular_mass(allsites) << std::endl;
  std::cout << "Unit cell volume: " << st.cell.volume << std::endl;
  std::cout << "Density: " << molecular_mass(allsites) / st.cell.volume / 6.02e23 * 1e30 * 1e-6 << std::endl;
}