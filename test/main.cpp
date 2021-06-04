#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <math.h>
#include <cctype>

#define PI 3.1415926535897932

using namespace std;

// String Manipulations
string strip(const string &inpt)
{
    auto start_it = inpt.begin();
    auto end_it = inpt.rbegin();
    while (isspace(*start_it))
        ++start_it;
    while (isspace(*end_it))
        ++end_it;
    return string(start_it, end_it.base());
}

void SplitString(string s, vector<string> &v){
	string temp = "";
	for(int i=0;i<s.length();++i){
		if(s[i]==' '){
      if(temp!=""){
			  v.push_back(temp);
      }
			temp = "";
		}
		else{
			temp.push_back(s[i]);
		}
	}
  if(temp!=""){
    v.push_back(temp);
  }
}

// Vector Manipulations
void PrintStringVector(vector<string> v){
	for(int i=0;i<v.size();++i)
		cout << v[i] << endl;
}

void PrintDoubleVector(vector<double> v){
	for(int i=0;i<v.size();++i)
		cout << v[i] << endl;
}

void init_2Dvector(vector<vector<double> >& vec, int x_dim, int y_dim){
    vec.resize(x_dim);
    for(int i=0; i < vec.size(); i++)
         vec[i].resize(y_dim);
}

void Print2DVector(vector< vector<double> > myArray) {
  int width = 3, height = 3;
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
        cout << myArray[i][j] << "\t";
    }
    cout << endl;
  }
}

void MatMul3x3(vector< vector<double> > &Mat_1, vector< vector<double> > &Mat_2, vector< vector<double> > &Mat_3) {
  for(int i=0; i<3;i++) {
    for(int j=0;j<3;j++) {
      Mat_3[i][j] = 0;
      for(int k=0;k<3;k++) {
        Mat_3[i][j] += Mat_1[i][k] * Mat_2[k][j];
      }
    }
  }
}

// Extraction of values
void AssignValue(string tag, string line, double &var) {
  if (line.find(tag) != string::npos) {
    vector<string> v_temp;
    SplitString(line, v_temp);
    var = stod(v_temp[1]);
    cout << tag << " ";
    cout << var << endl;
  }
}

// Conversion matrixes
void CreateFracToCartMatrix(double a, double b, double c, double alpha, double beta, double gamma, vector< vector<double> > &frac_to_cart) {
  const double deg_to_rad = PI/180;
  alpha = deg_to_rad * alpha;
  beta = deg_to_rad * beta;
  gamma = deg_to_rad * gamma;
  double n = ( cos(alpha) - cos(gamma)*cos(beta) ) / sin(gamma);
  frac_to_cart[0][0] = a;
  frac_to_cart[0][1] = b * cos(gamma);
  frac_to_cart[0][2] = c * cos(beta);
  frac_to_cart[1][0] = 0;
  frac_to_cart[1][1] = b * sin(gamma);
  frac_to_cart[1][2] = c * n;
  frac_to_cart[2][0] = 0;
  frac_to_cart[2][1] = 0;
  frac_to_cart[2][2] = c * sqrt( pow(sin(beta),2) - pow(n,2) );
}

void CreateCartToFracMatrix(double a, double b, double c, double alpha, double beta, double gamma, vector< vector<double> > &cart_to_frac) {
  const double deg_to_rad = PI/180;
  alpha = deg_to_rad * alpha;
  beta = deg_to_rad * beta;
  gamma = deg_to_rad * gamma;
  double n = ( cos(alpha) - cos(gamma)*cos(beta) ) / sin(gamma);
  double omega = a * b * c * sqrt(1 - pow(cos(alpha),2) - pow(cos(beta),2) - pow(cos(gamma),2) + 2*cos(alpha)*cos(beta)*cos(gamma) );
  cart_to_frac[0][0] = 1/a;
  cart_to_frac[0][1] = - cos(gamma) / (a * sin(gamma));
  cart_to_frac[0][2] = b*c * (cos(alpha)*cos(gamma) - cos(beta)) / (omega * sin(gamma));
  cart_to_frac[1][0] = 0;
  cart_to_frac[1][1] = 1 / (b * sin(gamma));
  cart_to_frac[1][2] = a*c * (cos(beta)*cos(gamma) - cos(alpha)) / (omega * sin(gamma));
  cart_to_frac[2][0] = 0;
  cart_to_frac[2][1] = 0;
  cart_to_frac[2][2] = 1 / (c * sqrt( pow(sin(beta),2) - pow(n,2) ));
}

// Conversion of one
double ConvertFracToCart(double x_frac, double y_frac, double z_frac, vector<double> row) {
  int length = row.size();
  if (length==3) {
    return(x_frac*row[0] + y_frac*row[1] + z_frac*row[2]);
  } else {
    cout << "Error: not the right length" << endl;
    return(NAN);
  }
}

void CreateCartCoord(vector< vector<double> > frac_to_cart, map<string, vector<string> > &structure_dict) {
  vector<double> frac_to_cart_row_0(frac_to_cart[0].begin(), frac_to_cart[0].end());
  vector<double> frac_to_cart_row_1(frac_to_cart[1].begin(), frac_to_cart[1].end());
  vector<double> frac_to_cart_row_2(frac_to_cart[2].begin(), frac_to_cart[2].end());
  structure_dict["_atom_site_cart_x"] = {};
  structure_dict["_atom_site_cart_y"] = {};
  structure_dict["_atom_site_cart_z"] = {};
  for (size_t i=0; i<structure_dict["_atom_site_fract_x"].size(); ++i) {
    double x = ConvertFracToCart(stod(structure_dict["_atom_site_fract_x"][i]),stod(structure_dict["_atom_site_fract_y"][i]),stod(structure_dict["_atom_site_fract_z"][i]), frac_to_cart_row_0);
    double y = ConvertFracToCart(stod(structure_dict["_atom_site_fract_x"][i]),stod(structure_dict["_atom_site_fract_y"][i]),stod(structure_dict["_atom_site_fract_z"][i]), frac_to_cart_row_1);
    double z = ConvertFracToCart(stod(structure_dict["_atom_site_fract_x"][i]),stod(structure_dict["_atom_site_fract_y"][i]),stod(structure_dict["_atom_site_fract_z"][i]), frac_to_cart_row_2);
    structure_dict["_atom_site_cart_x"].push_back(to_string(x));
    structure_dict["_atom_site_cart_y"].push_back(to_string(y));
    structure_dict["_atom_site_cart_z"].push_back(to_string(z));
  }
}

// Read a cif and put main infos in a map
void ReadCif(char *path, map<string, vector<string> > &structure_dict) {
  ifstream MyFile(path);
  string myText;
  vector<string> L;
  vector<int> cell_index;
  vector<int> loop_index;
  int i =0;
  while (getline (MyFile, myText)) {
    if (myText.find("_cell") != string::npos) {
      cell_index.push_back(i);
    }
    if (myText.find("loop_") != string::npos) {
      loop_index.push_back(i);
    }
    // Save lines into a string vector
    L.push_back(myText);
    i++;
  }
  loop_index.push_back(i);
  MyFile.close();

  // cell info read in the cif file
  vector<string> cellCifInfo;
  for (size_t i = 0; i < cell_index.size(); ++i ) {
    cellCifInfo.push_back(L[cell_index[i]]);
  }

  double a, b, c;
  double alpha, beta, gamma;
  // split the lines and get values for a,b,c, alpha, beta, gamma
  for (size_t i = 0; i < cellCifInfo.size(); ++i ) {
    // cout << cellCifInfo[i] << endl;
    AssignValue("_cell_length_a", cellCifInfo[i], a);
    AssignValue("_cell_length_b", cellCifInfo[i], b);
    AssignValue("_cell_length_c", cellCifInfo[i], c);
    AssignValue("_cell_angle_alpha", cellCifInfo[i], alpha);
    AssignValue("_cell_angle_beta", cellCifInfo[i], beta);
    AssignValue("_cell_angle_gamma", cellCifInfo[i], gamma);
  }
  // Extract unitcell vector and angles to build conversion matrix (fract to cart)
  vector< vector<double> > frac_to_cart;
  init_2Dvector(frac_to_cart,3,3);
  CreateFracToCartMatrix(a, b, c, alpha, beta, gamma, frac_to_cart);
  vector< vector<double> > cart_to_frac;
  init_2Dvector(cart_to_frac,3,3);
  CreateCartToFracMatrix(a, b, c, alpha, beta, gamma, cart_to_frac);
  cout << "Fractional to cartesian matrix" << endl;
  Print2DVector(frac_to_cart);
  // Print2DVector(cart_to_frac);
  // vector< vector<double> > In;
  // init_2Dvector(In,3,3);
  // MatMul3x3(frac_to_cart, cart_to_frac, In);
  // Print2DVector(In);

  // Atom site loop
  // identify the right loop using index
  int begin_atom, end_atom;
  for (size_t k = 0; k < loop_index.size()-1; ++k) {
    if (L[loop_index[k]+1].find("_atom") != string::npos) {
      begin_atom = loop_index[k];
      end_atom = loop_index[k+1];
    }
  }
  // extract using the index
  vector<string> atomCifInfo;
  atomCifInfo = vector<string>(L.begin()+begin_atom+1,L.begin()+end_atom);
  // TODO check that _atom_site_fract_x is in the columns
  vector<string> columns_values;
  int first_indexes = 0;
  while (atomCifInfo[first_indexes].find("_atom") != string::npos) {
    string col_name = strip(atomCifInfo[first_indexes]);
    columns_values.push_back(col_name);
    structure_dict[col_name] = {};
    first_indexes++;
  }
  // loop over the Information about atoms contained in the cif
  for (size_t i = first_indexes; i < atomCifInfo.size(); ++i ) {
    vector<string> split_row_temp;
    SplitString(atomCifInfo[i], split_row_temp);
    int k = 0;
    for (size_t j = 0; j < columns_values.size(); ++j ) {
      structure_dict[columns_values[j]].push_back(split_row_temp[k]);
      k++;
    }
  }
  // Calculate cartesian coordinates
  CreateCartCoord(frac_to_cart, structure_dict);
  
  // PrintStringVector(structure_dict["_atom_site_type_symbol"]);
  // PrintStringVector(structure_dict["_atom_site_fract_x"]);
  // PrintStringVector(structure_dict["_atom_site_fract_y"]);
  // PrintStringVector(structure_dict["_atom_site_fract_z"]);
  // PrintStringVector(structure_dict["_atom_site_charge"]);
  // PrintStringVector(structure_dict["_atom_site_cart_x"]);
  // PrintStringVector(structure_dict["_atom_site_cart_y"]);
  // PrintStringVector(structure_dict["_atom_site_cart_z"]);
  // TODO if atom_site_type_symbol does not exist...
  // TODO check that the main columns used are here and ways to replace them
}

void PeriodicBoundaryCondition(vector<double> X,vector<double> Y,vector<double> Z) {

}

int main() {
  const double R = 8.31446261815324e-3;
  float T;
  cout << "Hello!\n";
  cout << "Temperature:\n";
  // cin >> T;
  T = 298.0;
  cout << R*T;
  cout << " kJ/mol\n";
  map<string, vector<string> > structure_dict;
  ReadCif((char *)"KAXQIL_clean.cif", structure_dict);
  // PrintStringVector(structure_dict["_atom_site_cart_x"]);
  // vector<string> Vect;
  // string str = " a  b    c       d";
  // SplitString(, Vect);
  // PrintVector(Vect);
  // string str_strip = strip(str);
  // cout << str << endl;
  return 0;
}