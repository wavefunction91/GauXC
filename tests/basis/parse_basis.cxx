#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <fstream>
#include <iostream>
#include <map>

namespace GauXC {

BasisSet<double> parse_basis( const Molecule& mol,
                              std::string     fname,
                              SphericalType   sph    ) {


  std::ifstream infile( fname );
  if( ! infile.good() ) throw std::runtime_error(fname + " not found!");

  std::vector<std::vector<std::string>> basis_records;
  {
    std::string line;

    while( std::getline( infile, line ) ) {

      if( line.find("!") != std::string::npos ) continue;
      if( line.size() == 0 ) continue;

      // New record
      if( line.find("****") != std::string::npos ) {
        basis_records.emplace_back();
        continue;
      }

      basis_records.back().emplace_back( line );

    }
  }

  std::map<std::string,int>       atomic_number_map;
  std::map<int, BasisSet<double>> basis_shells;
  for( const auto& record : basis_records ) {
    std::string atom_line = record[0];
    std::string atom_symb = atom_line.substr(0,2);
    if( atom_symb[1] == ' ' ) atom_symb = atom_symb[0];
    std::transform( atom_symb.begin(), atom_symb.end(), atom_symb.begin(),
                    [](auto a){ return std::toupper(a); } );
    
    int Z = atomic_number_map.at(atom_symb);

    BasisSet<double> atom_basis;
    for( auto rec_it = record.begin()+1; rec_it != record.end(); ) {
      std::string type_line = *rec_it; rec_it++; // Read type line

      bool gencon = false;
      int nprim = ...;
      std::vector<double> alpha(nprim);
      std::vector<double> coeff_primary(nprim), coeff_secondary(nprim);
      for( int i = 0; i < nprim; ++i ) {
        std::string prim_line = *rec_it; rec_it++; // Read prim line
        auto prim_tokens = tokenize( prim_line, ' ');

        alpha[i]         = std::stod( prim_tokens.at(0) );
        coeff_primary[i] = std::stod( prim_tokens.at(1) );
        if( gencon )
          coeff_secondary[i] = std::stod( prim_tokens.at(2) );

      }
      
      using prim_array = Shell<double>::prim_array;
      using cart_array = Shell<double>::cart_array;

      prim_array alpha_arr, coeff_primary_arr, coeff_secondary_arr;
      std::copy( alpha.begin(), alpha.end(), alpha_arr.begin() );
      std::copy( coeff_primary.begin(), coeff_primary.end(), 
                 coeff_primary_arr.begin() );
      if( gencon )
        std::copy( coeff_secondary.begin(), coeff_secondary.end(), 
                   coeff_secondary_arr.begin() );

      SphericalType sph_use = l > 1 ? sph : SphericalType(false);
      atom_basis.emplace_back( Shell<double>(
        PrimSize(nprim), AngularMomentum(l), sph_use,
        alpha_arr, coeff_primary_arr, {0., 0., 0.}
      ));

      if( gencon )
        atom_basis.emplace_back( Shell<double>(
          PrimSize(nprim), AngularMomentum(1), SphericalType(false),
          alpha_arr, coeff_secondary_arr, {0., 0., 0.}
        ));
    }

  }

  return BasisSet<double>();
}

}
