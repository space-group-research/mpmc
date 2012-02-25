/*  File:   Atom.cpp
    Author: Brant Tudor
    Date: 6/7/2010

    Space Research Group
    University of South Florida
 
    Class implementation for representation of an atom
  
    Holds local (relative) coordinates and atomic number of an atom.
    Additionally, member functions provide ability to represent atomic
    number as an elemental symbol. An id field is provided in order to
    label specific atoms of interest. May add atomic weight (ability
    to handle isotopes) in the future.
 
*/

#include "Atom.h"
#include <iostream>

Atom::Atom():atomicNumber(1),x(0.0),y(0.0),z(0.0),id(0) {}
Atom::Atom( int aN ):atomicNumber(aN),x(0.0),y(0.0),z(0.0),id(0) {}
Atom::Atom( int aN, double X, double Y, double Z ):atomicNumber(aN),x(X),y(Y),z(Z),id(0) {}
Atom::Atom( int aN, double X, double Y, double Z, int ID ):atomicNumber(aN),x(X),y(Y),z(Z),id(ID) {}

Atom::Atom(const Atom& orig)
{
    id = orig.id;
    atomicNumber = orig.atomicNumber;
    x = orig.x;
    y = orig.y;
    z = orig.z;
}

Atom::~Atom() {}


void Atom::set_location( double X, double Y, double Z )
{
    set_x( X );
    set_y( Y );
    set_z( Z );
}


std::string Atom::get_symbol()
{
    std::string symbols[] = { "",
      "H",                                                                                                                                                                                      "He",
      "Li", "Be",                                                                                                                                                 "B",  "C",  "N",  "O",  "F",  "Ne",
      "Na", "Mg",                                                                                                                                                 "Al", "Si", "P",  "S",  "Cl", "Ar",
      "K",  "Ca", "Sc",                                                                                     "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
      "Rb", "Sr", "Y",                                                                                      "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
      "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
      "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"
    };

    return symbols[ atomicNumber ];
}
