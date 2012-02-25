/*  File:   Atom.h
    Author: Brant Tudor
    Date: 6/7/2010

    Space Research Group
    University of South Florida
 
    Class definition for representation of an atom
  
    Holds local (relative) coordinates and atomic number of an atom.
    Additionally, member functions provide ability to represent atomic
    number as an elemental symbol. An id field is provided in order to
    label specific atoms of interest. May add atomic weight (ability
    to handle isotopes) in the future.
 
*/

#ifndef _ATOM_H
#define	_ATOM_H

#include <string>

class Atom {
public:
    Atom();
    Atom(int);
    Atom(int aN, double x, double y, double z );
    Atom(int aN, double x, double y, double z, int ID );
    Atom(const Atom& orig);
    virtual ~Atom();

    void set_id(int i) { id = i; };
    int  get_id() { return id; };

    void set_atomic_number( int aN ) { atomicNumber = aN;   };
    int  get_atomic_number()         { return atomicNumber; };

    void   set_x( double X ) { x = X;    };
    double get_x()           { return x; };

    void   set_y( double Y ) { y = Y;    };
    double get_y()           { return y; };

    void   set_z( double Z ) { z = Z;    };
    double get_z()           { return z; };

    void set_location( double, double, double );

    std::string get_symbol();

private:
    int id;             // optional numeric identifier
    int atomicNumber;   // Atom Type
    double x, y, z;     // Spatial Location
};

#endif	/* _ATOM_H */
