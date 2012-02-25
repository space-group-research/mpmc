/*  File:   Molecule.h
    Author: Brant Tudor
    Date: 6/7/2010

    Space Research Group
    University of South Florida
 
    Class definition for representation and geometic manipulation of a molecule
  
    The molecule has an x, y, and z position vector that represents the coordinates
    of the molecule in space. Additionally, there is a linked list of atoms,
    each of which has x, y, z coordinates relative to the origin. The origin
    should therefore be the center of the molecule. When asked to produce absolute
    coordinates for each atom (e.g. when producing an XYZ file), the class simply
    adds the position vector of the molecule to the relative coordinates of each atom.
 
*/

#ifndef _MOLECULE_H
#define	_MOLECULE_H

#include "Atom.h"
#include <iostream>


class Molecule {
private:
    class AtomNode
    {
    public:
        Atom      atom;
        AtomNode *next;
    
        AtomNode(int aN, double X, double Y, double Z, int ID ):next(0)
        {
            atom.set_atomic_number(aN);
            atom.set_id( ID );
            atom.set_location( X, Y, Z );
        };
    };

public:
    Molecule():x(0), y(0), z(0), atomList(0,0,0,0,0){};
    Molecule(Molecule& orig);
    virtual ~Molecule();

    Molecule &operator=( Molecule &right );

    // Set absolute position of the molecule
    void set_location( double X, double Y, double Z );
	
    void add_atom( int aN, double X, double Y, double Z );
    void add_atom( int aN, double X, double Y, double Z, int ID );
    int atom_count();

    // Get individual coordinates of the molecule
    double get_x() { return x; };
    double get_y() { return y; };
    double get_z() { return z; };

    // Get individual coordinates of an atom (atomic coordinates
    // are relative to the center of the molecule).
    double get_absolute_x( Atom atom );
    double get_absolute_y( Atom atom );
    double get_absolute_z( Atom atom );

    void displayXYZ();
    void displayMolpro();

    // Move the entire molecule
    void translate( double dX, double dY, double dZ ){ set_location( get_x()+dX, get_y()+dY, get_z()+dZ ); };

    // Rotate the molecule about the indicated axis
    void rotateX( double theta );
    void rotateY( double theta );
    void rotateZ( double theta );

private:
    double x, y, z;
    AtomNode atomList;
};

#endif	/* _MOLECULE_H */

