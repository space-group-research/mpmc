/*  File:   Molecule.cpp
    Author: Brant Tudor
    Date: 6/7/2010

    Space Research Group
    University of South Florida
 
    Class implementation for representation and geometic manipulation of a molecule
  
    The molecule has an x, y, and z position vector that represents the coordinates
    of the molecule in space. Additionally, there is a linked list of atoms,
    each of which has x, y, z coordinates relative to the origin. The origin
    should therefore be the center of the molecule. When asked to produce absolute
    coordinates for each atom (e.g. when producing an XYZ file), the class simply
    adds the position vector of the molecule to the relative coordinates of each atom.
 
*/


#include "Molecule.h"
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;


Molecule::Molecule(Molecule& orig):atomList(0,0,0,0,0)
{
    // Copy the position coordinates of the molecule
    x = orig.x;
    y = orig.y;
    z = orig.z;

    // Copy the list of atoms by traversing the atomList and
    // allocating a new node for each node in the original list.

    AtomNode *currentOrigNode = &(orig.atomList);
    AtomNode *currentNewNode = &atomList;
    
    while( currentOrigNode->next )
    {
        // Go to the next atom of the original
        currentOrigNode = currentOrigNode->next;
        // Save the values in the original atom
        int aN = currentOrigNode->atom.get_atomic_number();
        double X = currentOrigNode->atom.get_x();
        double Y = currentOrigNode->atom.get_y();
        double Z = currentOrigNode->atom.get_z();
        int ID   = currentOrigNode->atom.get_id();
        // Use these values to create a new atom in the copy
        currentNewNode->next = new AtomNode(aN, X, Y, Z, ID);
        // Advance the copy pointer
        currentNewNode = currentNewNode->next;
    }
}

Molecule::~Molecule()
{
    AtomNode *currentNode = atomList.next;
    AtomNode *lastNode;

    // Advance through the list of atoms and delete them
    while( currentNode )
    {
        lastNode = currentNode;
        currentNode = currentNode->next;
        delete lastNode;
    }
}

Molecule &Molecule::operator=( Molecule& right )
{
    // Avoid self-assignment
    if( &right == this ) return *this;
	
    // Copy the position coordinates of the molecule
    x = right.x;
    y = right.y;
    z = right.z;

    // Copy the list of atoms by traversing the atomList and
    // allocating a new node for each node in the original list.

    AtomNode *currentOrigNode = &(right.atomList);
    AtomNode *currentNewNode = &atomList;
    
    while( currentOrigNode->next )
    {
        // Go to the next atom of the original
        currentOrigNode = currentOrigNode->next;
        // Save the values in the original atom
        int aN = currentOrigNode->atom.get_atomic_number();
        double X = currentOrigNode->atom.get_x();
        double Y = currentOrigNode->atom.get_y();
        double Z = currentOrigNode->atom.get_z();
        int ID   = currentOrigNode->atom.get_id();
        // Use these values to create a new atom in the copy
        currentNewNode->next = new AtomNode(aN, X, Y, Z, ID);
        // Advance the copy pointer
        currentNewNode = currentNewNode->next;
    }

    return *this;
}

void Molecule::set_location( double X, double Y, double Z )
{
    x = X;
    y = Y;
    z = Z;
}


// Get absolute coordinates of an atom (atom coords are stored 
// relative to the center of the molecule
double Molecule::get_absolute_x( Atom atom )
{
    return get_x() + atom.get_x();
}
double Molecule::get_absolute_y( Atom atom )
{
    return get_y() + atom.get_y();
}
double Molecule::get_absolute_z( Atom atom )
{
    return get_z() + atom.get_z();
}



void Molecule::add_atom( int aN, double X, double Y, double Z ) { add_atom( aN, X, Y, Z, 0 ); }
void Molecule::add_atom( int aN, double X, double Y, double Z, int ID )
{
    AtomNode *currentNode = &atomList;

    // Advance currentNode to the end of the list.
    while( currentNode->next )
        currentNode = currentNode->next;

    // Allocate memory for a new Atom
    currentNode->next = new AtomNode(aN, X, Y, Z, ID);
}

void Molecule::displayXYZ()
{
    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        Atom atom = currentNode->atom;
        double X = get_absolute_x( atom );
        double Y = get_absolute_y( atom );
        double Z = get_absolute_z( atom );
        cout << setw(5) << left << atom.get_symbol() << right << setw(15) << setprecision(7) << fixed << X << setw(15) << Y << setw(15) << Z << endl;

        currentNode = currentNode->next;
    }
}


void Molecule::displayMolpro()
{
    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        Atom atom = currentNode->atom;
        double X = get_absolute_x( atom );
        double Y = get_absolute_y( atom );
        double Z = get_absolute_z( atom );
        cout << atom.get_symbol() << setw(7) << left << atom.get_id() << right << setw(15) << setprecision(7) << fixed << X << setw(15) << Y << setw(15) << Z << endl;

        currentNode = currentNode->next;
    }
}


int Molecule::atom_count()
{
    int count = 0;
    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        ++count;
        currentNode = currentNode->next;
    }

    return count;
}




// Functions for rotation about the indicated axis.
// NOTE: The individual atomic coordinates are stored as if each molecular center is the
// origin. Consequently, rotations are performed around the origin. Positions of the
// molecules are stored in the x,y,z vector outside and independent of any atomic data
// structure. These are then added to the atomic coords to yield absolute coordinates.

void Molecule::rotateX( double theta )
{
    // Calculate constants used during rotation
    const double phi = theta * 0.017453292520;  // phi = theta * PI/180
    const double SIN = sin(phi);
    const double COS = cos(phi);

    // Traverse the molecule list, applying the transformation matrix to
    // the coordinates of each atom individually

    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        Atom *atom = &(currentNode->atom);

        // Get the current location of the atom (relative to the molecular center)
        double X = atom->get_x();
        double Y = atom->get_y();
        double Z = atom->get_z();

        // Apply the rotation matrix and set the new location accordingly
        atom->set_location( X, Y*COS - Z*SIN, Y*SIN + Z*COS );

        // Advance to the next atom
        currentNode = currentNode->next;
    }
}

void Molecule::rotateY( double theta )
{
    // Calculate constants used during rotation
    const double phi = theta * 0.017453292520;  // phi = theta * PI/180
    const double SIN = sin(phi);
    const double COS = cos(phi);

    // Traverse the molecule list, applying the transformation matrix to
    // the coordinates of each atom individually

    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        Atom *atom = &(currentNode->atom);

        // Get the current location of the atom (relative to the molecular center)
        double X = atom->get_x();
        double Y = atom->get_y();
        double Z = atom->get_z();

        // Apply the rotation matrix and set the new location accordingly
        atom->set_location( X*COS + Z*SIN, Y, Z*COS - X*SIN );

        // Advance to the next atom
        currentNode = currentNode->next;
    }
}

void Molecule::rotateZ( double theta )
{
    // Calculate constants used during rotation
    const double phi = theta * 0.017453292520;  // phi = theta * PI/180
    const double SIN = sin(phi);
    const double COS = cos(phi);

    // Traverse the molecule list, applying the transformation matrix to
    // the coordinates of each atom individually

    AtomNode *currentNode = atomList.next;

    // Advance through the list of atoms
    while( currentNode )
    {
        Atom *atom = &(currentNode->atom);

        // Get the current location of the atom (relative to the molecular center)
        double X = atom->get_x();
        double Y = atom->get_y();
        double Z = atom->get_z();

        // Apply the rotation matrix and set the new location accordingly
        atom->set_location( X*COS - Y*SIN, X*SIN + Y*COS, Z );

        // Advance to the next atom
        currentNode = currentNode->next;
    }
}
