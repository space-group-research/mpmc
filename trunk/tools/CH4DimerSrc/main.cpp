/*  File: main.cpp
    Author: Brant Tudor
    Date: 6/7/2010

    Space Research Group
    University of South Florida
 
    This program mathematically constructs a tetrahedral methane molecule,
    duplicates it, then transforms the molecules according to inputs passed on the
    command line. All the data necessary to construct an XYZ file is then output
    to the screen, suitable for redirection (piping) to an XYZ file.
 
*/



#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include "Molecule.h"

using namespace std;

void printInstructions(char*);



int main(int argc, char** argv) {

    double r = 0.0;         // separation distance between molecules
    double bondLength = 0.0;
    double angleZ1 = 0.0;   // angle of rotation about Z-axis rotation to be applied to 1st molecule
    double angleX1 = 0.0;   // "                     " X-axis "                                    "
    double angleY1 = 0.0;   // "                     " Y-axis "                                    "
    double angleZ2 = 0.0;   // "                     " Z-axis rotation to be applied to 2nd molecule
    double angleX2 = 0.0;   // "                     " X-axis "                                    "
    double angleY2 = 0.0;   // "                     " Y-axis "                                    "

    typedef enum { xyz, molpro } Mode;
    Mode mode = xyz;
	
    // Check that we have recieved all eight arguments
    if( argc < 9 ) printInstructions( argv[0] );

    // Check to make sure the first two arguments are positive numbers, and that
    // the last four arguments are positive or negative numbers.
    if(     isdigit(argv[1][0]) &&
            isdigit(argv[2][0]) &&
        (   isdigit(argv[3][0]) || ( (argv[3][0]=='-') && isdigit(argv[3][1]))   ) &&
        (   isdigit(argv[4][0]) || ( (argv[4][0]=='-') && isdigit(argv[4][1]))   ) &&
        (   isdigit(argv[5][0]) || ( (argv[5][0]=='-') && isdigit(argv[5][1]))   ) &&
        (   isdigit(argv[6][0]) || ( (argv[6][0]=='-') && isdigit(argv[6][1]))   ) &&
        (   isdigit(argv[7][0]) || ( (argv[7][0]=='-') && isdigit(argv[7][1]))   ) &&
        (   isdigit(argv[8][0]) || ( (argv[8][0]=='-') && isdigit(argv[8][1]))   )       )
	{
        // set bond length
        bondLength = atof( argv[1] );
        // set C-C separation distance
        r = atof( argv[2] );

        // set angle of rotation about the z-axis for molecule 1
        angleZ1 = atof( argv[3] );
        // set angle of rotation about the x-axis for molecule 1
        angleX1 = atof( argv[4] );
        // set angle of rotation about the Y-axis for molecule 1
        angleY1 = atof( argv[5] );
        // ditto for molecule 2
        angleZ2 = atof( argv[6] );
        angleX2 = atof( argv[7] );
        angleY2 = atof( argv[8] );

        if(   argc>=10  &&  (strcmp(argv[9],"-mol")==0)   )
            mode = molpro;
    }
    else
        printInstructions( argv[0] );


    const int C = 6;
    const int H = 1;


    // Construct identical methane molecule, labeling each atom uniquely
    Molecule m1, m2;

    // Add carbon at origin
    m1.add_atom( C, 0.0, 0.0, 0.0, 1 );
    m2.add_atom( C, 0.0, 0.0, 0.0, 6 );

    // Add a hydrogen on +x-axis
    m1.add_atom( H,   bondLength, 0.0, 0.0, 2 );
    m2.add_atom( H,   bondLength, 0.0, 0.0, 7 );

    double interiorAngle = (double)2.0 * atan(sqrt((double)2.0));
    double X = bondLength*cos(interiorAngle);            // x coordinate for all other H atoms
    double YZBondLength = bondLength*sin(interiorAngle); // Projected 2-D distance from origin (on the YZ-plane)

    // add a hydrogen opposite the YZ plane from the first, at the 12:00 position
    m1.add_atom( H, X, 0.0, YZBondLength, 3 );
    m2.add_atom( H, X, 0.0, YZBondLength, 8 );
    // add a hydrogen opposite the YZ plane from the first, at the  8:00 position
    m1.add_atom( H, X, YZBondLength*cos((double)210.0*M_PI/(double)180.0), YZBondLength*sin((double)210.0*M_PI/(double)180.0),  4 );
    m2.add_atom( H, X, YZBondLength*cos((double)210.0*M_PI/(double)180.0), YZBondLength*sin((double)210.0*M_PI/(double)180.0),  9 );
    // add a hydrogen opposite the YZ plane from the first, at the  4:00 position
    m1.add_atom( H, X, YZBondLength*cos((double)330.0*M_PI/(double)180.0), YZBondLength*sin((double)330.0*M_PI/(double)180.0),  5 );
    m2.add_atom( H, X, YZBondLength*cos((double)330.0*M_PI/(double)180.0), YZBondLength*sin((double)330.0*M_PI/(double)180.0), 10 );


    // Separate the molecules, along the x-axis
    m1.translate( -r/(double)2.0, 0, 0 );
    m2.translate(  r/(double)2.0, 0, 0 );
	
    
    // Apply rotations

    m1.rotateZ( angleZ1 );
    m1.rotateX( angleX1 );
    m1.rotateY( angleY1 );

    m2.rotateZ( angleZ2 );
    m2.rotateX( angleX2 );
    m2.rotateY( angleY2 );


    
     // OUTPUT DATA TO SCREEN
    /////////////////////////

    // molpro header
    if( mode==molpro )
        cout << "geometry={" << endl;

    // Atom count & comment line
    cout << m1.atom_count() + m2.atom_count() << endl;
    cout << "CH4 Dimer [Bond Length: " << bondLength << " Angstroms][C-C Separation: " << r << " Angstroms]" << endl;

    // molecular geometry
    if( mode==molpro )
    {
        m1.displayMolpro();
        m2.displayMolpro();
    }
    else if( mode==xyz )
    {
        m1.displayXYZ();
        m2.displayXYZ();
    }

    // molpro closing brace
    if( mode==molpro)
        cout << "}" << endl;


   
    return (EXIT_SUCCESS);
}


void printInstructions( char *executablePath )
{
    char * executableName = executablePath;
    char * i = executablePath;
    #ifndef _WIN32
    const char file_separator = '/';
    #else
    const char file_separator = '\\';
    #endif

    // Skip to the last separator in the file
    while( *i != '\0' )
    {
        i++;
        if( (*i) == file_separator )
            executableName = i+1;
    }


    cout << "\nUSAGE: " << executableName << " "
        << "<b> <r> <rZ1> <rX1> <rY1> <rZ2> <rX2> <rY2> [-mol | -xyz]\n\n"
        << "b    - bond length (in Angstroms)\n"
        << "r    - distance between carbon centers (in Angstroms)\n"
        << "rZ1  - rotation about the Z-axis for the 1st molecule (in degrees)\n"
        << "rX1  - rotation about the X-axis for the 1st molecule (in degrees)\n"
        << "rY1  - rotation about the Y-axis for the 1st molecule (in degrees)\n"
        << "rZ2  - rotation about the Z-axis for the 2nd molecule (in degrees)\n"
        << "rX2  - rotation about the X-axis for the 2nd molecule (in degrees)\n"
        << "rY2  - rotation about the Y-axis for the 2nd molecule (in degrees)\n"
        << "-mol - add this option if you want geometry formatted for Molpro (default is xyz)\n"
        << "-xyz - add this option if you want geometry formatted for an .xyz file\n\n"
        << "NOTE: Rotations are applied in Z-X-Y order.\n"
        << endl;

    exit(1);
}
