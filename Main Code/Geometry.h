// Geometry.h
// ----------
// Header for the Geometry file.  The Geometry files read in the
// GEOMETRY section of the input file.


#ifndef Geometry_h
#define Geometry_h



#include <fstream>
using namespace std;



#include "Base/DynArrays.h"
// #include "Base/Facilitator.h"
#include "Base/Estring.h"
#include "Base/Position.h"
#include "Base/Constants.h"


// This is the maximum number of atoms a fragment can have.
// I think Eric's "12" shows up only in Geometry.h/.cpp, but it may
// also show up in other places as well, such as Model.cpp
#define MAX_ATOMS  30


typedef struct
{
	EPosition<double> Position[MAX_ATOMS];
	Eshort NumAtoms;
	Eshort AtomNames[MAX_ATOMS];
	Eshort SpeciesName;
	Eshort MetalCoord;
  CharString SpecName;
  CharString AtomName[MAX_ATOMS];
	Eshort MMFF94_Type[MAX_ATOMS];
	float Charge[MAX_ATOMS];
  double SurfaceCharge;
} Coordinates;

class Species;
class Geometry : public Facilitator
{
public:

    // The constructor requires you to pass in a species
    Geometry(Species &theSpecies, fstream &fin, fstream &fout);
    
    // This reads the input file
    int Read(fstream &fin);

private:
    // These two are responsible for reorienting the input geometry
    // into a more standard orientation.
    void ReOrient(EPosition<double> Basis[], EPosition<double> Coords[],
        Eshort Num);    
    void rotation_matrix_calc(double va[], double vb[], double basvec[3][3],
        double t[][3]);

public:
    // theCoordinates is referenced by other files, so must be public
    SeqList <Coordinates> theCoordinates;
    int NumberOfGeometries;
    
private:
    Species &theSpecies;

};
#endif

