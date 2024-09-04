#ifndef GRID_H
#define GRID_H


#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#include "Base/Constants.h"
#include "Base/Facilitator.h"
#include "Base/Matrix.h"
#include "Base/IntRandom.h"
#include "Base/DynArrays.h"
#include "Base/Estring.h"
#include "Base/Position.h"

#ifndef PatternDim
#define PatternDim  180
#endif


class Component;
class GridSite;
class GridWindow;
class Species;
class CharString;


extern int SymmetryPattern[PatternDim][3];
extern int PatternsNormalHollow[289];
extern int PatternsNormalBridge[289];
extern int PatternsNormalAtop[289];
extern Component* theNULLSpecies;
enum SiteCoordination{NONE = 0, Atop=1, Bridge=2, Hollow_3 = 3, Hollow=4, Hollow_6 = 6};

typedef struct {
  Eshort x,y;
} Patt;

// void Redisplay();
float FindDistance(GridSite* MySite, GridSite* NeighSite);
float FindDistance(GridSite* MySite, int index);

class Grid : public Facilitator
{

public:
    
	fstream gin;
    
	bool MainGrid;  // DWDWDW
    
    //Constructor
	Grid(fstream &fin, fstream &fout, int SetSize=0); 
  bool AllSite;
  bool Atop_a;
  bool Bridge_a;
  bool Hollow_a;
	Grid(fstream &fin, fstream &fout, int SetSize, double Cut,
		CharString &M, CharString &Face);

	// Initialize Grid
	void Init();
	// Read Input for Startup
	void Read(fstream &fin);
	void ReOpenImages(fstream &fin);
	// Resize
	void Resize(int size);

	void ShowAllOrientations();
	// COORDINATION : Finding Surface Coordination
	SiteCoordination FindCoordination(int i, int j);

	// DEPOPULATE : Depopulate the surface grid using fragments
	int DePopulate();

	double GetCoverage(Component*);
	//void WriteGrid();
	//void RetrieveGrid();
	//void Restart();
	// Pick A Random Site
	GridSite* PickRandom(Component *Type);
	GridSite* PickRandom();
  GridSite* PickRandomAtop();
	void PickRandomPair(GridSite* &A, GridSite* &B);

	// Change Surface
	void ChangeSurface(CharString &Face);
	void SetSurface(CharString TypeSurface);
	CharString& GetSurface() {return SurfaceType;}

	void GetSurfaceAtoms(int* Coverage, int &Total);
	void AddMetalAtom(GridSite* PassedSite);

	void PutDownMetalAtoms();
	void ConnectMetalAndOverlayer();
	void SetSiteMatrix();

	GridSite* ColorSite[2];
	
	// Allows access to grid private parts
	double SiteArea;
	double SurfaceCoordination;	

	int GetNumMetalAtoms() { return NumberOfMetalAtoms; }

	GridSite* LocateCentralSite(int Coord);
  GridSite* LocateDistinguishSites(CharString Coord);
	GridSite* LocateSurroundingSite(GridSite* AnySite, int Pos);

	double Cutoff;
	double GridLength[2];

	Component* Alloy_Type;
	double Alloy_Composition;

	int SubstrateDimension;
	int SubstratePattern[12][12];
	
	int AdsorbateDimension;
	int AdsorbatePattern[24][24];
	void AddOverLayerAdsorbates();

	int MetalLayer;
	int AdsorbateLayer;
    double MMdistance;
	int NumberOfMetalAtoms;
	int NumberOfBonds;
	CharString SurfaceType;
	IRandom RandomNumSurr;
	IRandom RandomNumX;
	IRandom RandomNumY;

	float SiteDistances[PatternDim][7];

	// Metal Type
	Component* Metal;

	// Metal Dimension
	int EDim[2];
	int OldEDim[2];

	double DeltaCoverage;

	// For 111 Surface, Searching patterns
    Patt BridgeTypeOne[300];
    Patt BridgeTypeTwo[300];
    Patt SNormalPattern[300];
    Patt SNormalPatternRotate180[300];
	void WriteCoordinates(int a);
  void WriteCoordinatesXYZ(int a);
 void WriteMetalCoordinates(GridSite *a);
 void WriteCoordinatesHan();
	// Surface Sites and their connections
	Matrix<GridSite*> Surface[4];
	void DefineQuadrants();
    int Quadrant[180][4];

    //MPMPMPMPMPMP
    string coords_file;
    
    int new_surface;

// private:
	int Coordinates_init;
  char Coordinates_File[100];
  char CoordinatesXYZ_File[100];
  char CoordinatesHan_File[100];
	int MetalCoordinates_init;
};

extern Grid *theGrid;

#endif		
