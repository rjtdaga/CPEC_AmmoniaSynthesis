// GridSite.h
// ----------



#ifndef ATOM_H
#define	ATOM_H


#include "Base/Position.h"
#include "Base/DynArrays.h"
#include "Base/Constants.h"
#include "Grid.h"
#include "Species.h"
class Component;

// Should be a pairwise choosing of sites

class GridSite
{
public:

	// CONSTRUCTORS

	GridSite();
	GridSite(GridSite &A); // copy constructor	
	GridSite(const EPosition<double> &at,Component* A); 
	GridSite(const double &x, const double &y, const double &z, Component* A);
 
  double Ads_Pos[100][3]; // Defines the position of the adsorbate on the gridsite
  CharString Ads_Atom[100]; // Defines the names of the atoms of the adsorbate on the gridsite
  int NumAtoms; // Number of atoms in the adsorbate molecule
  int SiteType;
	// SET POSITION : Mutator to Set Position
	void SetPosition(const EPosition<double> &at);
	void SetPosition(const double &x, const double &y, const double &z);
 
	// SET TYPE : Mutator to Set Type
	void SetType(Component* A);

	// GET TYPE : Inspector to Atom Type
	Component* GetType() const {return Type;} 

	// GET POSITION : Inspector to Atom Position
	EPosition<double> GetPosition() const;

	// SET ATTRIBUTES : Assignment Operator for Atom Attributes
	GridSite& operator=(GridSite &A);
	GridSite& operator=(Component* A);
	GridSite& operator=(const EPosition<double> &A);

	// BONDS : Metal - Surface Bonds
    GridSite* Bonds[13];
	Eshort BondSize;
 
 int Coordination;
 CharString Name[100];

	// BONDS : Metal - Metal Bonds
	SeqList<GridSite*> MetalBonds;

	// BONDS : Far Surrouding Bonds to overlayer
	GridSite *Surround[PatternDim];
	// Eshort SurroundTalkBack[PatternDim];
	Eshort SurroundSize;

	// ADD BOND : Mutator to add a bond
	void AddBond(GridSite* New);
	void AddMetalBond(GridSite* New);

	// POSITION : Spatial location of this type 
	EPosition<double> Location;
	double FindDistance(GridSite *Away, double GridLength[]);
	void FindUnitVector(GridSite *Away, double GridLength[], double Delta[]);

	Eshort GridPos[3];
	void SetGridPos(Eshort x, Eshort y, Eshort z);
	void GetGridPos(Eshort NewPos[]);
    
    
    Component* Type;
	Eshort OrientationPos;
	Eshort Orientation;
    

	// CURRENT TYPE : occupancy of this grid site
	double Angle;
  //double BindingEnergy;
private:

};

extern GridSite* thisSite;

#endif
