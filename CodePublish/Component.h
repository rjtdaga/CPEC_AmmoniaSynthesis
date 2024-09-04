#ifndef component_h
#define component_h



#include <fstream>
#include <cmath>
using namespace std;



#include "Base/Constants.h"
#include "Base/Estring.h"
#include "Base/Color.h"

class ElementaryRxn;



enum Orientation{Parallel = 2, Perpendicular = 1, NoOrientation = 0};

#ifndef MAX_NUM_SPECIES
#define MAX_NUM_SPECIES 25
#endif

class Component {

public:

    Component();
    CharString& 		GetName() {return SpeciesName;}; 
    void        		SetName(CharString Name);
    void 			addFragments(Component *StrongE,
        Component *WeakE);

    Component& operator=(Component A);
    // Attributes
    double Charge;
    double GetAtomizationEnergy() const {return Atomize;}
    double GetVDW()   const {return VDWRadius;} 
    double GetCore()  const {return CoreRadius;}	
    double GetMolWt() const {return MolWt;}
    void   SetMolWt(const double val);
    double* SetCore(const double val);
    double* SetVDW( const double val); // memory address pointer
    double* SetAtomizationEnergy(const double val);
    int  GetClosedShell() const {return ClosedShell;}
    void  SetClosedShell(const int val) {ClosedShell = val; return;}
    Color GetColor() { return MyColor; }
    void SetColor(double r, double g, double b);
    void OutputColor(double r, double g, double b);
    Orientation GetOrientation() const { return theOrientation; }
    void SetOrientation(CharString A);
    void SetElectrons(int a) { Electrons = a; }
    int  GetElectrons() {return Electrons;} 
    int CrossReference();
    double subAtomCoordinates[10];

    Eshort _id; // local copy
    Eshort _idhook; // hook for species listing only on surface 

    Component* GetBindingAtom();
    int isAtom();
    int Atom;

    Component* Strong;
    Component* Weak;

	// Simulation Variables 
    double SimEnergy[10];
    double SimBindingEnergy;
    double SimBindingEnergyCount;
    double SimVacancyEnergyCount[7];
    double SimVacancyEnergy;
    double VacancyCount[7];

    double NumberCount;

    double StableSiteEnergy[10];
    double BindingEnergy[7][2][10];

    double* SetPressure(double a) {Pressure = a; return &Pressure;}
    double  GetPressure() { return Pressure; }
    double* SetConcentration(double a) {Concentration = a; return &Concentration;}
    double  GetConcentration() { return Concentration; }

    double CalculateAdsorptionRate(double Temp);
    double TotalDesorbed;

    
    bool OnSurface; // DWDWDW


    Orientation   theOrientation;
    CharString 	SpeciesName;

    // Macroscopics
    int Electrons;

    double Pressure;
    double Concentration;
    int NumberOfAtoms;
    int NumAtoms;
    int EquilCov;
    CharString AtomName[100];
    // Attributes

    Color MyColor;
    double VDWRadius;		// VDW Radius
    double MolWt;			// Molecular Weight
    double CoreRadius;		// Core Radius
    int ClosedShell;		// Radical or Closed Shell

    ElementaryRxn* DiffusionReaction;
    Eshort Consider;
    float FavorableDist[MAX_NUM_SPECIES][7][10];

private:    
    double Atomize;  		// Energy to Atomize
};

extern Component* theNULLSpecies;
extern Component* thisComponent;
 
#endif

