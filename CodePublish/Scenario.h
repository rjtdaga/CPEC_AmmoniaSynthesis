#include "Base/Constants.h"
#include "Base/Random.h"
#include "GridSite.h"
#include "Grid.h"
#ifndef Scenario_h
#define Scenario_h
#include "ElementaryRxn.h" 

class Component;
class ScenarioAction;

#ifndef SITE_NUM_SCENARIOS
#define SITE_NUM_SCENARIOS 100
#endif

class Scenario
{
    
public:

    GridSite *ReactantSites[2];
	GridSite *ProductSites[2];
    Component *Reactants[2];
    Component *Products[2];
    bool Possible; // DWDWDW
    double ForwardRate;
    double ReverseRate;
    double ForwardBarrier;
	double ForwardPreexponential;
    double ReverseBarrier;
    double Preexponential;
	double ReversePreexponential;
	Eshort OrigReactOrientation[2];
    
    // DW: Is this ever even used??
    bool GetPossibility(int Surr, int coord); // DWDWDW
    
	bool Screen(GridSite *PassedSite); // DWDWDW
    
	void SetUpScenario(const ScenarioAction &S);	
	Scenario();
	 //~Scenario();
	Scenario(const Scenario &S);

	virtual void doit();
	virtual void undoit();

	int Equilibrate();

	void SetReaction(ElementaryRxn *thisRxn); 
	Process GetProcess();  
    ElementaryRxn *theReaction;
virtual ~Scenario();

};

class ScenarioAction {
	public:

	ScenarioAction() {}; // constructor
	ScenarioAction(GridSite* R[2], GridSite* P[2], Component* r[2],
			Component* p[2], double Rate) {
	  ReactantSites[0] = R[0];
	  ReactantSites[1] = R[1];
	  ProductSites[0] = P[0];
	  ProductSites[1] = P[1];
	  ForwardRate = Rate;
	} 
	ScenarioAction(const ScenarioAction &SA) { // copy constructor
	  ReactantSites[0] = SA.ReactantSites[0];
	  ReactantSites[1] = SA.ReactantSites[1];
	  ProductSites[0] = SA.ProductSites[0];
	  ProductSites[1] = SA.ProductSites[1];
	  Rxn = SA.Rxn;
	  ForwardRate = SA.ForwardRate; 
	}
	ScenarioAction& operator=(const ScenarioAction &SA) {
	  ReactantSites[0] = SA.ReactantSites[0];
	  ReactantSites[1] = SA.ReactantSites[1];
	  ProductSites[0] = SA.ProductSites[0];
	  ProductSites[1] = SA.ProductSites[1];
	  ForwardRate = SA.ForwardRate;
	  Rxn = SA.Rxn;
	  return *this;
	} 
	void SetUpScenarioAction(const Scenario *SA) {
	  ReactantSites[0] = SA->ReactantSites[0];
	  ReactantSites[1] = SA->ReactantSites[1];
	  ProductSites[0] = SA->ProductSites[0];
	  ProductSites[1] = SA->ProductSites[1];
	  Rxn = SA->theReaction;
	  ForwardRate = SA->ForwardRate;
	}
        GridSite *ReactantSites[2];
	GridSite *ProductSites[2];
	ElementaryRxn* Rxn;
        double ForwardRate;
};
extern Scenario* thisScenario;
extern Scenario* thisEquilScenario;
#endif

