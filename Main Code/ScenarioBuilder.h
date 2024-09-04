
#include "Base/Facilitator.h"
#include "Base/Constants.h"
#include "Grid.h"

class Scenario;
class Component;
class ElementaryRxn;
#ifndef ScenarioBuilder_h
#define ScenarioBuilder_h

class ScenarioBuilder : public Facilitator
{
public:

    int Calculate(Scenario *S, Scenario* FoundScenarios[]);
    void CalculatePairwisePreferredDistance();
    
    ScenarioBuilder(fstream &fin, fstream &fout);
    
    double SurfaceRate(double, double);
    
//    void FindBestSingleSites(GridSite *AnalSite,
//        Component* ProductA, Component* ProductB,
//        GridSite* BestSiteA[], double BestEnergyA[],
//        GridSite* BestSiteB[], double BestEnergyB[],
//        double cutoff, bool Forced); // DWDWDW
        
    bool FindBestSingleSites(GridSite * AnalSite,
    Component * ProductA, Component * ProductB, GridSite * BestSiteA,
    double BestEnergyA, GridSite * BestSiteB,
    double BestEnergyB, double cutoff, bool Forced); 
    
    bool FindBestSingleSite(GridSite *AnalSite,
        Component* ProductA, 
        GridSite* &BestSiteA, double &BestEnergyA,
        double cutoff, bool Forced); // DWDWDW
    
    int CalculateRxnScenarios(GridSite* AnalSite, ElementaryRxn* RxnToConsider,
        GridSite* ReactSiteA[], double ReactSiteA_E[],
        GridSite* ProdSiteA[],  double ProdSiteA_E[],
        GridSite* ReactSiteB[], double ReactSiteB_E[],
        GridSite* ProdSiteB[],  double ProdSiteB_E[]);
    
//    GridSite* FindBestPairBasedSite(GridSite* Center, /* Center Site */
//        GridSite* SiteA /* Given to watch effect of local environment */,
//        int qa /* quadrant to exclude from consideration */,
//        Component* NeighType, /* set surroundings to */
//        double &BestBEnergy, double &BestAEnergy);
        
    GridSite *FindBestPairBasedSite(GridSite * Center,     /* Center Site */
    GridSite * SiteA /* Given to watch effect of local environment */ ,
    Component * NeighType,      /* set surroundings to */
    double &BestBEnergy, double &BestAEnergy);
    
    double CalculateRate(ElementaryRxn* RxnToConsider, double RA, double RB,
        double PA,  double PB, double &Barrier);
    
    int QuadrantIndex;
};

extern ScenarioBuilder *theScenarioBuilder;

#endif
