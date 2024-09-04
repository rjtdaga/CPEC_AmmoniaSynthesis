#ifndef Simulation_h
#define Simulation_h

#include <sys/types.h>
#include "Component.h"
#include "Base/Constants.h"
#include "Base/Facilitator.h"
#include "Scenario.h"
#include "Model.h"
#include "Base/DynArrays.h"

class GridSite;
class Grid;
class Reactions;
class Component;
class ElementaryRxn;
#ifndef MAXN
#define MAXN 30
#endif

class Simulation : public Facilitator
{
public:
    Simulation(fstream &fin, fstream &fout);
    
    void Read(fstream &fin);			// Initialization of parameters
    void Read_Analysis(fstream &fin);
    void FindTime();
    bool Solvent;
    bool Gas;
    
    // For multi-level equilibrated reactions
    
    
    void SimulationSummary();
    bool UnknownBE_BOC;
    void Simulate(fstream &fout);		// Simulate Surface Evolution
    
    void Startup(fstream &fout);
    
    void ConsistentBarrierModification(fstream &fout);
    
    void FindUnknownBE_BOC();
    
    ElementaryRxn* PickAReaction(); 
    
    ElementaryRxn* PickEquilReaction(); 
    
    bool UserSet;
    
    double StableSiteEnergy(Component * PassedComponent);
    
    bool NoInteraction;
    
    void Scale();
    
    int Univ;
    
    void CalculateStaticVariables();
    
    int Initialize_Matrices();
    
    int Initialize_Equil_Matrices();
    
    void WriteEmptySites();
    
    int UpdateEnvironment();
    
    void Define_SimRxns();

    void CalculateSumOfRates(int &n);
    
    void VTS_KineticMC(int numsteps=1);
    
    void EquilSteps(int numsteps = 1);

    void CalculateVTSReactionRates();
    
    void CalculateEquilReactionRates();
    
    void AnalyzeReactions(double);
    
    int Populate(int *Diff, fstream &fout); // Populate Surface

    int SuperFastDiffuseEquilibration(double num=1.);
    
    void WriteCoverage();
    
    void WriteDesorbed();
    
    void WriteReactions();
    
    void WriteBindingEnergies();
    
    void WriteBarrier();
    
    ElementaryRxn* PickRandomReaction();
    
    ElementaryRxn* ReactionPossible(ElementaryRxn* );
    
    ElementaryRxn* throttleRxn;
    
    double GetTemperature() { return Temperature; }
    
    double GetTime() { return Time; } 
    
    double ChargePattern[10];
    
    double TimePeriod[10];
    
    int NumberChargeType;    
    
    int FindChargeIndex(double SC);
    
    double SurfaceCharge;
    
    bool Ads_Arrhenius;
    bool Ads_CollisionTheory;
    
private:
    double TotalReactionRate;
    
    double TotalEquilReactionRate;

    bool Population; // DWDWDW
    
    int RemoveGas;   // Boolean to remove gas after surface init  
    
    int ConstConc; // Boolean to keep the concentration of solvent species constant
    
    double StopTime; // Time to Stop Program
    
    double DelayTime; // Time to equilbrate through normal simulation
                      // then to remove partial pressures of gas
                      // components 
    
    double Temperature;		// Surface Temperature
    double Beta;			// Heating Rate
    double Time;			// Simulation Time
    int GasOff;


    double running_time_step;
    double input_time_step;

    double TimeStep;
    //----------------------------------------------------------------------------
    // Variables specifically for Dybeck parameters to be taken from input file:
    int CycleTime;
    int MinExecSuperbasin;
    int MaxRxnQueueSize;
    int QEMaxDiff;
    int k_buff;
    
    //---------------------------------------------------------------------------
    double equilibration_factor; // A factor by which to increase
                                 // equilibration events
    
    int Coverage_Init;
    
    int Desorbed_Init;
    char Desorbed_File[100];
    int Reaction_Init;
    char Reaction_File[100];
    int BindingEnergy_Init;
    char BindingEnergy_File[100];
    char EmptySites_File[100];
    int Barrier_Init;
    char Barrier_File[100];

    
    int NumberOfSimReactions;
    int NumberOfEquilRxns;
    ElementaryRxn **SimReactions;
    ElementaryRxn **EquilReactions;
    
    bool Debug;
    
    double* Barrier;
    double* Prefactor;
    void unThrottleAllReactions();
    double TotalSystemEnergy;
    double RateTotal;
    int ScenarioTotalThisPass;
    SeqList<ScenarioAction> All[100];  
    SeqList<ScenarioAction> EquilAll[100];  
    
    double VTS_TimeStep;
    double FTS_TimeStep;
    double VTS_Chance;
    double *EventsByRate;
    double RT;
    ScenarioAction SaveScenario;
    
    char Coverage_File[100];
};

extern Simulation *theSimulation;	// Global Pointer to this Simulation

extern int NumberOfSpecies;

#endif

