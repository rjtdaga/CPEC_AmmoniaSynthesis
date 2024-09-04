#include <cmath>
#include <cassert>
#include <ctime>
#include <string>
#include "rank.h"
#include <malloc.h>
#include "Simulation.h"
#include "ScenarioBuilder.h"
#include "Scenario.h"
#include "Reactions.h"
#include "ElementaryRxn.h"
#include "Species.h"
#include "Model.h"
#include "Grid.h"
#include "Component.h"
#include "Base/Random.h"
#include "Base/IntRandom.h"
#include "Base/Matrix.h"
#include "Interactions.h"
#include "ShortRoutines.h"
#include "Base/CommonBlock.h"
#include "Base/DynArrays.h"
#include "Base/Constants.h"
#include "SystemHeaders.h"

Simulation *theSimulation = NULL;
double double_limit = 1.e307;
GridSite **Adsorbates;
int SimulationStep = 0;
ElementaryRxn *throttleRxn = NULL;
double Tao = 0.0;
bool nonQEEventThisCycle = 0;
bool isnonQEexecuted = false;
ElementaryRxn *lastNonQEReaction = NULL;
int SwitchCharge = 0;
int CycleStep = 0;
SeqList <int> dummy[MAX_NUM_REACTIONS];
double CycleTimeLog[1000];
int SearchEvent = 0;
time_t TrackTime;
time_t OldTime;
time_t NetTime;
double PeriodTime = 0.;
int oldChargeIndex = 0;
int ChargeIndex = 0;
Simulation::Simulation(fstream &fin, fstream &fout) : Facilitator(fin, fout)
{
    Adsorbates = new GridSite *[theGrid->EDim[0] * theGrid->EDim[1] + 1];
    Population = false; // DWDWDW
    thisScenario = new Scenario();
    thisEquilScenario = new Scenario();
    theSimulation = this;
    equilibration_factor = 1.;
    DelayTime = 0.;
    Time = 10. / double_limit;
    StopTime = 3600.0;           // must be before read, default time
    Temperature = 300.;
    UserSet = false; 
    NoInteraction = false;
    TimeStep = 0.;
    running_time_step = 1.0;
    RT = R_GasConst * Temperature;
    Beta = 0.;
    Debug = false; // DWDWDW
    Initialize((CharString) "Simulation", fin);
    Solvent = false;
    Gas = false;
    Ads_CollisionTheory = true;
    Ads_Arrhenius = false;
    UnknownBE_BOC = 0;
    //-------------------------------------------
    // Initialization of Dybeck Parameters
    CycleTime = 20;
    MinExecSuperbasin = 20;
    MaxRxnQueueSize = 200;
    QEMaxDiff = 20;
    k_buff = 1000;
    //---------------------------------------------
    // Reading the input file
    Read(fin);
    //-----------------------------------------------
    if (UnknownBE_BOC == 1)
    {
      FindUnknownBE_BOC();
//      theSpecies->theGeometries->theCoordinates.ResetToFront();
//      while (theSpecies->theGeometries->theCoordinates)
//      {
//        Coordinates SiteCoord = theSpecies->theGeometries->theCoordinates();
//        cout << SiteCoord.SpecName << " " << SiteCoord.MetalCoord << " " << SiteCoord.SurfaceCharge << endl;
//        ++(theSpecies->theGeometries->theCoordinates);
//      }
    }
    //exit(1);
    SurfaceCharge = ChargePattern[0];
    RemoveGas = true; // DWDWDW
    Coverage_Init = true; // DWDWDW
    Desorbed_Init = true; // DWDWDW
    Reaction_Init = true; // DWDWDW
    BindingEnergy_Init = true; // DWDWDW
    Barrier_Init = true; // DWDWDW
    Univ = 0;
    const char *d1 = "SurfaceCoverage";
    strcpy(Coverage_File, d1);
    strcat(Coverage_File, rnk);
    const char *d2 = "Desorbed";
    strcpy(Desorbed_File, d2);
    strcat(Desorbed_File, rnk);
    const char* d3 = "ReactionRates";
    strcpy(Reaction_File, d3);
    strcat(Reaction_File, rnk);
    const char *d4 = "BindingEnergies";
    strcpy(BindingEnergy_File, d4);
    strcat(BindingEnergy_File, rnk);
    const char *d5 = "Activation_Energies";
    strcpy(Barrier_File, d5);
    strcat(Barrier_File, rnk);
    const char *d6 = "EmptySites";
    strcpy(EmptySites_File, d6);
    strcat(EmptySites_File, rnk);
    GasOff = false; // DWDWDW
    Barrier = new double[Rxns->Size()];
    Prefactor = new double[Rxns->Size()];
    //cout << "Hi" << endl;
    SimReactions = new ElementaryRxn *[Rxns->Size()];
    EquilReactions = new ElementaryRxn *[Rxns->Size()];
    NumberOfSimReactions = 0;
    NumberOfEquilRxns = 0;
    Rxns->ResetToFront();
    while(*Rxns)
    {
        thisReaction = Rxns->GetPtr();
        //cout << thisReaction->GetName() << " " << thisReaction->_id << endl;
        if(thisReaction->GetProcess() != Diff && thisReaction->SurfaceCharge == SurfaceCharge && thisReaction->UserEquil == 0)
        {
            SimReactions[NumberOfSimReactions] = thisReaction;
            SimReactions[NumberOfSimReactions]->_idhook = NumberOfSimReactions;
            ++NumberOfSimReactions;
            //cout << thisReaction->GetName() << " " << NumberOfSimReactions-1 << " " << thisReaction->GetProcessName() << endl;
        }
        else if(thisReaction->GetProcess() != Diff && thisReaction->SurfaceCharge == SurfaceCharge && thisReaction->UserEquil == 1)
        {
            EquilReactions[NumberOfEquilRxns] = thisReaction;
            EquilReactions[NumberOfEquilRxns]->_idhook = NumberOfEquilRxns;
            ++NumberOfEquilRxns;
        }
        ++(*Rxns);
    }
    //exit(1);
    Initialize((CharString) "Analysis", fin);
    Read_Analysis(fin);
    if (TimePeriod[0] != 0)
    {
      for (int i = 0; i < NumberChargeType; ++i)
      {
        PeriodTime = PeriodTime + TimePeriod[i];
        //cout << TimePeriod[i] << " " << rank << endl;
      }
    }
    else
    {
      cout << "Period time is not mentioned!!" << endl;
    }
    EventsByRate = new double[NumberOfSimReactions];
    for(int i = 0; i < NumberOfSimReactions; ++i)
    {
        EventsByRate[i] = 0.;
    }
//    cout << "Hi" << endl;
//    GridSite *Central = theGrid->LocateDistinguishSites("Hollow");
//    Eshort a[3];
//    Central->SetType(Listing[0]);
//    Central->GetGridPos(a);
//    cout << theModel->CalculateBindingEnergy(Central) << endl;
//    exit(1);
    if (NumberOfSimReactions >= MAX_NUM_REACTIONS)
    {
      cout << "Number of simulation reaction exceeds Max number of reactions" << endl;
      exit(1);
    }
    ConsistentBarrierModification(fout);
    theModel->CheckInteractions();
    //--------------------------------------------------
//    int Layer = theGrid->AdsorbateLayer;
//    for (int i = 0; i < theGrid->EDim[0]; ++i)
//    {
//      for (int j = 0; j < theGrid->EDim[1]; ++j)
//      {
//        GridSite * thisSite = theGrid->Surface[Layer][i][j];
//        if (thisSite != NULL && thisSite->BondSize == 3)
//        {
//          thisSite->SetType(Listing[2]);
//          theModel->CalculateBindingEnergy(thisSite, true);
//          Listing[2]->NumberOfAtoms++;
//        }
//      }
//    }
//    theGrid->WriteCoordinatesXYZ(0);
//    exit(1);
    //--------------------------------------------------
    theScenarioBuilder->CalculatePairwisePreferredDistance();
    theGrid->AddOverLayerAdsorbates();
//    Rxns->ResetToFront();
//    while (*Rxns)
//    {
//      ElementaryRxn *thisReaction = Rxns->GetPtr();
//      cout << thisReaction->GetName() << " " << thisReaction->deltaH << " " << thisReaction->ActivationEnergy << endl;
//      ++(*Rxns);
//    }
//    exit(1);
    // Checking whether the surrounding sites are correct or not
    
}


bool RelaxProductState(GridSite * Products[]); // DWDWDW

void Simulation::Define_SimRxns()
{
  Rxns->ResetToFront();
  NumberOfSimReactions = 0;
  NumberOfEquilRxns = 0;
  ElementaryRxn * dummy; 
  while(*Rxns)
  {
    dummy = Rxns->GetPtr();
    if (dummy->SurfaceCharge == SurfaceCharge && dummy->GetProcess() != Diff && dummy->UserEquil == 0)
    {
      SimReactions[NumberOfSimReactions] = dummy;
      SimReactions[NumberOfSimReactions]->_idhook = NumberOfSimReactions;
      ++NumberOfSimReactions;
    }
    if (dummy->SurfaceCharge == SurfaceCharge && dummy->GetProcess() != Diff && dummy->UserEquil == 1)
    {
      EquilReactions[NumberOfEquilRxns] = dummy;
      EquilReactions[NumberOfEquilRxns]->_idhook = NumberOfEquilRxns;
      ++NumberOfEquilRxns;
    }
    ++(*Rxns);
  }
}

void Simulation::FindUnknownBE_BOC()
{
  for (int z = 0; z < NumberOfSpecies; ++z)
  {
    Component * A = Listing[z];
    for (int j = 0; j < NumberChargeType; ++j)
    {
      int non_zero[7];
      int MaxInd = 0;
      int n = 0;
      for (int i = 1; i < 7; ++i)
      {
        double BE = A->BindingEnergy[i][0][j];
        if (BE > 0. && i!=5)
        {
          non_zero[i] = 1;
          ++n;
          if (MaxInd < i)
          {
            MaxInd = i;
          }
        }
        else
        {
          non_zero[i] = 0;
        }
      }
      double AvgE = 0.;
      for (int i = 1; i <= MaxInd; ++i)
      {
        if (non_zero[i] == 1)
        {
          AvgE += A->BindingEnergy[i][0][j]/(2-1.0/i)/n;
        }
      }
      for (int i = 1; i < MaxInd; ++i)
      {
        if (non_zero[i] == 0)
        {
          A->BindingEnergy[i][0][j] = AvgE*(2-1.0/i);
          Coordinates SiteCoord;
          Coordinates Dummy;
          bool Found = 0;
          theSpecies->theGeometries->theCoordinates.ResetToFront();
          while(theSpecies->theGeometries->theCoordinates)
          {
            SiteCoord = theSpecies->theGeometries->theCoordinates();
            if (SiteCoord.SurfaceCharge == ChargePattern[j] && A->_id == SiteCoord.SpeciesName && SiteCoord.MetalCoord == i)
            {
              Found == 1;
            }
            if (SiteCoord.SurfaceCharge == ChargePattern[j] && A->_id == SiteCoord.SpeciesName)
            {
              Dummy = SiteCoord;
            }
            ++(theSpecies->theGeometries->theCoordinates);
          }
          if (Found == 0)
          {
            Dummy.MetalCoord = i;
            theSpecies->theGeometries->theCoordinates.Add(Dummy);
          }
        }
      }
    }
  }
}


void CheckCoverage()
{
    int *NumA = new int[NumberOfSpecies];
    int Total = 0;

    for(int j = 0; j < NumberOfSpecies; ++j)
        NumA[j] = 0;
    theGrid->GetSurfaceAtoms(NumA, Total);
    for(int i = 0; i < NumberOfSpecies; ++i)
    {
        if(NumA[i] != Listing[i]->NumberOfAtoms)
        {
            cout << "Garbage" << endl;
        }
    }
    delete[]NumA;
    return;
}

int Simulation::FindChargeIndex(double SC)
{
  int got = 0;
  int ind = 0;
  while (ind < NumberChargeType && got == 0)
  {
    if (Listing[0]->BindingEnergy[5][0][ind] == SC)
    {
      got = 1;
    }
    else
      ++ind;
  }
  if (got == 0)
  {
    cout << "No such charge exists" << endl;
    exit(1);
  }
  return ind;
}

void Simulation::AnalyzeReactions(double InitialTime)
{

    double DeltaT = Time - InitialTime;

    for(int i = 0; i < NumberOfSimReactions; ++i)
    {
        ElementaryRxn *PickedReaction = SimReactions[i];

        switch (PickedReaction->GetProcess())
        {
            case Ad1:
            case Ad2:
            case EleyAd:
            case EleyAdEleyDes:
            case ReactFromGas:
                PickedReaction->MacroscopicBarrier =
                    -RT * log(PickedReaction->MacroscopicBarrier /
                    (PickedReaction->NumberOfScenarios *
                        PickedReaction->GetPreexponential() *
                        PickedReaction->Reactant[0]->
                        CalculateAdsorptionRate(GetTemperature())));
                break;
            default:
                PickedReaction->MacroscopicBarrier =
                    -RT * log(PickedReaction->MacroscopicBarrier /
                    (PickedReaction->NumberOfScenarios *
                        PickedReaction->GetPreexponential()));
                break;
        }
        PickedReaction->AverageSimulationRate =
            PickedReaction->AverageSimulationRate / DeltaT;
    }
    WriteReactions();
    WriteBarrier();
    return;
}

double alpha = 0.;
double beta = 0.;
double SpeedUp;
    
void Simulation::Simulate(fstream &fout)
{
    //Startup populates the grid with initial coverages, as provided in the input file and execute diffusion
    Startup(fout);
    // Write the coordinates of adsorbates in .xyz, minimal and Hansen format
    //theGrid->WriteCoordinates(0);
    theGrid->WriteCoordinatesXYZ(0);
    //theGrid->WriteCoordinatesHan();
    // Initialize the alpha and beta variables for TPD simulations as 0
    alpha = 0.0;
    beta = 0.0;
    
    // The maximum time step is the input delta t given in the input file
    TimeStep = running_time_step;
    
    // Creating the header of the Equilibrium ratio file
//    const char *d1 = "Equilibrium_Ratio";
//    char Equilibrium_File[100];
//    strcpy(Equilibrium_File, d1);
//    strcat(Equilibrium_File, rnk);
//    ofstream covout(Equilibrium_File, ios::app);
//    covout << "Time";
//    int rxn_cntr = 0; 
//    ElementaryRxn *thisReaction;
//    Rxns->ResetToFront();
//    while (*Rxns)
//    {
//      thisReaction = Rxns->GetPtr();
//      if (thisReaction->isReverse == false && thisReaction->GetProcess() != Diff)
//      {
//        covout << setw(10) << thisReaction->GetName();
//      }
//      rxn_cntr++;
//      ++(*Rxns);
//    }
//    covout << endl;
//    covout.close();
    
    // Removing the list with last 200 events
    for (int i = 0; i < NumberOfSimReactions; i++)
    {
      dummy[i].RemoveAll();
    }
    
    // Initialization of variables
//    double Forw_rate[rxn_cntr];
//    double Back_rate[rxn_cntr];
    int TimePeriod_rat;
    //Univ = 1;
    // Setting up the KMC loop
    time(&OldTime);
    time(&TrackTime);
    NetTime = TrackTime - OldTime;
    while(Time < StopTime)
    {
        ++SimulationStep;
        ++CycleStep;
        if (CycleStep > CycleTime)
          CycleStep = 0;
        ifstream dein("StopFile.dat");
        int StopParm = 0;
        dein >> StopParm >> SpeedUp;
        // Making sure speed up makes sense
        if(SpeedUp <= 0.0)
        {
            SpeedUp = 1.0;       // 0 is illegal
        }
        // If the stop parameter in stopfile says to stop, stop
        if(StopParm == 1)
        {
            exit(1);
        }
        dein.close();
        // Executing the KMC step
        //FindTime();
        //cout << "Before EquilSteps: " << NetTime << endl;
        EquilSteps(20);
        //FindTime();
        //cout << "After EquilSteps: " << NetTime << endl;
        //cout << SimulationStep << endl;
        VTS_KineticMC();
        //FindTime();
        //cout << "After VTS_KineticMC: " << NetTime << endl;
        //cout << "Done VTS_Kinetic" << endl;
//        WriteCoverage();
//        WriteDesorbed();
//        WriteBindingEnergies();
//        theGrid->WriteCoordinates(SimulationStep);
//        theGrid->WriteCoordinatesXYZ(SimulationStep);
        //theGrid->WriteCoordinatesHan();
        //SimulationSummary();
//--------------------------------------------------------------------------------
        //Finding the equilibrium ratio
//        int rxn_cntr = 0;
//        Rxns->ResetToFront();
//        ofstream covout(Equilibrium_File, ios::app);
//        covout << Time;
//        ElementaryRxn *thisReaction;
//        while (*Rxns)
//        {
//          thisReaction = Rxns->GetPtr();
//          if (thisReaction->isReverse == false && thisReaction->GetProcess() != Diff)
//          {
//            ElementaryRxn *Forw = thisReaction;
//            ElementaryRxn *Back = thisReaction->GetReverseReaction();
//            Forw_rate[rxn_cntr] += Forw->Rate*VTS_TimeStep;
//            Back_rate[rxn_cntr] += Back->Rate*VTS_TimeStep;
//            double rat = 0;
//            if (Forw_rate[rxn_cntr] != 0)
//            {
//              if (Forw_rate[rxn_cntr] - Back_rate[rxn_cntr] >= 0)
//                rat = (Forw_rate[rxn_cntr] - Back_rate[rxn_cntr])/Forw_rate[rxn_cntr];
//              else
//                rat = (Forw_rate[rxn_cntr] - Back_rate[rxn_cntr])/Back_rate[rxn_cntr];
//              
//            }
//            covout << setw(20) << rat;
//          }
//          rxn_cntr++;
//          ++(*Rxns);
//        }
//        covout << endl;
//        covout.close();
//-----------------------------------------------------------------------------------
//      Testing the update in finding charge index
//        Time = Time + 0.05;
//------------------------------------------------------------------------------------
        //cout << "Old Charge Index: " << oldChargeIndex << " Charge Index is: " << ChargeIndex << endl;
        oldChargeIndex = ChargeIndex;
        TimePeriod_rat = Time/PeriodTime;
        double rem_time = Time - TimePeriod_rat*PeriodTime;
        double cum_time = 0.;
        for (int i = 0; i < NumberChargeType; ++i)
        {
          cum_time = cum_time + TimePeriod[i];
          if (cum_time > rem_time)
          {
            ChargeIndex = i;
            break;
          }
        }
        SurfaceCharge = ChargePattern[ChargeIndex];
        //cout << "Time: " << Time << " Charge is: " << SurfaceCharge << endl;
        if (oldChargeIndex != ChargeIndex)
        {
          //cout << "Updating coordinates at charge: " << ChargePattern[oldChargeIndex] << " to: " << ChargePattern[ChargeIndex] << endl;
          CycleStep = 0;
          for (int i = 0; i < NumberOfSimReactions; i++)
          {
            dummy[i].RemoveAll();
          }
          unThrottleAllReactions();
          Scale();
          Define_SimRxns();
          //cout << "Done updating coordinates at charge: " << ChargePattern[oldChargeIndex] << " to: " << ChargePattern[ChargeIndex] << endl;
          //theGrid->WriteCoordinatesXYZ(++SimulationStep);
          theModel->CheckInteractions();
          //theGrid->WriteCoordinatesXYZ(++SimulationStep);
        }
        //cout << SurfaceCharge << endl;
    }
    Time = 0.0;
    return;
}

void Simulation::FindTime()
{
  time(&TrackTime);
  NetTime = TrackTime - OldTime;
  OldTime = TrackTime;
  return;
}

void Simulation::CalculateVTSReactionRates()
{
    ElementaryRxn *PickedReaction;
    ScenarioAction theBasics;
    for(int ii = 0; ii < NumberOfSimReactions; ++ii)
    {
        PickedReaction = SimReactions[ii];
        if(theReactions->ReactionPossible(PickedReaction) != NULL)
        {
            PickedReaction->Possible = true; // DWDWDW
        }
        else
        {
            PickedReaction->Possible = false; // DWDWDW
        }
    }
    ElementaryRxn *RxnToConsider;
    int Layer = theGrid->AdsorbateLayer;        // z dimension
    GridSite *AnalSite;

    // Reactant and Product Sites and Energies
    GridSite *ProdSiteA[SITE_NUM_SCENARIOS];
    double ProdSiteA_E[SITE_NUM_SCENARIOS];
    GridSite *ProdSiteB[SITE_NUM_SCENARIOS];
    double ProdSiteB_E[SITE_NUM_SCENARIOS];
    GridSite *ReactSiteA[SITE_NUM_SCENARIOS];
    double ReactSiteA_E[SITE_NUM_SCENARIOS];
    GridSite *ReactSiteB[SITE_NUM_SCENARIOS];
    double ReactSiteB_E[SITE_NUM_SCENARIOS];
    double Barrier;
    TotalReactionRate = 0.0;
//    for(int ri = 0; ri < NumberOfSimReactions; ++ri)
//    {
//      cout << SimReactions[ri]->GetName() << " " << SimReactions[ri]->GetProcessName() << endl;
//    }
//    exit(1);
    int NumbReact[NumberOfSimReactions];
    for (int i = 0; i<NumberOfSimReactions; i++)
      NumbReact[i] = 0;
    // x dimension
    //cout << "Hi ";
    for(int i = 0; i < theGrid->EDim[0]; ++i)
    {
        // y dimension
        for(int j = 0; j < theGrid->EDim[1]; ++j)
        {
            AnalSite = theGrid->Surface[Layer][i][j];   // Chosen Site
            if(AnalSite != NULL && AnalSite->SurroundSize > 0)
            {
                for(int ri = 0; ri < NumberOfSimReactions; ++ri)
                {
                    RxnToConsider = SimReactions[ri];
                    thisScenario->SetReaction(RxnToConsider);
                    
//                    if (SimReactions[10] == RxnToConsider)
//                    {
//                      int a = 3;
//                      cout << RxnToConsider->GetName() << " " << thisScenario->Screen(AnalSite) << " " << Listing[a]->GetName() << " " << AnalSite->GetType()->GetName() << endl;
//                    }
                    if(RxnToConsider->Possible && thisScenario->Screen(AnalSite))
                    {    
                        int Cntr = 0;   // Number Of Reaction Scenarios Found
                        // Find the Best Product Sites and Energies
                        //cout << "Before Rxn Scenarios " << RxnToConsider->GetName() << endl;
                        Cntr = theScenarioBuilder->
                            CalculateRxnScenarios(AnalSite, RxnToConsider,
                                ReactSiteA, ReactSiteA_E, ProdSiteA,
                                ProdSiteA_E, ReactSiteB, ReactSiteB_E,
                                ProdSiteB, ProdSiteB_E);
                        NumbReact[ri] = NumbReact[ri] + Cntr;
                        //cout << "Done Rxn Scenarios" << endl;
                        //if (SurfaceCharge == -0.17)
                        //  cout << RxnToConsider->GetName() << " " << Cntr << endl;
                        // Calculate Rate for each Possible Reaction found
                        for(int ni = 0; ni < Cntr; ++ni)
                        {
                            double RxnRate =
                                theScenarioBuilder->
                                CalculateRate(RxnToConsider,
                                    ReactSiteA_E[ni], ReactSiteB_E[ni],
                                    ProdSiteA_E[ni], ProdSiteB_E[ni], Barrier);
//                            if (SimReactions[10] == RxnToConsider)
//                              cout << SimulationStep << " " << RxnRate << endl;
                            //cout << ni << " " << RxnRate << endl;
                            if(RxnRate > 0.0)
                            {
                                thisScenario->ReactantSites[0] =
                                    ReactSiteA[ni];
                                thisScenario->ReactantSites[1] =
                                    ReactSiteB[ni];
                                thisScenario->ProductSites[0] =
                                    ProdSiteA[ni];
                                thisScenario->ProductSites[1] =
                                    ProdSiteB[ni];
                                thisScenario->ForwardRate = RxnRate;
                                RxnToConsider->NumberOfScenarios++;
                                RxnToConsider->Rate += RxnRate;
                                TotalReactionRate += RxnRate;
                                theBasics.SetUpScenarioAction(thisScenario);
                                All[RxnToConsider->_idhook].Add(theBasics);
                            }
                        } // end: for(int ni=0 ; ni<Cntr ; ++ni)
                    }
                }
            }
        }
    }
    //cout << "Hey" << endl;
//    for(int ri = 0; ri < NumberOfSimReactions; ++ri)
//    {
//      RxnToConsider = SimReactions[ri];
//      cout << RxnToConsider->GetName() << " " << NumbReact[ri] << " " << RxnToConsider->Rate << endl; 
//    }
    for(int ri = 0; ri < NumberOfSimReactions; ++ri)
    {
      RxnToConsider = SimReactions[ri];
      if (RxnToConsider->NumberOfScenarios > 0)
      {
        RxnToConsider->AvActivationEnergy = RxnToConsider->AvActivationEnergy/RxnToConsider->NumberOfScenarios;
      }
    }
    return;
}

void Simulation::CalculateEquilReactionRates()
{
    ElementaryRxn *PickedReaction;
    ScenarioAction theBasics;
    for(int ii = 0; ii < NumberOfEquilRxns; ++ii)
    {
        PickedReaction = EquilReactions[ii];
        if(theReactions->ReactionPossible(PickedReaction) != NULL)
        {
            PickedReaction->Possible = true; // DWDWDW
        }
        else
        {
            PickedReaction->Possible = false; // DWDWDW
        }
    }
    ElementaryRxn *RxnToConsider;
    int Layer = theGrid->AdsorbateLayer;        // z dimension
    GridSite *AnalSite;

    // Reactant and Product Sites and Energies
    GridSite *ProdSiteA[SITE_NUM_SCENARIOS];
    double ProdSiteA_E[SITE_NUM_SCENARIOS];
    GridSite *ProdSiteB[SITE_NUM_SCENARIOS];
    double ProdSiteB_E[SITE_NUM_SCENARIOS];
    GridSite *ReactSiteA[SITE_NUM_SCENARIOS];
    double ReactSiteA_E[SITE_NUM_SCENARIOS];
    GridSite *ReactSiteB[SITE_NUM_SCENARIOS];
    double ReactSiteB_E[SITE_NUM_SCENARIOS];
    double Barrier;
    
    TotalEquilReactionRate = 0.0;
    int NumbReact[NumberOfEquilRxns];
    for (int i = 0; i<NumberOfEquilRxns; i++)
      NumbReact[i] = 0;
    
    // x dimension
    for(int i = 0; i < theGrid->EDim[0]; ++i)
    {
        // y dimension
        for(int j = 0; j < theGrid->EDim[1]; ++j)
        {
            AnalSite = theGrid->Surface[Layer][i][j];   // Chosen Site
            if(AnalSite != NULL && AnalSite->SurroundSize > 0)
            {
                for(int ri = 0; ri < NumberOfEquilRxns; ++ri)
                {
                    RxnToConsider = EquilReactions[ri];
                    thisEquilScenario->SetReaction(RxnToConsider);
                    
                    if(RxnToConsider->Possible && thisEquilScenario->Screen(AnalSite))
                    {    
                        int Cntr = 0;   // Number Of Reaction Scenarios Found
                        
                        // Find the Best Product Sites and Energies
                        Cntr = theScenarioBuilder->
                            CalculateRxnScenarios(AnalSite, RxnToConsider,
                                ReactSiteA, ReactSiteA_E, ProdSiteA,
                                ProdSiteA_E, ReactSiteB, ReactSiteB_E,
                                ProdSiteB, ProdSiteB_E);
                        NumbReact[ri] = NumbReact[ri] + Cntr;
                        
                        // Calculate Rate for each Possible Reaction found
                        for(int ni = 0; ni < Cntr; ++ni)
                        {
                            double RxnRate = theScenarioBuilder->
                                CalculateRate(RxnToConsider,
                                    ReactSiteA_E[ni], ReactSiteB_E[ni],
                                    ProdSiteA_E[ni], ProdSiteB_E[ni], Barrier);
                            
                            if(RxnRate > 0.0)
                            {
                                thisEquilScenario->ReactantSites[0] =
                                    ReactSiteA[ni];
                                thisEquilScenario->ReactantSites[1] =
                                    ReactSiteB[ni];
                                thisEquilScenario->ProductSites[0] =
                                    ProdSiteA[ni];
                                thisEquilScenario->ProductSites[1] =
                                    ProdSiteB[ni];
                                thisEquilScenario->ForwardRate = RxnRate;
                                RxnToConsider->NumberOfScenarios++;
                                RxnToConsider->Rate += RxnRate;
                                TotalEquilReactionRate += RxnRate;
                                theBasics.SetUpScenarioAction(thisEquilScenario);
                                EquilAll[RxnToConsider->_idhook].Add(theBasics);
                            }
                        } // end: for(int ni=0 ; ni<Cntr ; ++ni)
                    }
                }
            }
        }
    }
    
    return;
}

void Simulation::CalculateStaticVariables()
{
    double SC = SurfaceCharge;
    int nn = 0;
    int ind = FindChargeIndex(SC);

    CalculateSumOfRates(nn);    // Static Contribution

    for(int i = 0; i < NumberOfSpecies; ++i)
    {
        thisComponent = Listing[i];
        double Ncount = thisComponent->NumberCount;

        if(Ncount > 0)
        {
            double temp = -RT * log(thisComponent->SimBindingEnergyCount / (Ncount));
            if(!isnan(temp))
            {
                thisComponent->SimBindingEnergy = temp;
            }
        }
        double Vcount = 0;
        double VEnergyCount = 0;

        for(int j = 0; j < 7; ++j)
        {
            Vcount += thisComponent->VacancyCount[j];
            VEnergyCount += thisComponent->SimVacancyEnergyCount[j];
        }
        if(Vcount > 0)
        {
            if(isnan(VEnergyCount) || VEnergyCount == HUGE_VAL)
                VEnergyCount = double_limit;    // max for double
            double temp = RT * log(VEnergyCount / Vcount);
            thisComponent->SimVacancyEnergy = temp;
        }
        else
        {
            thisComponent->SimVacancyEnergy = thisComponent->StableSiteEnergy[ind];
        }
    }
    return;
}

ElementaryRxn *Simulation::PickAReaction()
{
    double select = 0;
    ElementaryRxn *PickedReaction = NULL;
    double ARandom = FindRandomNumber() * TotalReactionRate;
    for(int j = 0; j < NumberOfSimReactions && PickedReaction == NULL; ++j)
    {
        ElementaryRxn *NextReaction = SimReactions[j];
        select += NextReaction->Rate;
        if(select >= ARandom)
        {
            PickedReaction = NextReaction;
        }
    }
    return PickedReaction;
}

ElementaryRxn *Simulation::PickEquilReaction()
{
    double select = 0;
    ElementaryRxn *PickedReaction = NULL;
    double ARandom = FindRandomNumber() * TotalEquilReactionRate;
    for(int j = 0; j < NumberOfEquilRxns && PickedReaction == NULL; ++j)
    {
        ElementaryRxn *NextReaction = EquilReactions[j];
        
        select += NextReaction->Rate;
        if(select >= ARandom)
        {
            PickedReaction = NextReaction;
        }
    }
    return PickedReaction;
}

void Simulation::VTS_KineticMC(int NumSteps)
{
    // Initializing the VTS_TimeStep
    VTS_TimeStep = running_time_step;
    int ind = FindChargeIndex(SurfaceCharge);
    // Initializing the Events file for output
    const char *d1 = "Events";
    char Events_file[100];
    strcpy(Events_file, d1);
    strcat(Events_file, rnk);
    ofstream evout(Events_file, ios::app);
    bool Event = true;
    for(int k = 0; k < NumSteps; ++k)
    {
        // Initialize the EventOccured variable
        bool EventOccured = false;
        // Initialize various matrices involved in defining the events
        Initialize_Matrices();
        // Evecute diffusion
        SuperFastDiffuseEquilibration();
        // Finds the favorable and unfavorable sites for reaction event count
        //CalculateStaticVariables();
        // Finds the reaction rate of all possible scenarios and store in All array
        //cout << "Entering Reaction Rate" << endl;
        //cout << "Before VTS rates" << endl;
        CalculateVTSReactionRates();
        //cout << "Done VTS rates" << endl;
        //cout << "Done Reaction Rate" << endl;
        // Sum of scenarios of the picked reaction to assist in randomly picking of a scenario
        double select = 0;
        // The picked reaction randomly
        ElementaryRxn *PickedReaction = NULL;
        double ARandom;
        // Pick a reaction to be executed
        PickedReaction = PickAReaction();
        //cout << PickedReaction->GetName() << endl;
        //cout << "Hi " << endl;
        if(PickedReaction != NULL)
        {
            // The id hook of the picked reaction
            int Ri = PickedReaction->_idhook;
            // Come to the front of the scenarios of the picked reaction
            All[Ri].ResetToFront();
            // Selecting a scenario randomly
            select = 0.;
            ARandom = FindRandomNumber() * PickedReaction->Rate;
            //cout << ARandom << endl;
            ScenarioAction PickedScenario;
            
            PickedScenario.ForwardRate = 0.0;
            ScenarioAction NextScenario;
            
            while(All[Ri] && PickedScenario.ForwardRate == 0.0)
            {
                NextScenario = All[Ri].Get();
                select += NextScenario.ForwardRate;
                if(select >= ARandom)
                {
                    PickedScenario = NextScenario;
                }
                ++All[Ri];
            }
            //cout << PickedScenario.Rxn->GetName() << endl;
            if(PickedScenario.ForwardRate > 0.0)
            {
                ARandom = FindRandomNumber();
                // Calculating the variable time step
                VTS_TimeStep = -log(ARandom) / (TotalReactionRate);
                Event = true;
                // Rechecking the possibility of the occurence of the event
                if(VTS_TimeStep > running_time_step)
                {
                    ARandom = FindRandomNumber();
                    VTS_Chance = 1.0 - exp(-TotalReactionRate * running_time_step);
                    if(VTS_Chance < ARandom)
                    {
                        Event = false;
                    }
                    VTS_TimeStep = running_time_step;
                }
                else
                {
                    VTS_Chance = 1.0;
                }
                if (NumberChargeType > 1)
                {
                  int CycleNum = Time/PeriodTime;
                  double rem_time = (double) (Time - (CycleNum)*PeriodTime);
                  double cum_time = 0.;
                  for (int i = 0; i <= ChargeIndex; ++i)
                  {
                    cum_time = cum_time + TimePeriod[i];
                  }
                  rem_time = cum_time - rem_time;
                  //cout << "Time: " << Time << " CycleNum: " << CycleNum << " rem time: " << rem_time << " cum_time: " << cum_time << endl;
                  //cout << CycleNum << " " << SwitchTime << " " << rem_time << " " << VTS_TimeStep << endl;
                  if (VTS_TimeStep > rem_time)
                  {
                    //cout << "KMC Fast forwarded: " << Time << " " << CycleNum*SwitchTime << " " << rem_time << " " << VTS_TimeStep << endl;
                    if (rem_time < 0.1*TimePeriod[ChargeIndex] || SearchEvent > -1)
                    {
                      Event = false;
                      if (rem_time > 1e-10)
                        VTS_TimeStep = rem_time;
                      else
                        VTS_TimeStep = 1e-10;
                      //cout << VTS_TimeStep << endl;
                      SearchEvent = 0;
                    }
                    else
                    {
                      Event = false;
                      VTS_TimeStep = 0.;
                      ++SearchEvent;
                      //cout << SearchEvent << endl;
                    }
                  }   
                }
                else
                {
                  if (VTS_TimeStep >= running_time_step && Event == false)
                  {
                    if (SearchEvent > 10)
                    {
                      Event = false;
                      SearchEvent = 0;
                    }
                    else
                    {
                      Event = false;
                      VTS_TimeStep = 0.;
                      ++SearchEvent;
                    }
                  }  
                }
                // Adding the time step in the Cycle time log
                // If event occurs, then:
                if(Event)
                {
                    // Dybeck's approach
                    // Increment the reaction counter of the picked reaction
                    if (1)
                    {
                      CycleTimeLog[CycleStep] = VTS_TimeStep;
  //                    for (int i = 0; i < NumberOfSimReactions; i++)
  //                    {
  //                      cout << SimReactions[i]->GetName() << " " << SimReactions[i]->RxnCounter << " " << SimReactions[i]->isQE << " " << SimReactions[i]->isExec << " " << SimReactions[i]->isTotExec << endl;
  //                    }
  //                    cout << endl;
                      if (dummy[PickedReaction->_idhook].Size() > 0)
                      {
                        dummy[PickedReaction->_idhook].ResetToFront();
                        dummy[PickedReaction->GetReverseReaction()->_idhook].ResetToFront();
                        PickedReaction->First += dummy[PickedReaction->_idhook].Get();
                        PickedReaction->GetReverseReaction()->First += dummy[PickedReaction->GetReverseReaction()->_idhook].Get();
                      }
                      //cout << "Hi" << endl;
                      
                      PickedReaction->RxnCounter += 1;
                      if (PickedReaction->RxnCounter > MinExecSuperbasin)
                        PickedReaction->isExec = 1;
                      else
                        PickedReaction->isExec = 0;
                      
                      if (PickedReaction->GetProcess() != Diff)
                      {
                        //cout << PickedReaction->_idhook << endl;
                        dummy[PickedReaction->_idhook].Add(1);
                        dummy[PickedReaction->GetReverseReaction()->_idhook].Add(-1);
                        //cout << PickedReaction->GetProcessName() << endl;
                        PickedReaction->Last +=1;
                        PickedReaction->GetReverseReaction()->Last -=1;
                        //cout << PickedReaction->GetProcessName() << endl;
                      }
                      else
                      {
                        dummy[PickedReaction->_idhook].Add(0);
                        PickedReaction->Last = 0;
                      }
                      //cout << "Hi" << endl;
                      if (dummy[PickedReaction->_idhook].Size() >= MaxRxnQueueSize)
                      {
                        PickedReaction->isTotExec = 1;
                        PickedReaction->GetReverseReaction()->isTotExec = 1;
                      }
                      while (dummy[PickedReaction->_idhook].Size() > MaxRxnQueueSize)
                      {
                        dummy[PickedReaction->_idhook].ResetToFront(); 
                        dummy[PickedReaction->_idhook].Remove();
                      }
                      while (dummy[PickedReaction->GetReverseReaction()->_idhook].Size() > MaxRxnQueueSize)
                      {
                        dummy[PickedReaction->GetReverseReaction()->_idhook].ResetToFront(); 
                        dummy[PickedReaction->GetReverseReaction()->_idhook].Remove();
                      }
                      for (int i = 0; i < NumberOfSimReactions; i++)
                      {
                        ElementaryRxn *thisReaction = SimReactions[i];
                        SimReactions[i]->Total_Rate += SimReactions[i]->Rate*VTS_TimeStep;
                        SimReactions[i]->CycleRateLog[CycleStep] = SimReactions[i]->Rate;
                      }
                      if (PickedReaction->isTotExec == 1)
                      {
                        if (1)// PickedReaction->Checked_dummy == 0)
                        {
                          int forw_cntr = 0;
                          int rev_cntr = 0;
                          dummy[PickedReaction->_idhook].ResetToFront();
                          while (dummy[PickedReaction->_idhook])
                          {
                            int m = dummy[PickedReaction->_idhook].Get();
                            if (m == 1)
                              ++forw_cntr;
                            else if (m == -1)
                              ++rev_cntr;
                            else
                            {
                              cout << "Why this value propped up? " << m << endl;
                              exit(1);
                            }
                            PickedReaction->Forw_Rev_Diff = forw_cntr - rev_cntr;
                            PickedReaction->Checked_dummy = 1;
                            ++(dummy[PickedReaction->_idhook]);
                          }
                        }
  //                      else
  //                      {
  //                        PickedReaction->Forw_Rev_Diff = PickedReaction->Forw_Rev_Diff - First + Last;
  //                      }
                        if (PickedReaction->Forw_Rev_Diff > QEMaxDiff || PickedReaction->Forw_Rev_Diff < -QEMaxDiff)
                        {
                          PickedReaction->isQE = 0;
                          PickedReaction->GetReverseReaction()->isQE = 0;
                        }
                        else
                        {
                          PickedReaction->isQE = 1;
                          PickedReaction->GetReverseReaction()->isQE = 1;
                        }
                      }
                      else
                      {
                        PickedReaction->isQE = 0;
                        PickedReaction->GetReverseReaction()->isQE = 0;
                      }
                      
  //                        cout << thisReaction->GetName() << " " << diff << " " << " " << dummy[thisReaction->_idhook].Size() << " " << 
  //                        thisReaction->isTotExec << " " << thisReaction->isQE << " " << thisReaction->Scale << " " << isnonQEexecuted << endl;
                        
                        //cout << endl;
                      // Checking whether a non quasi-equilibriated reaction has been executed to calculate the scaling factor in the new cycle
                      if (CycleStep == CycleTime)
                      {
                        Tao = 0;
                        if (isnonQEexecuted)
                        {
                          isnonQEexecuted = false;
                        }
                        else
                        {
                          // Calculating the time spent in basin
                          for (int i = 0; i <= CycleTime; i++)
                            Tao += CycleTimeLog[i];
                          // Find the rate of reaction weighted by time spent in a configuration
                          for (int i = 0; i < NumberOfSimReactions; i++)
                          {
                            ElementaryRxn *thisReaction = SimReactions[i];
                            thisReaction->CycleAveragedRate = 0;
                            for (int j = 0; j <= CycleTime; j++)
                            {
                              thisReaction->CycleAveragedRate += 1/Tao*CycleTimeLog[j]*thisReaction->CycleRateLog[j];
                            }                         
                          }
                          //cout << endl;
                          // Finding the scaling factor for each reaction
                          double temp_k_throttle_point = 0.0;
                          for (int i = 0; i < NumberOfSimReactions; i++)
                          {
                            ElementaryRxn *thisReaction = SimReactions[i];
                            if (thisReaction->isQE == 0)
                            {
                              if (thisReaction->CycleAveragedRate > temp_k_throttle_point)
                              {
                                temp_k_throttle_point = thisReaction->CycleAveragedRate;
                                throttleRxn = thisReaction;
                              }
                            }
                            else if (thisReaction->isExec == 0)
                            {
                              double tempk = thisReaction->CycleAveragedRate;
                              tempk = (tempk + thisReaction->GetReverseReaction()->CycleAveragedRate)/2;
                              if (tempk > k_buff*temp_k_throttle_point)
                              {
                                throttleRxn = thisReaction;
                                temp_k_throttle_point = tempk;
                              }
                            }
                          }
                          if (temp_k_throttle_point > 1e-309)
                          {
                            for (int i = 0; i < NumberOfSimReactions; i++)
                            {
                              ElementaryRxn *thisReaction = SimReactions[i];
                              if (thisReaction->isQE == 1 && thisReaction->isExec == 1 && thisReaction->isReverse == 0)
                              {
                                double dummy = thisReaction->Scale;
                                dummy = 2*temp_k_throttle_point*dummy/(thisReaction->CycleAveragedRate + thisReaction->GetReverseReaction()->CycleAveragedRate);
                                if (dummy > 1)
                                {
                                  thisReaction->Scale = 1;
                                  thisReaction->GetReverseReaction()->Scale = 1;
                                }
                                else
                                {
                                  thisReaction->Scale = dummy;
                                  thisReaction->GetReverseReaction()->Scale = dummy;
                                }
                              }
                            }
                          }
                        }
                      }
                      ElementaryRxn *thisReaction;
                      //Resetting rate constants if reaction executed is slow
                      if (lastNonQEReaction != NULL && PickedReaction != lastNonQEReaction && PickedReaction != lastNonQEReaction->GetReverseReaction() && PickedReaction->isQE == 0)
                      {
                        // cout << "NonQE reaction: " << PickedReaction->GetName() << endl;
                        // No need to scale the reactions in this cycle
                        nonQEEventThisCycle = 1;
                        // Scale, RxnCounter, isExec and Total_Rate for every reaction is reset
                        unThrottleAllReactions();
                      }
                      if (PickedReaction->isQE == 0)
                        lastNonQEReaction = PickedReaction;
                      
                      Scale();
                    }
//                    cout << SimulationStep << " " << Time << " " << VTS_TimeStep << " " << PickedReaction->GetName() << " RxnCntr: " << PickedReaction->RxnCounter << " Diff: " << PickedReaction->Forw_Rev_Diff << " isEquil: " << PickedReaction->isQE << " isExec: " << PickedReaction->isExec << " " << PickedReaction->Scale << " " << lastNonQEReaction->GetName() << endl;
//                    cout << "Number of CO on surface = " << Listing[2]->NumberCount << endl;
                    EventOccured = true;
                    thisScenario->SetUpScenario(PickedScenario);
                    assert(PickedScenario.Rxn == PickedReaction);
                    evout << Time << "\t\t" << SurfaceCharge << "\t\t" << thisScenario->theReaction->GetName().Ec_str() << "\t\tTimeStep: " << VTS_TimeStep << endl;
                    SearchEvent = 0;
                    thisScenario->doit();
                    theModel->CalculateBindingEnergy(PickedScenario.ProductSites[0]);
                    theModel->CalculateBindingEnergy(PickedScenario.ProductSites[1]);
                    alpha = 0.;
                    beta -= 0.05 / SpeedUp;
                    if(beta < 0.)
                        beta = 0.;
                    double EstRate = TotalReactionRate;
                    
//                    running_time_step = 1.0 / EstRate;
//                    double OldTimeStep = running_time_step;
//                    if(running_time_step > 0.2 * SpeedUp * OldTimeStep)
//                        running_time_step = 0.2 * SpeedUp * OldTimeStep;
//                    // Shorten the time step in preparation for faster events that
//                    // result from this reaction
//                    else if(running_time_step < 1.e-3 * OldTimeStep)
//                        running_time_step = 1.e-3 * OldTimeStep;        // no drastic movement
                }
//                else
//                {
//                    // Having shortened the time step
//                    running_time_step =
//                        (1.2 + SpeedUp*(alpha + 0.2*beta))*running_time_step;
//                    
//                    // we can now increase it at a faster rate
//                    // cautiously increase with beta
//                    if(beta < 0.2) beta += 0.01;
//                    else if(beta < 0.4) beta += 0.005;
//                    else if(beta < 0.6) beta += 0.0025;
//                    else if(beta < 0.8) beta += 0.00125;
//                    else if(beta > 1.0) beta = 1.0;
//                    else beta += 0.0001;
//                    if(alpha > 0.2)
//                    {
//                        alpha = 0.;
//                    }
//                    else
//                    {
//                        alpha += 0.02;
//                    }
//                }
//                if(running_time_step > input_time_step)
//                {
//                    running_time_step =
//                        exp((1. - beta) * log(input_time_step) +
//                            beta * log(running_time_step));
//                    // cautiously increase with beta
//                }
//                else if(running_time_step < 1.e-5 * input_time_step)
//                {               // Not truly interested in this time regime
//                    running_time_step =
//                        exp(0.5 * log(input_time_step) +
//                            0.5 * log(running_time_step));
//                }
            }
            else
            {
                int CycleNum = Time/PeriodTime;
                double rem_time = (double) (Time - (CycleNum)*PeriodTime);
                double cum_time = 0.;
                for (int i = 0; i <= ChargeIndex; ++i)
                {
                  cum_time = cum_time + TimePeriod[i];
                }
                rem_time = cum_time - rem_time;
                if (rem_time > running_time_step)
                  VTS_TimeStep = running_time_step;
                else
                  VTS_TimeStep = rem_time;
                //cout << "If forward rate = 0: " << VTS_TimeStep << endl;
            }
            for(int i = 0; i < NumberOfSimReactions; ++i)
            {
                SimReactions[i]->AverageSimulationRate +=
                    (VTS_TimeStep * SimReactions[i]->Rate);
                SimReactions[i]->MacroscopicBarrier +=
                    SimReactions[i]->Rate;
            }
        }
        //cout << "Hey" << endl;
        Time += VTS_TimeStep;
        //cout << "Final VTS_TimeStep: " << VTS_TimeStep << " Time: " << Time << endl;
        AnalyzeReactions(Time - VTS_TimeStep);
        Temperature += Beta * VTS_TimeStep;
        RT = R_GasConst * Temperature;
        if (Event)
        {
          WriteCoverage();
          WriteDesorbed();
          WriteBindingEnergies();
          theGrid->WriteCoordinatesXYZ(SimulationStep);
        }
    }
    evout.close();
}

void Simulation::EquilSteps(int NumSteps)
{
    if (NumberOfEquilRxns == 0)
    {
      return;
    }
    
    for(int k = 0; k < NumSteps; ++k)
    {
        // Initialize the EventOccured variable
        // Initialize various matrices involved in defining the events
//        FindTime();
//        cout << "At start: " << NetTime << endl;
        Initialize_Equil_Matrices();
        // Evecute diffusion
        SuperFastDiffuseEquilibration();
        // Finds the favorable and unfavorable sites for reaction event count
        int nn = 0;
        //CalculateSumOfRates(nn);
        // Finds the reaction rate of all possible scenarios and store in All array
//        FindTime();
//        cout << "After sum of rates: " << NetTime << endl;
        CalculateEquilReactionRates();
//        FindTime();
//        cout << "After Calculation of rates: " << NetTime << endl;
        // Sum of scenarios of the picked reaction to assist in randomly picking of a scenario
        double select = 0;
        // The picked reaction randomly
        ElementaryRxn *PickedReaction = NULL;
        double ARandom;
        // Pick a reaction to be executed
        PickedReaction = PickEquilReaction();
        if(PickedReaction != NULL)
        {
            // The id hook of the picked reaction
            int Ri = PickedReaction->_idhook;
            // Come to the front of the scenarios of the picked reaction
            EquilAll[Ri].ResetToFront();
            // Selecting a scenario randomly
            select = 0.;
            ARandom = FindRandomNumber() * PickedReaction->Rate;
            ScenarioAction PickedScenario;
            
            PickedScenario.ForwardRate = 0.0;
            ScenarioAction NextScenario;
            
            while(EquilAll[Ri] && PickedScenario.ForwardRate == 0.0)
            {
                NextScenario = EquilAll[Ri].Get();
                select += NextScenario.ForwardRate;
                if(select >= ARandom)
                {
                    PickedScenario = NextScenario;
                }
                ++EquilAll[Ri];
            }
            if(PickedScenario.ForwardRate > 0.0)
            {
                bool Event = true;
                thisEquilScenario->SetUpScenario(PickedScenario);
                assert(PickedScenario.Rxn == PickedReaction);
                thisEquilScenario->doit();
                theModel->CalculateBindingEnergy(PickedScenario.ProductSites[0]);
                theModel->CalculateBindingEnergy(PickedScenario.ProductSites[1]);
            }
        }
    }
}


void Simulation::unThrottleAllReactions()
{
  for (int i = 0; i < NumberOfSimReactions; i++)
  {
    ElementaryRxn *thisReaction = SimReactions[i];
    thisReaction->Scale = 1;
    thisReaction->RxnCounter = 0;
    thisReaction->Total_Rate = 0;
    thisReaction->isExec = 0;
  }
}

void Simulation::Scale()
{
  for (int i = 0; i < NumberOfSimReactions; i++)
  {
    ElementaryRxn *thisReaction = SimReactions[i];
    if (thisReaction->isQE)
    {
      thisReaction->SetUpdPreexponential(thisReaction->Scale*thisReaction->GetPreexponential());
    }
  }
}

void Simulation::Startup(fstream &fout)
{
    GasOff = false; // DWDWDW
    thisComponent = theNULLSpecies;
    thisSite = NULL;
    if(Population)
        theGrid->DePopulate();
    int TotalCoverage = 0;
    int cov;
    int *Coverage = new int[NumberOfSpecies];
    
    //---------------------------------------------------------
    //Experimenting with adsorbate orientations
//    Component* Spec = Listing[4];
//    cout << Spec->GetName() << endl;
//    Spec->NumberOfAtoms = 1;
//    GridSite* Site = theGrid->LocateCentralSite(1);
//    Site->SetType(Spec);
//    double a = theModel->CalculateBindingEnergy(Site, 1);
//    cout << a << " " << Listing[3]->GetClosedShell() << endl;
//    for (int i = 0; i < 6; ++i)
//    {
//      Site->Orientation = i;
//      theModel->CalculateBindingEnergy(Site, 1);
//      theGrid->WriteCoordinatesXYZ(i);
//    }
//    exit(1);
    //---------------------------------------------------------
    
    for(int i = 0; i < NumberOfSpecies; ++i)
    {
        thisComponent = Listing[i];
        cov = thisComponent->NumberOfAtoms;
        TotalCoverage += cov;
        Coverage[i] = cov;
       // cout << thisComponent->GetName() << " " << i << endl;
    }
    if(Population)
    {
        Populate(Coverage, fout);
        for (int j = 0; j < 10; j++)
        {
          SuperFastDiffuseEquilibration(0.25);
        }
    }
}

int Simulation::UpdateEnvironment()
{
    // If there is an applied partial pressure and it is to be removed
    if(!GasOff && RemoveGas)
    {
        // If the time is greater than the delay time
        if(Time > DelayTime)
        {
            theSpecies->theSpecies.ResetToFront();
            cout << "\n\n\n"
                << "As Requested, Removing ambient pressure of all species"
                << endl;
            
            while(theSpecies->theSpecies)
            {
                Component *newSpecies = theSpecies->theSpecies.GetPtr();
                
                newSpecies->SetPressure(0.0);
                ++theSpecies->theSpecies;
            }
            GasOff = true; // DWDWDW
            running_time_step = TimeStep;
        }
    }
    return false; // DWDWDW
}


int Simulation::Initialize_Matrices()
{
    /* Zero Matrices and clear lists */
    // Set elements to zero with initialize
    // A clean up operation
    for(int i = 0; i < NumberOfSpecies + 1; ++i)
    {
        Listing[i]->SimBindingEnergyCount = 0.;
        Listing[i]->NumberCount = 0.;
        for (int ii = 0; ii < 10; ii++)
            Listing[i]->SimEnergy[ii] = 0.0;
        for(int iii = 0; iii < 7; ++iii)
        {
            Listing[i]->SimVacancyEnergyCount[iii] = 0.;
            Listing[i]->VacancyCount[iii] = 0.;
        }
    }
    for(int i2 = 0; i2 < NumberOfSimReactions; ++i2)
    {
        All[i2].RemoveAll();
        thisReaction = SimReactions[i2];
        Barrier[thisReaction->_id] = 0.0;
        Prefactor[thisReaction->_id] = 0.0;
        thisReaction->Rate = 0.0;
        thisReaction->NumberOfScenarios = 0;
        thisReaction->AverageSimulationRate = 0.0;
        thisReaction->MacroscopicBarrier = 0.0;
        thisReaction->AvActivationEnergy = thisReaction->UserSetActivationEnergy;
    }
    
    return 0;
}

int Simulation::Initialize_Equil_Matrices()
{
    for(int i2 = 0; i2 < NumberOfEquilRxns; ++i2)
    {
        EquilAll[i2].RemoveAll();
        thisReaction = EquilReactions[i2];
        Barrier[thisReaction->_id] = 0.0;
        Prefactor[thisReaction->_id] = 0.0;
        thisReaction->Rate = 0.0;
        thisReaction->NumberOfScenarios = 0;
        thisReaction->AverageSimulationRate = 0.0;
        thisReaction->MacroscopicBarrier = 0.0;
    }
    return 0;
}

void Simulation::CalculateSumOfRates(int &NumberOfTimes)
{
    double SC = SurfaceCharge;
    /* Calculate and return the favorable sites for surface reaction */
    int Layer = theGrid->AdsorbateLayer;        // z dimension
    double SiteBinding = 0.;
    NumberOfTimes++;
    int index;
    GridSite *AnalSite;

    RateTotal = 0.;
    
    // DW 8/24/04: TestEnergy was unused, so commented it out
    //double TestEnergy;

    double BestEnergy;

    // DW 8/24/04: TestC was unused, so commented it out
    //Component *TestC;

    TotalSystemEnergy = 0.;
    GridSite *UnFavorableSite = NULL;
    double UnFavorableEnergy = 0.1;
    GridSite *FavorableSite[MAX_NUM_SPECIES];
    double FavorableEnergy[MAX_NUM_SPECIES];

    for(index = 0; index < NumberOfSpecies; ++index)
    {
        FavorableEnergy[index] = 0.;
        FavorableSite[index] = NULL;
    }
    int BindingEnergyTooLow = false; // DWDWDW

    for(int i = 0; i < theGrid->EDim[0]; ++i)
    {                           // x dimension
        for(int j = 0; j < theGrid->EDim[1]; ++j)
        {                       // y dimension
            AnalSite = theGrid->Surface[Layer][i][j];   // Chosen Site
            if(AnalSite != NULL && AnalSite->SurroundSize > 0)
            {                   // If Exists   
                Component *siteComponent = AnalSite->GetType();

                index = siteComponent->CrossReference();
                if(siteComponent != theNULLSpecies)
                {
                    SiteBinding = theModel->CalculateBindingEnergy(AnalSite, true); //DWDWDW
                    
                    // the equilibrated version is unstable for unknown reasons
                    //AnalSite->BindingEnergy = SiteBinding;
                    TotalSystemEnergy += SiteBinding;
                    if(SiteBinding < 0.)
                    {
                        BestEnergy = 0.;
                        for(int ni = 0; ni < AnalSite->SurroundSize; ++ni)
                        {
                            if(AnalSite->Surround[ni]->Type->_id)
                            {
                                //AnalSite->Surround[ni]->BindingEnergy = theModel->CalculateBindingEnergy(AnalSite->Surround[ni]);
                            }   // Try a different orientation for the neighbors
                        }
                        TotalSystemEnergy -= SiteBinding;
                        SiteBinding = theModel->CalculateBindingEnergy(AnalSite);
                        if(SiteBinding > 0.)
                        {
                            //AnalSite->BindingEnergy = SiteBinding;
                            TotalSystemEnergy += SiteBinding;
                        }
                        else
                        {
                            BindingEnergyTooLow = true; // DWDWDW
                            if(SiteBinding < UnFavorableEnergy)
                            {
                                UnFavorableEnergy = SiteBinding;
                                UnFavorableSite = AnalSite;
                            }
                        }
                    }
                    siteComponent->SimBindingEnergyCount += (double)exp(-SiteBinding / RT);
                    siteComponent->NumberCount += 1.;
                }
                else
                {
                    for(int ii = NumberOfSpecies; ii--;)
                    {
                        AnalSite->Type = Listing[ii];
                        Eshort NP[3];
                        AnalSite->GetGridPos(NP);
                        
                        
                        SiteBinding = theModel->CalculateBindingEnergy(AnalSite);
                        Univ = 0;
                        if(SiteBinding > 0.)
                        {
                            siteComponent = Listing[ii];
                            if(SiteBinding > FavorableEnergy[ii])
                            {
                                FavorableSite[ii] = AnalSite;
                            }
                            int Coord = AnalSite->BondSize - 1;

                            siteComponent->SimVacancyEnergyCount[Coord] += (double)exp(SiteBinding / RT);
                            siteComponent->VacancyCount[Coord] += 1.;
                        }
                        AnalSite->Type = theNULLSpecies;
                    }
                }
            }                   // If Exists
        }                       // y dimension
    }                           // x dimension
    if(NumberOfTimes < 10 && BindingEnergyTooLow)
    {
        Component *Type = UnFavorableSite->Type;

        if(FavorableSite[Type->_idhook] != NULL)
        {
            UnFavorableSite->Type = theNULLSpecies;
            FavorableSite[Type->_idhook]->Type = Type;
            CalculateSumOfRates(NumberOfTimes);
        }
    }
    return;
}


//int AnalyzeSurroundings(GridSite * FirstSite, double &BestE,
//    double &CenterBE, double SC)
//{
//    int i;
//    int SurrPass = true; // DWDWDW
//    double TestE = 0;
//    double E = 0;
//    double CTEST = 0;
//    Eshort OrigOrientation[PatternDim];
//
//    for(i = 0; i < FirstSite->SurroundSize; ++i)
//    {
//        OrigOrientation[i] = FirstSite->Surround[i]->Orientation;
//    }
//    for(i = 0; i < FirstSite->SurroundSize && SurrPass; ++i)
//    {
//        if(FirstSite->Surround[i]->Type->_id)
//        {
//            E += theModel->CalculateBindingEnergy(FirstSite->Surround[i]);
//            if(FirstSite->Surround[i] == FirstSite)
//                CTEST = E;
//            TestE += E;
//            if(TestE < 0.)
//            {
//                SurrPass = false; // DWDWDW
//            }
//        }
//    }
//    if(!SurrPass || TestE < BestE)
//    {
//        for(i = 0; i < FirstSite->SurroundSize; ++i)
//        {
//            FirstSite->Surround[i]->Orientation = OrigOrientation[i];
//        }
//    }
//    if(TestE > BestE && SurrPass)
//    {
//        BestE = TestE;
//        CenterBE = CTEST;
//    }
//    return SurrPass;
//}

//// DW: Looks like this is never used.
//bool RelaxProductState(GridSite * Products[]) // DWDWDW
//{
//
//    GridSite *FirstSite = Products[0];
//    GridSite *SecondSite = Products[1];
//    Component *FirstSiteType = NULL;
//    Component *SecondSiteType = NULL;
//    double MaxBE = 0;
//    double CenterABE = 0;
//    double CenterBBE = 0;
//    GridSite *NewFirstSite = NULL;
//    GridSite *NewSecondSite = NULL;
//    int i;
//    int SurrPass = true; // DWDWDW
//    int SurrPass2 = true; // DWDWDW
//
//    if(FirstSite != NULL && FirstSite->Type != theNULLSpecies)
//    {
//        FirstSiteType = FirstSite->Type;
//        SurrPass = AnalyzeSurroundings(FirstSite, MaxBE, CenterABE);
//        FirstSite->SetType(theNULLSpecies);
//        NewFirstSite = NULL;
//        for(i = 0; i < FirstSite->SurroundSize; ++i)
//        {
//            GridSite *Surround = FirstSite->Surround[i];
//
//            if(Surround->Type == theNULLSpecies)
//            {
//                Surround->SetType(FirstSiteType);
//                double OldMaxBE = MaxBE;
//                int Possible =
//                    AnalyzeSurroundings(Surround, MaxBE, CenterABE);
//                if(Possible && MaxBE > OldMaxBE)
//                {
//                    NewFirstSite = Surround;
//                    if(!SurrPass)
//                        SurrPass = true; // DWDWDW
//                }
//                Surround->SetType(theNULLSpecies);
//            }
//        }
//    }
//    if(SecondSite != NULL && SecondSite->Type != theNULLSpecies)
//    {
//        SecondSiteType = SecondSite->Type;
//        SurrPass2 = AnalyzeSurroundings(SecondSite, MaxBE, CenterBBE);
//        SecondSite->SetType(theNULLSpecies);
//        NewSecondSite = NULL;
//        for(i = 0; i < SecondSite->SurroundSize; ++i)
//        {
//            GridSite *Surround = SecondSite->Surround[i];
//
//            if(Surround->Type == theNULLSpecies)
//            {
//                Surround->SetType(SecondSiteType);
//                double OldMaxBE = MaxBE;
//                int Possible =
//                    AnalyzeSurroundings(Surround, MaxBE, CenterBBE);
//                if(Possible && MaxBE > OldMaxBE)
//                {
//                    NewSecondSite = Surround;
//                    if(!SurrPass2)
//                        SurrPass2 = true; // DWDWDW
//                }
//                Surround->SetType(theNULLSpecies);
//            }
//        }
//    }
//    if(SurrPass && SurrPass2)
//    {
//        if(NewFirstSite)
//        {
//            NewFirstSite->SetType(FirstSiteType);
//            NewFirstSite->BindingEnergy = CenterABE;
//        }
//        else if(FirstSite)
//        {
//            FirstSite->SetType(FirstSiteType);
//            FirstSite->BindingEnergy = CenterABE;
//        }
//        if(NewSecondSite)
//        {
//            NewSecondSite->SetType(SecondSiteType);
//            NewSecondSite->BindingEnergy = CenterBBE;
//        }
//        else if(SecondSite)
//        {
//            SecondSite->SetType(SecondSiteType);
//            SecondSite->BindingEnergy = CenterBBE;
//        }
//    }
//    else
//    {
//        SurrPass = false; // DWDWDW
//    }
//    return (bool)SurrPass;
//}
