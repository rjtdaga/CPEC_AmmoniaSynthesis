#include <string>
#include <cassert>
#include <cmath>
using namespace std;


#include "ScenarioBuilder.h"
#include "Scenario.h"
#include "Simulation.h"
#include "Reactions.h"
#include "Model.h"
#include "Grid.h"
#include "Simulation.h"
#include "ShortRoutines.h"
#include "Base/Constants.h"
#include "SystemHeaders.h"


ScenarioBuilder *theScenarioBuilder = NULL;

ScenarioBuilder::ScenarioBuilder(fstream &fin, fstream &fout) :
    Facilitator(fin, fout)
{
    theScenarioBuilder = this;
    if(theGrid->GetSurface() == "100")
    {
        QuadrantIndex = 4;
    }
    else
    {
        QuadrantIndex = 6;
    }
}

void ScenarioBuilder::CalculatePairwisePreferredDistance()
{
    double SC = theSimulation->SurfaceCharge;
    GridSite *AnalSite;
    GridSite *Neighbor;
    int MyCoord;
    int NeighCoord;
    double IdealEnergy = 0;
    double NetEnergy = 0;
    double MyEnergy;
    double NeighborEnergy;
    int ind = theSimulation->FindChargeIndex(SC);
    int k;
    for (int l = 0; l < theSimulation->NumberChargeType; ++l)
    {
        for(k = 1; k < 4; ++k)
        {
            if ((k == 1 && theGrid->Atop_a == 0) || (k == 2 && theGrid->Bridge_a == 0) || (k == 3 && theGrid->Hollow_a == 0))
              continue;
              
            if(theGrid->GetSurface() == "100" && k == 3)
            {
                k = 4;
            }
            if(theGrid->GetSurface().lower() == "graphene" && k == 3)
            {
                k = 6;
            }
            AnalSite = theGrid->LocateCentralSite(k);
            for(int i = 0; i < NumberOfSpecies; ++i)
            {
                AnalSite->Type = Listing[i];
                
                MyCoord = AnalSite->BondSize - 1;
                for(int j = 0; j < NumberOfSpecies; ++j)
                {
                    double BestEnergy = 0;
                    double BestDistance = 1.e10;
    
                    for(int ni = 0; ni < AnalSite->SurroundSize; ++ni)
                    {
                        Neighbor = AnalSite->Surround[ni];
                        if (theGrid->AllSite == 1)
                        {
                          if(Neighbor != AnalSite)
                          {
                              Neighbor->Type = Listing[j];
                              
                              NeighCoord = Neighbor->BondSize - 1;
                              double Distance = AnalSite->FindDistance(Neighbor,theGrid->GridLength);
                              if(Distance < AnalSite->Type->FavorableDist[Neighbor->Type->_idhook][MyCoord][l])
                              {
                                  IdealEnergy = AnalSite->Type->StableSiteEnergy[ind] + Neighbor->Type->StableSiteEnergy[ind];
                                  MyEnergy = theModel->CalculateBindingEnergy(AnalSite);
                                  NeighborEnergy = theModel->CalculateBindingEnergy(Neighbor);
                                  NetEnergy = (MyEnergy + NeighborEnergy);
                                  if(Distance > BestDistance)
                                  {
                                      if(NetEnergy > 1.02 * BestEnergy)
                                      {
                                          BestEnergy = NetEnergy;
                                          BestDistance = Distance;
                                      }
                                  }
                                  else
                                  {
                                      if(NetEnergy > 0.98 * BestEnergy)
                                      {
                                          BestEnergy = NetEnergy;
                                          BestDistance = Distance;
                                      }
                                  }
                              }
                              Neighbor->Type = theNULLSpecies;
                          }
                        }
                    }
                    AnalSite->Type->FavorableDist[Listing[j]->_idhook][MyCoord][l] = BestDistance;
                }
                AnalSite->Type = theNULLSpecies;
            }
        }
    }
    theGrid->DePopulate();
    return;
}

bool ScenarioBuilder::FindBestSingleSite(GridSite * AnalSite,
    Component * ProductA, GridSite * &BestSiteA,
    double &BestEnergyA, double cutoff, bool Forced) // DWDWDW
{
    double SC = theSimulation->SurfaceCharge;
    int SurrSize = AnalSite->SurroundSize;
    GridSite *Neighbor;
    GridSite *Neighbor2;
    double Cmp = 0;
    BestSiteA = NULL;
    BestEnergyA = -1.;
    bool Found = 0;
    bool FoundFinal = 0;
    double NeighCmp = 0.;
//    if (theSimulation->Univ == 13)
//      cout << AnalSite->GetType()->GetName() << endl;
    int a = 0;
    int ni;
    for(ni = SurrSize; ni--;)
    {
        Found = 1;
        Neighbor = AnalSite->Surround[ni];      // Neighboring Site
//        if (Neighbor == AnalSite)
//        {
//          if (theSimulation->Univ == 13)
//          {
//            cout << "Surr Num: " << ni << endl;
//            cout << "Distance: " << AnalSite->FindDistance(Neighbor, theGrid->GridLength) << endl;
//          }
//          a = ni;
//        }
        if(Neighbor->Type == theNULLSpecies)
        {                       // two separate sites, both empty
            if(AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff)
            {
                if(ProductA)
                {
                    if(Forced)
                    {
                        Neighbor->Type = ProductA;
                        Cmp = theModel->CalculateBindingEnergy(Neighbor);
                        if (Cmp > 0.)
                        {
                          for(int ni = SurrSize; ni--;)
                          {
                            Neighbor2 = AnalSite->Surround[ni];
                            if(Neighbor->Type != theNULLSpecies)
                            {
                              NeighCmp = theModel->CalculateBindingEnergy(Neighbor2, true);
                              if (NeighCmp < 0)
                              {
                                Found = 0;
                                break;
                              }
                            }
                          }
                        }
                        else
                        {
                          Found = 0;
                        }
                        if (Found == 1)
                        {
                          if (Cmp > BestEnergyA)
                          {
                            BestEnergyA = Cmp;
                            BestSiteA = Neighbor;
                            FoundFinal = 1;
                          }
                        }
                        Neighbor->Type = theNULLSpecies;
                    }
                    else
                    {
                        Neighbor->Type = ProductA;
                        Cmp = theModel->CalculateBindingEnergy(Neighbor);
//                        if (theSimulation->Univ == 13)
//                        {
//                          cout << Cmp << " " << BestEnergyA << " " << Neighbor->BondSize << endl;
//                        }
                        if (Cmp > 0.)
                        {
                          for(int ni = SurrSize; ni--;)
                          {
                            Neighbor2 = AnalSite->Surround[ni];
                            if(Neighbor->Type != theNULLSpecies)
                            {
                              NeighCmp = theModel->CalculateBindingEnergy(Neighbor2);
//                              if (theSimulation->Univ == 13)
//                              {
//                                cout << NeighCmp << " ";
//                              }
                              if (NeighCmp < 0)
                              {
                                Found = 0;
                                break;
                              }
                            }
                          }
                        }
                        else
                        {
                          Found = 0;
                        }
//                        if (theSimulation->Univ == 13)
//                        {
//                          cout << endl;
//                        }
                        if (Found == 1)
                        {
                          if (Cmp > BestEnergyA)
                          {
                            BestEnergyA = Cmp;
                            BestSiteA = Neighbor;
                            Eshort c[3];
                            Neighbor->GetGridPos(c);
                            FoundFinal = 1;
                          }
                        }
                        Neighbor->Type = theNULLSpecies;
                    }
                }
            }
        }
    }
    return FoundFinal;
}

int ScenarioBuilder::CalculateRxnScenarios(GridSite * AnalSite,
    ElementaryRxn * RxnToConsider, GridSite * ReactSiteA[],
    double ReactSiteA_E[], GridSite * ProdSiteA[], double ProdSiteA_E[],
    GridSite * ReactSiteB[], double ReactSiteB_E[], GridSite * ProdSiteB[],
    double ProdSiteB_E[])
{
    //cout << "Hi ";
    /* Calling function must set reaction of thisScenario */
    double SC = theSimulation->SurfaceCharge;
    //cout << RxnToConsider->GetName() << endl;
    GridSite *BestSiteA[6];
    GridSite *BestSiteB[6];
    double BestEnergyA[6];
    double BestEnergyB[6];
    GridSite *BestSiteA_Backup[6];
    double BestEnergyA_Backup[6];
    int Cntr = 0;
    int qa;
    int ni;
    double cutoff0 = 3.01 * theGrid->MMdistance;
    double cutoff1 = 1.01 * theGrid->MMdistance;
    double cutoff2 = 0.99 * theGrid->MMdistance;
    int SurrSize = AnalSite->SurroundSize;
    GridSite *Neighbor = NULL;
    bool Found;
    switch (RxnToConsider->GetProcess())
    {
        case ElecTransfer:
        {
            //cout << "Entering ElecTransfer" << endl;
            Eshort tor = AnalSite->Orientation;
            ReactSiteA[0] = AnalSite;
            ReactSiteA_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteB[0] = NULL;
            ReactSiteB_E[0] = 0.;
            ProdSiteB[0] = NULL;
            ProdSiteB_E[0] = 0.;
            
            AnalSite->SetType(theNULLSpecies);
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[0],
            ProdSiteA_E[0], cutoff2, false);
            
            AnalSite->Type = RxnToConsider->Reactant[0];
            AnalSite->Orientation = tor;
            //AnalSite->BindingEnergy = theModel->CalculateBindingEnergy(AnalSite, true);
            theModel->CalculateBindingEnergy(AnalSite, true);
            if (Found)
              Cntr = 1;
            break;
        }
        case Des1:             // type A
        {
            //cout << "Entering Des1" << endl;
            ReactSiteA[0] = AnalSite;
            ReactSiteA_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteB[0] = NULL;
            ReactSiteB_E[0] = 0.;
            ProdSiteA[0] = NULL;
            ProdSiteA_E[0] = 0.;
            ProdSiteB[0] = NULL;
            ProdSiteB_E[0] = 0.;
            Cntr = 1;
            break;
        }
        case DissToGas:        // type A
        {
            Eshort tor = AnalSite->Orientation;
            ReactSiteA[0] = AnalSite;
            ReactSiteA_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteB[0] = NULL;
            ReactSiteB_E[0] = 0.;
            ProdSiteA[0] = NULL;
            ProdSiteA_E[0] = 0.;
            
            AnalSite->SetType(theNULLSpecies);
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[0],
            ProdSiteB_E[0], cutoff2, false);
            
            AnalSite->Type = RxnToConsider->Reactant[0];
            AnalSite->Orientation = tor;
            //AnalSite->BindingEnergy = theModel->CalculateBindingEnergy(AnalSite, true);
            theModel->CalculateBindingEnergy(AnalSite, true);
            if (Found)
              Cntr = 1;
            break;
        }
        case Ad1:              // type V
        {
            //cout << "Entering Ad1" << endl;
            AnalSite->SetType(RxnToConsider->Product[0]);
            double Cmp;
            double NeighCmp;
            Cmp = theModel->CalculateBindingEnergy(AnalSite);
            ProdSiteA_E[0] = Cmp;
            bool Pass;
            Pass = 1;
            if (Cmp > 0)
            {
              for (int i = 0; i < AnalSite->SurroundSize; ++i)
              {
                GridSite *Neighbor = AnalSite->Surround[i];
                if (Neighbor->GetType() != theNULLSpecies)
                {
                  NeighCmp = theModel->CalculateBindingEnergy(Neighbor, true);
                  if (NeighCmp < 0)
                  {
                    ProdSiteA_E[0] = -1;
                    Pass = 0;
                    break;
                  }
                }
              }
            }
            else
            {
              Pass = 0;
            }
            if (Pass == 1)
            {
              ProdSiteA[0] = AnalSite;
            }
            AnalSite->SetType(theNULLSpecies);
            ProdSiteB[0] = NULL;
            ProdSiteB_E[0] = 0.;
            ReactSiteB[0] = NULL;
            ReactSiteB_E[0] = 0.;
            ReactSiteA[0] = NULL;
            ReactSiteA_E[0] = 0.;
            if(ProdSiteA_E[0] > 0.)
                Cntr = 1;
            break;
        }
        case ReactFromGas:     
        {
            Eshort tor = AnalSite->Orientation;
            ReactSiteB[0] = AnalSite;
            ReactSiteB_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteA[0] = NULL;
            ReactSiteA_E[0] = 0.;
            ProdSiteB[0] = NULL;
            ProdSiteB_E[0] = 0.;
            
            AnalSite->SetType(theNULLSpecies);
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[0],
            ProdSiteA_E[0], cutoff2, false);
            
            AnalSite->Type = RxnToConsider->Reactant[1];
            AnalSite->Orientation = tor;
            //AnalSite->BindingEnergy = theModel->CalculateBindingEnergy(AnalSite, true);
            theModel->CalculateBindingEnergy(AnalSite, true);
            if (Found)
              Cntr = 1;
            break;
        }
        case Ad2:              // type VV
        {
            ReactSiteA[Cntr] = NULL;
            ReactSiteA_E[Cntr] = 0.;
            ReactSiteB[Cntr] = NULL;
            ReactSiteB_E[Cntr] = 0.;
            
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[0],
            ProdSiteB_E[0], cutoff0, false);
            if (Found == 1)
            {
              ProdSiteB[0]->Type = RxnToConsider->Product[1];
              theModel->CalculateBindingEnergy(ProdSiteB[0]);
              Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[0],
              ProdSiteA_E[0], cutoff0, false);
              ProdSiteB[0]->Type = theNULLSpecies;
            }
            
            if (Found == 1)
              Cntr = 1;
            break;
        }
        case Ad2ToGas:              // type VV
        {
            AnalSite->SetType(RxnToConsider->Product[1]);
            double Cmp;
            double NeighCmp;
            Cmp = theModel->CalculateBindingEnergy(AnalSite);
            ProdSiteB_E[0] = Cmp;
            bool Pass;
            Pass = 1;
            if (Cmp > 0)
            {
              for (int i = 0; i < AnalSite->SurroundSize; ++i)
              {
                GridSite *Neighbor = AnalSite->Surround[i];
                if (Neighbor->GetType() != theNULLSpecies)
                {
                  NeighCmp = theModel->CalculateBindingEnergy(Neighbor, true);
                  if (NeighCmp < 0)
                  {
                    ProdSiteB_E[0] = -1;
                    Pass = 0;
                    break;
                  }
                }
              }
            }
            else
            {
              Pass = 0;
            }
            if (Pass == 1)
            {
              ProdSiteB[0] = AnalSite;
            }
            AnalSite->SetType(theNULLSpecies);
            ProdSiteA[0] = NULL;
            ProdSiteA_E[0] = 0.;
            ReactSiteB[0] = NULL;
            ReactSiteB_E[0] = 0.;
            ReactSiteA[0] = AnalSite;
            ReactSiteA_E[0] = 0.;
            if(ProdSiteB_E[0] > 0.)
                Cntr = 1;
            break;
        }
        case Des2ToGas:              // type VV
        {
            ReactSiteB[0] = AnalSite;
            ReactSiteB_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteA[0] = NULL;
            ReactSiteA_E[0] = 0.;
            ProdSiteA[0] = NULL;
            ProdSiteA_E[0] = 0.;
            ProdSiteB[0] = NULL;
            ProdSiteB_E[0] = 0.;
            Cntr = 1;
            break;
        }
        case Des2:             // type AA
        {
            // Recode to include only the neighboring sites connected to the same metal atom
            //cout << "Entering Des2" << endl;
            for(ni = SurrSize; ni--;)
            {
                Neighbor = AnalSite->Surround[ni];
                if(Neighbor->Type == RxnToConsider->Reactant[1] &&
                    Neighbor != AnalSite && AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff0)
                {
                    ReactSiteA[Cntr] = AnalSite;
                    ReactSiteA_E[Cntr] = theModel->CalculateBindingEnergy(AnalSite,true);
                    ProdSiteA[Cntr] = NULL;
                    ProdSiteA_E[Cntr] = 0.;
                    ProdSiteB[Cntr] = NULL;
                    ProdSiteB_E[Cntr] = 0.;
                    ReactSiteB[Cntr] = Neighbor;
                    ReactSiteB_E[Cntr] = theModel->CalculateBindingEnergy(Neighbor,true);
                    ++Cntr;
                    assert(Cntr < SITE_NUM_SCENARIOS);
                }
            }

            break;
        }
        case EleyAdEleyDes:    // type VA
        {
            Eshort tor = AnalSite->Orientation;
            ReactSiteB[0] = AnalSite;
            ReactSiteB_E[0] = theModel->CalculateBindingEnergy(AnalSite,true);
            ReactSiteA[0] = NULL;
            ReactSiteA_E[0] = 0.;
            ProdSiteA[0] = NULL;
            ProdSiteA_E[0] = 0.;
            
            AnalSite->SetType(theNULLSpecies);
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[0],
            ProdSiteB_E[0], cutoff2, false);
            
            AnalSite->Type = RxnToConsider->Reactant[1];
            AnalSite->Orientation = tor;
            //AnalSite->BindingEnergy = theModel->CalculateBindingEnergy(AnalSite, true);
            theModel->CalculateBindingEnergy(AnalSite, true);
            if (Found)
              Cntr = 1;
            break;
        }
        case EleyAd:           // type VA
        {
            ReactSiteA[Cntr] = NULL;
            ReactSiteA_E[Cntr] = 0.;
            ReactSiteB[Cntr] = AnalSite;
            ReactSiteB_E[Cntr] = theModel->CalculateBindingEnergy(AnalSite, true);
            AnalSite->Type = theNULLSpecies;
            int tor1 = AnalSite->Orientation;
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[Cntr],
            ProdSiteB_E[Cntr], cutoff2, false);
            if (Found == 1)
            {
              ProdSiteB[Cntr]->Type = RxnToConsider->Product[1];
              theModel->CalculateBindingEnergy(ProdSiteB[Cntr]);
              Found = FindBestSingleSite(ProdSiteB[Cntr], RxnToConsider->Product[0], ProdSiteA[Cntr],
              ProdSiteA_E[Cntr], cutoff1, false);
              ProdSiteB[Cntr]->Type = theNULLSpecies;
            }
            
            if (Found == 1)
              Cntr = 1;
            AnalSite->Orientation = tor1;
            AnalSite->Type = RxnToConsider->Reactant[1];
            theModel->CalculateBindingEnergy(AnalSite, true);
            break;
        }
        case CompAd:           // type VA
        {
            ReactSiteA[Cntr] = NULL;
            ReactSiteA_E[Cntr] = 0.;
            ReactSiteB[Cntr] = NULL;
            ReactSiteB_E[Cntr] = 0.;
            
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[0],
            ProdSiteB_E[0], cutoff2, false);
            if (Found == 1)
            {
              ProdSiteB[0]->Type = RxnToConsider->Product[1];
              theModel->CalculateBindingEnergy(ProdSiteB[0]);
              Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[0],
              ProdSiteA_E[0], cutoff1, false);
              ProdSiteB[0]->Type = theNULLSpecies;
            }
            if (Found == 1)
              Cntr = 1;
            break;
        }
        case CompDes:          // type AA
        {
            for(ni = SurrSize; ni--;)
            {
                Neighbor = AnalSite->Surround[ni];
                if(Neighbor->Type == RxnToConsider->Reactant[1] &&
                    Neighbor != AnalSite && AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff0)
                {
                    ProdSiteB[Cntr] = NULL;
                    ProdSiteB_E[Cntr] = 0.;
                    ProdSiteA[Cntr] = NULL;
                    ProdSiteA_E[Cntr] = 0.;
                    ReactSiteA[Cntr] = AnalSite;
                    ReactSiteA_E[Cntr] = theModel->CalculateBindingEnergy(AnalSite, true);
                    ReactSiteB[Cntr] = Neighbor;
                    ReactSiteB_E[Cntr] = theModel->CalculateBindingEnergy(Neighbor, true);
                    ++Cntr;
                    assert(Cntr < SITE_NUM_SCENARIOS);
                }
            }
            break;
        }
        case EleyDes:          // type AA
        {
            bool Found;
            AnalSite->Type = theNULLSpecies;    // should be equal to reactant[0]
            for(ni = SurrSize; ni--;)
            {
                Neighbor = AnalSite->Surround[ni];
                if(Neighbor->Type == RxnToConsider->Reactant[1] &&
                    Neighbor != AnalSite && AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff0)
                {
                    // Check for a favorable product state
                    Neighbor->Type = theNULLSpecies;
                    Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1],
                        BestSiteB[0], BestEnergyB[0], cutoff1, true); // DWDWDW
                    if(Found == 1 && BestEnergyB[0] > 0.)
                    {
                        BestSiteB[0]->Type = theNULLSpecies;
                        ProdSiteB[Cntr] = BestSiteB[0];
                        ProdSiteB_E[Cntr] = BestEnergyB[0];
                        ProdSiteA[Cntr] = NULL;
                        ProdSiteA_E[Cntr] = 0.;
                        ReactSiteA[Cntr] = AnalSite;
                        ReactSiteA_E[Cntr] = theModel->CalculateBindingEnergy(AnalSite, true);
                        ReactSiteB[Cntr] = Neighbor;
                        ReactSiteB_E[Cntr] = theModel->CalculateBindingEnergy(Neighbor, true);
                        ++Cntr;
                        assert(Cntr < SITE_NUM_SCENARIOS);
                    }
                    Neighbor->Type = RxnToConsider->Reactant[1];
                }
            }
            AnalSite->Type = RxnToConsider->Reactant[0];
            break;
        }
        case React:            // type AA
        {
            bool Found;
            double Cmp1, Cmp2;
            Cmp1 = theModel->CalculateBindingEnergy(AnalSite, true);
            int tor1 = AnalSite->Orientation;
            AnalSite->Type = theNULLSpecies;    // should be equal to reactant[0]
            for(ni = SurrSize; ni--;)
            {
                Neighbor = AnalSite->Surround[ni];
                if(Neighbor->Type == RxnToConsider->Reactant[1] &&
                    Neighbor != AnalSite && AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff0)
                {
                    // Check for a favorable product state
                    int tor2 = Neighbor->Orientation;
                    Cmp2 = theModel->CalculateBindingEnergy(Neighbor, true);
                    Neighbor->Type = theNULLSpecies;
                    Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0],
                        BestSiteA[0], BestEnergyA[0], cutoff0, true); // DWDWDW
                    if(Found == 1 && BestEnergyA[0] > 0.)
                    {
                        BestSiteA[0]->Type = theNULLSpecies;
                        ProdSiteA[Cntr] = BestSiteA[0];
                        ProdSiteA_E[Cntr] = BestEnergyA[0];
                        ProdSiteB[Cntr] = NULL;
                        ProdSiteB_E[Cntr] = 0.;
                        ReactSiteA[Cntr] = AnalSite;
                        ReactSiteA_E[Cntr] = Cmp1;
                        ReactSiteB[Cntr] = Neighbor;
                        ReactSiteB_E[Cntr] = Cmp2;
                        ++Cntr;
                        assert(Cntr < SITE_NUM_SCENARIOS);
                    }
                    Neighbor->Orientation = tor2;
                    Neighbor->Type = RxnToConsider->Reactant[1];
                    theModel->CalculateBindingEnergy(Neighbor, true);
                }
            }
            AnalSite->Orientation = tor1;
            AnalSite->Type = RxnToConsider->Reactant[0];
            theModel->CalculateBindingEnergy(AnalSite, true);
            break;
        }
        case Diss:             // type AV
        {
            //cout << "hi" << " ";
//            if (theSimulation->SurfaceCharge == 0 && AnalSite->GetType()->GetName() == "NH3")
//            {
//              theSimulation->Univ = 7;
//            }
            Eshort tor = AnalSite->Orientation;
            ReactSiteA[Cntr] = AnalSite;
//            if (theSimulation->SurfaceCharge == 0 && AnalSite->GetType()->GetName() == "NH3")
//            {
//              cout << theSimulation->Univ << endl;
//            }
            ReactSiteA_E[Cntr] = theModel->CalculateBindingEnergy(AnalSite, true);;
            ReactSiteB[Cntr] = NULL;
            ReactSiteB_E[Cntr] = 0.;
            AnalSite->Type = theNULLSpecies;
            theSimulation->Univ = 13;
            Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[1], ProdSiteB[0],
                  ProdSiteB_E[0], cutoff2, false);
            theSimulation->Univ = 0;
            if (Found)
            {
              ProdSiteB[0]->Type = RxnToConsider->Product[1];
              theModel->CalculateBindingEnergy(ProdSiteB[0]);
              Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[0],
                    ProdSiteA_E[0], cutoff0, false);
              ProdSiteB[0]->Type = theNULLSpecies;
            }
            
            if (Found == 1)
            {
              Cntr = 1;
//              if (theSimulation->SurfaceCharge == 0)
//              {
//                cout << RxnToConsider->GetName() << " " << RxnToConsider->Product[0]->GetName() << " " << ProdSiteA_E[0] << " " << RxnToConsider->Product[1]->GetName() << " " << ProdSiteB_E[0] << endl;
//              }
            }
            theSimulation->Univ = -1;
            AnalSite->SetType(RxnToConsider->Reactant[0]);
            AnalSite->Orientation = tor;
            theModel->CalculateBindingEnergy(AnalSite, true);
            //cout << "hey" << endl;
            break; 
        }
        case Disp:             // type AA
        {
            bool Found;
            int tor1 = AnalSite->Orientation;
            double Cmp = theModel->CalculateBindingEnergy(AnalSite, true);
            AnalSite->Type = theNULLSpecies;    // should be equal to reactant[0]
            for(ni = SurrSize; ni--;)
            {
                Neighbor = AnalSite->Surround[ni];
                if(Neighbor->Type == RxnToConsider->Reactant[1] &&
                    Neighbor != AnalSite && AnalSite->FindDistance(Neighbor, theGrid->GridLength) < cutoff1)
                {
                    // Check for a favorable product state
                    double NeighCmp = theModel->CalculateBindingEnergy(Neighbor, true);
                    int tor2 = Neighbor->Orientation;
                    Neighbor->Type = theNULLSpecies;
                    Found = FindBestSingleSite(AnalSite, RxnToConsider->Product[0], ProdSiteA[Cntr],
                              ProdSiteA_E[Cntr], cutoff2, false);
                    if (Found == 1)
                    {
                      ProdSiteA[Cntr]->Type = RxnToConsider->Reactant[0];
                      Found = FindBestSingleSite(Neighbor, RxnToConsider->Product[1], ProdSiteB[Cntr],
                                ProdSiteB_E[Cntr], cutoff2, false);
                      ProdSiteA[Cntr]->Type = theNULLSpecies;
                    }
                    if(Found == 1)
                    {
                        ReactSiteA[Cntr] = AnalSite;
                        ReactSiteA_E[Cntr] = Cmp;
                        ReactSiteB[Cntr] = Neighbor;
                        ReactSiteB_E[Cntr] = NeighCmp;
                        ++Cntr;
                        assert(Cntr < SITE_NUM_SCENARIOS);
                    }
                    Neighbor->Orientation = tor2;
                    Neighbor->Type = RxnToConsider->Reactant[1];
                    theModel->CalculateBindingEnergy(Neighbor, true);
                }
            }
            AnalSite->Orientation = tor1;
            AnalSite->Type = RxnToConsider->Reactant[0];
            theModel->CalculateBindingEnergy(AnalSite, true);
            break;
        }
        default: break;  // DW 8/24/04
    }
    //cout << "Hey" << endl;
    return Cntr;
}

double ScenarioBuilder::CalculateRate(ElementaryRxn * RxnToConsider,
    double RA, double RB, double PA, double PB, double &Barrier)
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    double QA, QB, QC, QAB, QBC, Resist, D;
    double eps = 0.001;
    double delH;
    double ForwardBarrier = 0.;
    double ReverseBarrier = 0.;
    double ForwardRate = 0.;
    Process ThisProcess = RxnToConsider->GetProcess();
    double Preexponential = RxnToConsider->GetUpdPreexponential();
    double EaCorrection = RxnToConsider->GetActivationCorrection();
    double Factor;
    bool UserSet = theSimulation->UserSet;
    if (theSimulation->NoInteraction == false){
        switch (ThisProcess)
        {
            case DissToGas:
                
                D = RxnToConsider->BondDissociationEnergy;
                delH = RA + D - PB;
                ForwardBarrier = delH + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case ReactFromGas:

                Factor = RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());
                
                Preexponential = Preexponential * Factor;
                
                if(Preexponential < 1.e-50)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = RB - PA + D;
                if(PA < eps)
                    break;
                ForwardBarrier = delH + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Des1:

                ForwardBarrier = RA + EaCorrection;
                //cout << RxnToConsider->GetName() << " " << ForwardBarrier << " " << RA << endl;
                ReverseBarrier = EaCorrection;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
//                if (1)
//                {
//                  cout << RxnToConsider->GetName() << " Charge: " << theSimulation->SurfaceCharge << " Ea: " << ForwardBarrier << " RA: " << RA << " EACorrection: " << EaCorrection << " Preexp: " << Preexponential << endl;
//                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
//                cout << RxnToConsider->GetName() << " Rate: " << ForwardRate << " Barrier: " << ForwardBarrier << endl;
                break;

            case Ad1:
                Factor = RxnToConsider->Reactant[0]->CalculateAdsorptionRate(theSimulation->GetTemperature());
                Preexponential = Preexponential * Factor;
                if(Preexponential < 1.e-50)
                    break;
                if(PA < eps)
                    break;
                
                ReverseBarrier = PA + EaCorrection;
                ForwardBarrier = EaCorrection;
//                if (RxnToConsider->Reactant[0]->GetName() == "N2")
//                {
//                  cout << " Barrier: " << ForwardBarrier << " PA: " << PA << " BE: " << RxnToConsider->Reactant[0]->BindingEnergy[1][0][ind];
//                }
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
//                if (RxnToConsider->Reactant[0]->GetName() == "N2")
//                {
//                  cout << " Preexp: " << Preexponential << " Barrier: " << ForwardBarrier << " Rate: " << ForwardRate << endl;
//                }
                break;

            case EleyAdEleyDes:

                QA = RB;
                QBC = RxnToConsider->Reactant[0]->StableSiteEnergy[ind];
                if(QBC < 0.)
                    QBC = 0.;
                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                QAB = RxnToConsider->Product[0]->StableSiteEnergy[ind];
                QC = PB;
                if(QC < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = QA + D - QC;
                if(D >= 0.)
                {
                    Resist = (QC * QAB) / (QC + QAB);
                }
                else
                {
                    Resist = (QA * QBC) / (QA + QBC);
                }
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case CompAd:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[1]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                QA = PA;
                if(QA < eps)
                    break;
                QB = PB;
                if(QB < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = D - QA - QB;
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case CompDes:

                QA = RA;
                QB = RB;

                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = QA + QB + D;
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case EleyAd:

                QA = RB;
                QBC = RxnToConsider->Reactant[0]->StableSiteEnergy[ind];
                if(QBC < 0.)
                    QBC = 0.;
                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                QAB = PA;
                if(QAB < eps)
                    break;
                QC = PB;
                if(QC < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = QA + D - QC - QAB;
                if(D >= 0.)
                {
                    Resist = (QC * QAB) / (QC + QAB);
                }
                else
                {
                    Resist = (QA * QBC) / (QA + QBC);
                }
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case EleyDes:

                QA = RB;
                QBC = RA;

                QAB = RxnToConsider->Product[0]->StableSiteEnergy[ind];
                if(QAB < 0.)
                    QAB = 0.;
                QC = PB;
                if(QC < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = QA + QBC + D - QC;
                if(D >= 0.)
                {
                    Resist = (QC * QAB) / (QC + QAB);
                }
                else
                {
                    Resist = (QA * QBC) / (QA + QBC);
                }
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case Disp:

                QA = RB;
                QBC = RA;

                QAB = PA;
                if(QAB < eps)
                    break;
                QC = PB;
                if(QC < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = QA + QBC + D - QC - QAB;
                if(D >= 0.)
                {
                    Resist = (QC * QAB) / (QC + QAB);
                }
                else
                {
                    Resist = (QA * QBC) / (QA + QBC);
                }
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case React:

//                if (RxnToConsider->_id == 11)
//                {
//                  cout << RA << " " << RB << " " << PA << endl;
//                }
                QA = RA;
                QB = RB;
                QC = PA;
                if(QC < eps)
                {
                    break;
                }
                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = QA + QB + D - QC;    // negative is right for D, i think, check Reactions.c for sign convention
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }      
                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential,ForwardBarrier);
//                if (RxnToConsider->GetName() == "H * + NH * = NH2 *")
//                {
//                  cout << RxnToConsider->GetName() << " Forward Rate: " << ForwardRate << " Preexp: " << Preexponential << " Barrier: " << ForwardBarrier << endl;
//                }
                break;

            case Diss:

                QC = RA;

                QA = PB;
                if(QA < eps)
                    break;

                QB = PA;
                if(QB < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = D + QC - QA - QB;
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                
//                if (theSimulation->SurfaceCharge == 0)
//                  cout << RxnToConsider->GetName() << " Preexp: " << Preexponential << " Barrier: " << ForwardBarrier << " Forward Rate: " << ForwardRate << 
//                  " D: " << D << 
//                  " QC: " << QC << " " << RxnToConsider->Reactant[0]->StableSiteEnergy[ind] << " QA: " << QA << " " << RxnToConsider->Product[1]->StableSiteEnergy[ind] << " QB: " << QB << " " << RxnToConsider->Product[0]->StableSiteEnergy[ind] << endl;
                // No statistical factor for dissociation, either forward or reverse
                break;

            case Des2:
                QA = RA;
                QB = RB;
                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = QA + QB + D;
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case Ad2:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                QA = PA;
                if(QA < eps)
                    break;
                QB = PB;
                if(QB < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                Resist = (QA * QB) / (QA + QB);
                delH = D - QA - QB;
                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
//                if (RxnToConsider->GetName() == "N2 = N * + N *")
//                  cout << RxnToConsider->GetName() << " D: " <<  D << " QA: " << QA << " QB: " << QB << " dH: " << delH << " Barrier: " << ForwardBarrier << " Preexp: " << Preexponential << endl;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case Ad2ToGas:
                Preexponential = Preexponential * RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());
                QB = PB;
                if(QB < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = D - QB;
                ForwardBarrier = delH + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case Des2ToGas:
                Preexponential = Preexponential * RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());
                QB = RB;
                if(QB < eps)
                    break;
                D = RxnToConsider->BondDissociationEnergy;
                delH = D + QB;
                ForwardBarrier = delH + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case ElecTransfer:
            case Diff:

                QA = RA;
                QB = PA;
                Resist = (QA * QB) / (QA + QB);     // these equations work surprisingly well for diffusion 
                delH = QA - QB;

                ForwardBarrier = 0.5 * (delH + Resist) + EaCorrection;
                ReverseBarrier = ForwardBarrier - delH;
                if (UserSet == false)
                {
                  if(ForwardBarrier < 0.)
                  {
                      ForwardBarrier = 0.;
                      ReverseBarrier = delH;
                  }
                  if(ReverseBarrier < 0.)
                  {
                    ForwardBarrier -= ReverseBarrier;
                    ReverseBarrier = 0;
                  }
                }
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            default:
                cout << "\nUndefined method in scenario!" << endl;
        }
    }
    else if (theSimulation->NoInteraction == true)
    {
    RxnToConsider->ActivationEnergy = RxnToConsider->UserSetActivationEnergy;
    ForwardBarrier = RxnToConsider->ActivationEnergy;
    switch (ThisProcess)
        {
            case DissToGas:
                
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case ReactFromGas:

                Preexponential = Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Des1:

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Ad1:

                Preexponential = Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());
                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case EleyAdEleyDes:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case EleyAd:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case EleyDes:

                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case Disp:

                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case React:

                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Diss:

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Des2:

                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case CompDes:

                ForwardRate =
                    RxnToConsider->StatFactor * SurfaceRate(Preexponential,
                    ForwardBarrier);
                break;

            case Des2ToGas:

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Ad2ToGas:

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            case Ad2:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case CompAd:

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[0]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                Preexponential =
                    Preexponential *
                    RxnToConsider->Reactant[1]->
                    CalculateAdsorptionRate(theSimulation->GetTemperature());

                ForwardRate = RxnToConsider->StatFactor * SurfaceRate(Preexponential, ForwardBarrier);      // Ad2 has a statistical factor of 0.5
                break;

            case ElecTransfer:
            case Diff:

                ForwardRate = SurfaceRate(Preexponential, ForwardBarrier);
                break;

            default:
                cout << "\nUndefined method in scenario!" << endl;}
    }
    RxnToConsider->AvActivationEnergy = RxnToConsider->AvActivationEnergy + ForwardBarrier;
//    if (RxnToConsider->GetName() == "N2 = N * + N *")
//    {
//      cout << RxnToConsider->GetName() << " " << delH << " " << ForwardBarrier << endl;
//      exit(1);
//    }
    Barrier = ForwardBarrier;
    return ForwardRate;
}

double ScenarioBuilder::SurfaceRate(double A, double E)
{
    // DW 8/24/04: ddd is probably for debugging purposes.
    // added the void line to make it used in this function.
    double ddd = theSimulation->GetTemperature();
    (void)ddd;
    return A * exp(-E / (R_GasConst * theSimulation->GetTemperature()));
}

