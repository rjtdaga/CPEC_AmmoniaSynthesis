#include "Scenario.h"
#include "ScenarioBuilder.h"
#include "Base/Random.h"
#include "Species.h"
#include "ShortRoutines.h"
#include "Model.h"

Scenario *thisScenario;
Scenario *thisEquilScenario;

Scenario::Scenario()
{
    ReactantSites[0] = NULL;
    ReactantSites[1] = NULL;
    ProductSites[0] = NULL;
    ProductSites[1] = NULL;
    Reactants[0] = theNULLSpecies;
    Reactants[1] = theNULLSpecies;
    Products[0] = theNULLSpecies;
    Products[1] = theNULLSpecies;
    theReaction = NULL;
    Possible = false; // DWDWDW
    ForwardRate = 0.;
    ReverseRate = 0.;
    ForwardBarrier = 0.;
    ReverseBarrier = 0.;
    Preexponential = 0.;
    ForwardPreexponential = 0.;
}

// copy Constructor
Scenario::Scenario(const Scenario & S)
{

    ReactantSites[0] = S.ReactantSites[0];
    ReactantSites[1] = S.ReactantSites[1];
    ProductSites[0] = S.ProductSites[0];
    ProductSites[1] = S.ProductSites[1];
    Reactants[0] = S.Reactants[0];
    Reactants[1] = S.Reactants[1];
    Products[0] = S.Products[0];
    Products[1] = S.Products[1];
    theReaction = S.theReaction;
    Possible = S.Possible;
    ForwardRate = S.ForwardRate;
    ReverseRate = S.ReverseRate;
    ForwardBarrier = S.ForwardBarrier;
    ReverseBarrier = S.ReverseBarrier;
    Preexponential = S.Preexponential;
    ForwardPreexponential = S.ForwardPreexponential;
    ReversePreexponential = S.ReversePreexponential;
}

void Scenario::SetReaction(ElementaryRxn * thisRxn)
{
    if(thisRxn)
    {
        theReaction = thisRxn;
        theReaction->GetParticipators(Reactants, Products);
    }
    return;
}

Process Scenario::GetProcess()
{
    if(theReaction)
    {
        return theReaction->GetProcess();
    }
    else
        return NoProc;
}

void Scenario::SetUpScenario(const ScenarioAction & S)
{

    ReactantSites[0] = S.ReactantSites[0];
    ReactantSites[1] = S.ReactantSites[1];
    ProductSites[0] = S.ProductSites[0];
    ProductSites[1] = S.ProductSites[1];
    SetReaction(S.Rxn);
    ForwardRate = S.ForwardRate;
}


bool Scenario::Screen(GridSite * PassedSite)
{
    // Conditions that involve only the first site

    switch (GetProcess())
    {
        case NoProc:
        {
            return false;  // DWDWDW
        }
        case Ad1:
        {
            if(PassedSite->Type == theNULLSpecies && PassedSite->BondSize == 1)
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50)
                  {
                      return true; // DWDWDW
                  }
                }
                if (theSimulation->Solvent == true)
                {
                  if(Reactants[0]->GetConcentration() > 1.e-50)
                  {
                      return true; // DWDWDW
                  }
                }
            }
            return false; // DWDWDW
        }
        case ReactFromGas:
        {
            if(PassedSite->Type == Reactants[1])
            {
              if (theSimulation->Solvent == true)
              {
                if (Reactants[0]->GetConcentration() > 1.e-50)
                {
                  return true; // DWDWDW
                }
              }
              if (theSimulation->Gas == true)
              {
                if (Reactants[0]->GetPressure() > 1.e-50)
                {
                  return true; // DWDWDW
                }
              }
            }
            return false; // DWDWDW
        }
        case Ad2:
        {
            if(PassedSite->Type == theNULLSpecies
             && PassedSite->BondSize == 1 )
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50)
                    return true; // DWDWDW
                }
                if (theSimulation->Solvent == true)
                {
                  if(Reactants[0]->GetConcentration() > 1.e-50)
                    return true; // DWDWDW
                }
            }
            return false; // DWDWDW
        }
        case Ad2ToGas:
        {
            if(PassedSite->Type == theNULLSpecies)
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50)
                    return true; // DWDWDW
                }
                if (theSimulation->Solvent == true)
                {
                  //cout << Reactants[0]->GetConcentration() << endl;
                  if(Reactants[0]->GetConcentration() > 1.e-50)
                    return true; // DWDWDW
                }
            }
            return false; // DWDWDW
        }
        case Des2ToGas:
        {
            if(PassedSite->Type == Reactants[1])
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50)
                    return true; // DWDWDW
                }
                if (theSimulation->Solvent == true)
                {
                  if(Reactants[0]->GetConcentration() > 1.e-50)
                    return true; // DWDWDW
                }
            }
            return false; // DWDWDW
        }
        case ElecTransfer:
        case Diff:
        {
            if(PassedSite->Type == Reactants[0])
            {
                return true; // DWDWDW
            }
            return false; // DWDWDW
        }
        case Diss:
        {
            if(PassedSite->Type == Reactants[0])
                return true; // DWDWDW
            return false; // DWDWDW
        }
        case CompAd:
        {
            if(PassedSite->Type == theNULLSpecies)
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50 && Reactants[1]->GetPressure() > 1.e-50)
                    return true; // DWDWDW
                }
                if (theSimulation->Solvent == true)
                {
                  if(Reactants[0]->GetConcentration() > 1.e-50 && Reactants[1]->GetConcentration() > 1.e-50)
                    return true; // DWDWDW
                }
            }
            return false; // DWDWDW
        }
        case EleyAd:
        case EleyAdEleyDes:
        {
            if(PassedSite->Type == Reactants[1])
            {
                if (theSimulation->Gas == true)
                {
                  if(Reactants[0]->GetPressure() > 1.e-50)
                    return true; // DWDWDW
                }
                if (theSimulation->Solvent == true)
                {
                  if(Reactants[0]->GetConcentration() > 1.e-50)
                    return true; // DWDWDW
                }
            }
            return false; // DWDWDW
        }
        case CompDes:
        case Des2:
        case React:
        case Disp:
        case EleyDes:
        case Des1:
        case DissToGas:
        {
            if(PassedSite->Type == Reactants[0])
            {
                return true; // DWDWDW
            }
            return false; // DWDWDW
        }
        default:
        {
            return false; // DWDWDW
        }
    }
}

int Scenario::Equilibrate()
{

    double Keq = ForwardRate / ReverseRate;

    if(Keq > 1.)
    {
        doit();
        return 1;
    }
    double Prob = FindRandomNumber();

    if(Prob < Keq)
    {
        doit();
        return 1;
    }
    return 0;
}

void Scenario::doit()
{
    double SC = theSimulation->SurfaceCharge;
    thisComponent = NULL;
    //cout << theReaction->GetName() << endl;
    switch (GetProcess())
    {
        case Des1:
        {
            Reactants[0]->TotalDesorbed += 1;
            break;
        }
        case EleyDes:
        case Des2:
        case Des2ToGas:
        case Ad2ToGas:
        case DissToGas:
        case EleyAdEleyDes:
        {
            Products[0]->TotalDesorbed += 1;
            break;
        }
        case CompDes:
        {
            Products[0]->TotalDesorbed += 1;
            Products[1]->TotalDesorbed += 1;
            break;
        }

        default: break;  // DW 8/24/04
    }
    if(ReactantSites[0])
    {
        thisComponent = ReactantSites[0]->GetType();
//        if (thisComponent->GetName() == "NH3")
//        {
//          cout << thisComponent->GetName() << " BE: " << theModel->CalculateBindingEnergy(ReactantSites[0], true) << " Barrier: " << ForwardBarrier << endl;
//        }
        if(thisComponent != theNULLSpecies)
        {
            OrigReactOrientation[0] = ReactantSites[0]->Orientation;
            thisComponent->NumberOfAtoms -= 1;
            ReactantSites[0]->SetType(theNULLSpecies);
        }
        
    }
    if(ReactantSites[1])
    {
        thisComponent = ReactantSites[1]->GetType();
        if(thisComponent != theNULLSpecies)
        {
            OrigReactOrientation[1] = ReactantSites[1]->Orientation;
            thisComponent->NumberOfAtoms--;
            ReactantSites[1]->SetType(theNULLSpecies);  // Do not attempt to reset binding energy
        }
    }
    //cout << Products[0]->GetName() << endl;
    if(ProductSites[0])
    {
//        Eshort c[3];
//        ProductSites[0]->GetGridPos(c);
//        cout << c[0] << " " << c[1] << endl;
//        cout << "Entered" << endl;
        ProductSites[0]->SetType(Products[0]);
        if(Products[0] != theNULLSpecies)
        {
            //theSimulation->Univ = 5;
            //cout << Products[0]->GetName() << " ";
            double a = theModel->CalculateBindingEnergy(ProductSites[0]);
            //cout << a << endl;
            //theSimulation->Univ = 0;
            Products[0]->NumberOfAtoms += 1;
        }
    }
    if(ProductSites[1])
    {
        ProductSites[1]->SetType(Products[1]);
        if(Products[1] != theNULLSpecies)
        {
            //cout << Products[1]->GetName() << " ";
            double b = theModel->CalculateBindingEnergy(ProductSites[1]);
            //cout << b << endl;
            Products[1]->NumberOfAtoms += 1;
        }
    }
    return;
}


void Scenario::undoit()
{

    thisComponent = NULL;
    switch (GetProcess())
    {
        case Des1:
        {
            Reactants[0]->TotalDesorbed -= 1;
            break;
        }
        case EleyDes:
        case Des2:
        case Des2ToGas:
        case Ad2ToGas:
        case DissToGas:
        case EleyAdEleyDes:
        {
            Products[0]->TotalDesorbed -= 1;
            break;
        }
        case CompDes:
        {
            Products[0]->TotalDesorbed -= 1;
            Products[1]->TotalDesorbed -= 1;
            break;
        }
        default: break;  // DW 8/24/04
    }
    if(ProductSites[0])
    {
        ProductSites[0]->SetType(theNULLSpecies);
        if(Products[0] != theNULLSpecies)
        {
            Products[0]->NumberOfAtoms -= 1;
        }
    }
    if(ProductSites[1])
    {
        ProductSites[1]->SetType(theNULLSpecies);
        if(Products[1] != theNULLSpecies)
        {
            Products[1]->NumberOfAtoms -= 1;
        }
    }
    if(ReactantSites[0])
    {
        thisComponent = ReactantSites[0]->GetType();
        if(thisComponent == theNULLSpecies)
        {
            thisComponent->NumberOfAtoms += 1;
            ReactantSites[0]->SetType(Reactants[0]);
        }
        else
            ReactantSites[0]->Orientation = OrigReactOrientation[0];
    }
    if(ReactantSites[1])
    {
        thisComponent = ReactantSites[1]->GetType();
        if(thisComponent == theNULLSpecies)
        {
            thisComponent->NumberOfAtoms++;
            ReactantSites[1]->SetType(Reactants[1]);    // Do not attempt to reset binding energy
        }
        else
            ReactantSites[1]->Orientation = OrigReactOrientation[1];
    }
    return;
}

Scenario::~Scenario()
{
} 
