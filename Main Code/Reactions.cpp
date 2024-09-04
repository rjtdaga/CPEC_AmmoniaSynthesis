

#include <iomanip>
#include <cstdlib>
#include <cassert>
using namespace std;


#include "Component.h"
#include "Reactions.h"
#include "Base/IntRandom.h"
#include "ElementaryRxn.h"
#include "Species.h"
#include "Geometry.h"
Reactions *theReactions = NULL;

SeqList<ElementaryRxn> *Rxns = NULL;

Reactions::Reactions(fstream &fin, fstream &fout) : Facilitator(fin, fout)
{
    // Initialize Input and Output File
    ::theReactions = this; // playing with global indicator (fast solution) 
    Initialize((CharString) "Reactions", fin); // File Facilitator
    // Read from File
    NumberOfReactions = Read(fin);
    Rxns = &theReactions;
//    Rxns->ResetToFront();
//    while (*Rxns)
//    {
//      ElementaryRxn *thisReaction = Rxns->GetPtr();
//      cout << thisReaction->GetName() << " " << thisReaction->GetProcessName() << " " << thisReaction->UserSetActivationEnergy << " " << thisReaction->GetPreexponential() << endl;
//      ++(*Rxns);
//    }
//    exit(1);
    Listing = new Component *[theSpecies->GetSize() + 1];
    Initialize((CharString) "SIM_SPECIES", fin);
    CharString End = "END";
    CharString Next;
    
    while(fin >> Next && Next.lower() != End.lower())
    {
        Component *KeyComponent = theSpecies->get(Next);
        KeyComponent->Consider = true;  // DWDWDW
    }
    NumberOfSpecies = FindSurfaceSpecies();
    
    theSpecies->theGeometries->theCoordinates.ResetToFront();
    SeqList < Coordinates > NewCoordinates;
    Coordinates SiteCoords;
    Component *thisType;
    
    while(theSpecies->theGeometries->theCoordinates)
    {
        SiteCoords = theSpecies->theGeometries->theCoordinates();
        bool Continue = true;   // DWDWDW  // DWDWDW
        
        if(theSpecies->get(SiteCoords.SpeciesName)->Consider)
        {
//            for(int i = 0; i < SiteCoords.NumAtoms && Continue; ++i)
//            {
//                thisType = theSpecies->get(SiteCoords.AtomNames[i]);
//                if(thisType->_idhook >= 0)
//                {
//                    SiteCoords.AtomNames[i] = thisType->_idhook;
//                }
//                else
//                {
//                    Continue = false;   // DWDWDW
//                }
//            }
            Continue = true;
        }
        else
        {
            Continue = false;   // DWDWDW
        }
        
        if(Continue)
        {
            NewCoordinates.Add(SiteCoords);
        }
        ++theSpecies->theGeometries->theCoordinates;
    }
    
    theSpecies->theGeometries->theCoordinates.RemoveAll();
    NewCoordinates.ResetToFront();
    while(NewCoordinates)
    {
        theSpecies->theGeometries->theCoordinates.Add(NewCoordinates());
        //cout << NewCoordinates().SpecName << " " << NewCoordinates().NumAtoms << " " << NewCoordinates().MetalCoord << endl;
        ++NewCoordinates;
    }
    //exit(1);
}


void Reactions::SetBondDissociationEnergies()
{
    Rxns->ResetToFront();
    int setid = 0;
    
    while(*Rxns)
    {
        ElementaryRxn *A = Rxns->GetPtr();
        int ind = theSimulation->FindChargeIndex(A->SurfaceCharge);
        Component *RPtrs[2];
        
        RPtrs[0] = theNULLSpecies;
        RPtrs[1] = theNULLSpecies;
        
        Component *PPtrs[2];
        
        PPtrs[0] = theNULLSpecies;
        PPtrs[0] = theNULLSpecies;
        
        A->GetParticipators(RPtrs, PPtrs);
        Process Type;
        Type = A->GetProcess();
        
        // Set the _id field.  Increment this to the next value.
        A->_id = setid;
        ++setid;
        switch (Type)
        {
            case Ad1:
            {
              A->BondDissociationEnergy = 0.;
              break;
            }
            case Ad2:
            {
              A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] + PPtrs[1]->StableSiteEnergy[ind];
              break;
            }
            case Des1:
            {
              A->BondDissociationEnergy = 0.;
              break;
            }
            case Des2:
            {
              A->BondDissociationEnergy = A->deltaH - RPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
              break;
            }
            case Diff:
            {
              A->BondDissociationEnergy = 0.;
              break;
            }
            case Diss:
            {
              A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[0]->StableSiteEnergy[ind];
              break;
            }
            case React:
            {
              A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] - RPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
              break;
            }
            case Disp:
            {
              A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
              break;
            }
            case EleyAd:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case EleyDes:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case EleyAdEleyDes:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case DissToGas:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[1]->StableSiteEnergy[ind] - RPtrs[0]->StableSiteEnergy[ind];
                break;
            }
            case ReactFromGas:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case ElecTransfer:
            {
                A->BondDissociationEnergy = 0.;
                break;
            }
            case NonSurface_e:
            {
                A->BondDissociationEnergy = 0.;
                break;
            }
            case Des2ToGas:
            {
                A->BondDissociationEnergy = A->deltaH - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case Ad2ToGas:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case CompDes:
            {
                A->BondDissociationEnergy = A->deltaH - RPtrs[0]->StableSiteEnergy[ind] - RPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            case CompAd:
            {
                A->BondDissociationEnergy = A->deltaH + PPtrs[0]->StableSiteEnergy[ind] + PPtrs[1]->StableSiteEnergy[ind];
                break;
            }
            default:
            {
                exit(1);
            }
        }
        ++(*Rxns);
    }
    return;
}

int Reactions::Read(fstream &fin)
{
    fin.seekg(inp);             // Seek species section of file
    CharString Next;
    CharString End = "END";
    CharString Reactant[2];
    CharString Product[2];
    CharString LHS, RHS;
    
    Component *R[2], *P[2];
    
    R[0] = theNULLSpecies;
    R[1] = theNULLSpecies;
    P[0] = theNULLSpecies;
    P[1] = theNULLSpecies;
    
    Process ThisProcess;
    ElementaryRxn NewReaction;
    ElementaryRxn RevReaction;
    ElementaryRxn *NewReactionPtr = NULL;
    ElementaryRxn *RevReactionPtr = NULL;
    int NumberFound = 0;
    bool ReactantTwo = false;   // DWDWDW  // DWDWDW
    bool ProductTwo = false;    // DWDWDW  // DWDWDW
    bool ReactantGas = false;   // DWDWDW  // DWDWDW
    bool ProductGas = false;    // DWDWDW  // DWDWDW
    double SurfCharg;
    streampos linepos;
    char c;
    
    // DW 8/24/04 This was tellp()
    linepos = fin.tellg();
    fin >> Next;
    while(Next.lower() != End.lower())
    {
        // create space for the new reactions
        fin >> SurfCharg;
        NumberFound = 0;
        while (fin >> Next && Next.lower() != "surfacecharge" && Next.lower() != End.lower())
        {
          int NumReactantGas = 0;
          int NumProductGas = 0;
          NewReactionPtr = new ElementaryRxn;
          RevReactionPtr = new ElementaryRxn;
          assert(NewReactionPtr && RevReactionPtr);
          
          // Copy the data from the newly created on to the local one
          NewReaction = *NewReactionPtr;
          RevReaction = *NewReactionPtr;
          NewReaction.SurfaceCharge = SurfCharg;
          RevReaction.SurfaceCharge = SurfCharg;
          // clear out the names
          Reactant[0] = Reactant[1] = theNULLSpecies->GetName();
          Product[0] = Product[1] = theNULLSpecies->GetName();
          
          fin.seekg(linepos);     // Back up
          
          // skip over the name of the rxn until you get to the '#'
          while(fin.get(c) && (c != '#'))
          {
              ;
          }
          
          // DW 8/24/04 This was tellp()
          // Save the line position at the start of the actual reaction
          linepos = fin.tellg();
          
          Next = "";
          LHS = "";
          RHS = "";
          
          while(fin.get(c) && (c != '='))
          {
              // Everything up to the '=' is the Reaction LHS
              LHS = LHS + c;
          }
          LHS.TrimTrailingBlanks();
          LHS.TrimLeadingBlanks();
          
          while(fin.get(c) && (c != '#'))
          {
              // Everything up to the '#' is the Reaction RHS
              RHS = RHS + c;
          }
          RHS.TrimTrailingBlanks();
          RHS.TrimLeadingBlanks();
          
          Next = LHS + " = " + RHS;
          NewReaction.SetName(Next);      // Reaction Name
          
          Next = RHS + " = " + LHS;
          RevReaction.SetName(Next);      // Reaction Name
          
          
          // Back up to re-get the reaction pieces.
          fin.seekg(linepos);     
          
          // Get Reactants
          
          fin >> Next;
          Reactant[0] = Next;
          R[0] = theSpecies->add(Reactant[0]);
          
          fin >> Next;
          
          // A * is used to denote a surface state.  resolve this character
          if(Next == "*")
          {
              ReactantGas = false;   // DWDWDW // if a *, then a surface state
              R[0]->OnSurface = true;     // DWDWDW 
              fin >> Next;
          }
          else
          {
              ReactantGas = true; // DWDWDW  // if not a *, then a gas state
              NumReactantGas = 1;
          }
          
          if(Next == "+")
          {
              // if two reactants, then a recombination
              ReactantTwo = true; // DWDWDW 
              fin >> Next;
              Reactant[1] = Next;
              R[1] = theSpecies->add(Reactant[1]);
              fin >> Next;
              if(Next == "*")
              {                   // a surface state
                  R[1]->OnSurface = true; // DWDWDW
                  fin >> Next;
              }
              else
              {
                ReactantGas = true;
                NumReactantGas = NumReactantGas + 1;
              }
          }
          else
          {
              ReactantTwo = false;        // DWDWDW 
              Reactant[1] = theNULLSpecies->GetName();
              R[1] = theSpecies->add(Reactant[1]);
          }
          
          // No matter what, Next should be '=' since we've dealt with
          // all possibilities for reactants.
          assert(Next == "=");
          
          fin >> Next;
          Product[0] = Next;
          P[0] = theSpecies->add(Product[0]);
          fin >> Next;
          if(Next == "*")
          {
              ProductGas = false; // DWDWDW
              P[0]->OnSurface = true;     // DWDWDW
              fin >> Next;
          }
          else
          {
              ProductGas = true;  // DWDWDW
              NumProductGas = 1;
          }
          
          if(Next == "+")
          {
              ProductTwo = true;  // DWDWDW
              fin >> Next;
              Product[1] = Next;
              P[1] = theSpecies->add(Product[1]);
              fin >> Next;
              if(Next == "*")
              {
                  P[1]->OnSurface = true; // DWDWDW
                  fin >> Next;
              }
              else
              {
                ProductGas = true;
                NumProductGas = NumProductGas + 1;
              }
          }
          else
          {
              ProductTwo = false; // DWDWDW 
              Product[1] = theNULLSpecies->GetName();
              P[1] = theSpecies->add(Product[1]);
          }
          
          // Now should be the '#' following the reaction
          assert(Next == "#");
          
          double Pre, Pre2, Ea, dH;
          int Importance;
          bool Equil;
          fin >> Pre >> Pre2 >> Ea >> dH >> Equil;
          NewReaction.SetPreexponential(Pre);
          NewReaction.UserEquil = Equil;
          NewReaction.SetUpdPreexponential(Pre);
          NewReaction.isReverse = false;
          NewReaction.deltaH = dH;
          RevReaction.SetPreexponential(Pre2);
          RevReaction.UserEquil = Equil;
          RevReaction.SetUpdPreexponential(Pre2);
          if (Ea - dH > 0)
          {
            NewReaction.UserSetActivationEnergy = Ea;
            RevReaction.UserSetActivationEnergy = Ea - dH;
          }
          else
          {
            NewReaction.UserSetActivationEnergy = dH;
            RevReaction.UserSetActivationEnergy = 0.;
          }
          RevReaction.deltaH = -dH;
          RevReaction.isReverse = true;
          while(fin.get(c) && (c != '\0' && c != '\n'))
          {
              Next = Next + c;
          }
  
          if(ProductTwo == false && ReactantTwo == false) // DWDWDW
          {
              if(NumProductGas == 0 && NumReactantGas == 0) // DWDWDW
              {
                 if (Reactant[0] == Product[0])
                 {
                   ThisProcess = Diff; // A* ---> A*
                 }
                 else if (R[0]->MolWt == P[0]->MolWt)
                 {
                   ThisProcess = ElecTransfer; // A* + e- ---> B*
                 }
                 else
                 {
                   cout << "The Molecular weights do not add up!!" << endl;
                   exit(1);
                 }
              }
              else if(NumProductGas == 0 && NumReactantGas == 1) // DWDWDW 
              {
                  ThisProcess = Ad1; // A ---> A*
              }
              else if (NumReactantGas == 0 && NumProductGas == 1)
              {
                  ThisProcess = Des1; // A* ---> A
              }
              else if (NumReactantGas == 1 && NumProductGas == 1)
              {
                  ThisProcess = NonSurface_e; // A + e- <---> A-
              }
          }
          else if(ProductTwo == true && ReactantTwo == false) // DWDWDW
          {
              if(NumReactantGas == 0 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = Diss; // A* ---> B* + C*
              }
              else if(NumReactantGas == 1 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = Ad2; // A ---> B* + C*
              }
              else if(NumReactantGas == 0 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = DissToGas; // A* ---> B* + C
              }
              else if(NumReactantGas == 1 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = Ad2ToGas; //  C ---> A + B*
              }
              else if(NumReactantGas == 1 && NumProductGas == 2)
              {
                  exit(1); // This would be a gas phase rxn: A ---> B + C
              }
          }
          else if(ProductTwo == false && ReactantTwo == true) // DWDWDW
          {
              if(NumReactantGas == 0 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = React; // A* + B* ---> C*
              }
              else if(NumReactantGas == 0 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = Des2; //  A* + B* ---> C
              }
              else if(NumReactantGas == 1 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = Des2ToGas; //  A + B* ---> C
              }
              else if(NumReactantGas == 1 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = ReactFromGas; //  A + B* ---> C*
              }
              else if(NumReactantGas == 2 && NumProductGas == 1)
              {
                  exit(1); // This would be a gas phase rxn: A + B ---> C
              }
          }
          else if(ProductTwo == true && ReactantTwo == true) // DWDWDW
          {
              if(NumReactantGas == 0 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = Disp; // A* + B* ---> C* + D*
              }
              else if(NumReactantGas == 1 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = EleyAdEleyDes;  // A + B* ---> C + D*
              }
              else if(NumReactantGas == 1 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = EleyAd; // A + B* ---> C* + D*
              }
              else if(NumReactantGas == 0 && NumProductGas == 1) // DWDWDW
              {
                  ThisProcess = EleyDes; // A* + B* ---> C + D*
              }
              else if(NumReactantGas == 2 && NumProductGas == 0) // DWDWDW
              {
                  ThisProcess = CompAd; // A + B ---> C* + D*
              }
              else if(NumReactantGas == 0 && NumProductGas == 2) // DWDWDW
              {
                  ThisProcess = CompDes; // A* + B* ---> C + D
              }
          }
          
          // DW 8/24/04 was tellp();
          linepos = fin.tellg();
          NewReaction.SetParticipators(R, P);
          NewReaction.SetProcess(ThisProcess);
          if(ThisProcess != Diff)
          {
              RevReaction.SetParticipators(P, R);
              RevReaction.SetProcess(NewReaction.GetReverseProcess());
              add(NewReaction);
              add(RevReaction);
              has(NewReaction)->SetReverseReaction(has(RevReaction));
              has(RevReaction)->SetReverseReaction(has(NewReaction));
              NumberFound += 2;
          }
          else
          {
              R[0]->DiffusionReaction = add(NewReaction);
              has(NewReaction)->SetReverseReaction(has(NewReaction));
              NumberFound++;
          }
        }
        NumberOfReactions = NumberFound;
    }
    return NumberFound;
}

ElementaryRxn *Reactions::add(ElementaryRxn New)
{
    ElementaryRxn *Find = has(New);

    if(Find == NULL)
    {
        theReactions.Add(New);
        return has(New);
    }
    return Find;
}

ElementaryRxn *Reactions::has(ElementaryRxn A)
{

    Component *RPtrs[2];

    RPtrs[0] = theNULLSpecies;
    RPtrs[1] = theNULLSpecies;
    Component *PPtrs[2];

    PPtrs[0] = theNULLSpecies;
    PPtrs[0] = theNULLSpecies;
    Process Type;

    A.GetParticipators(RPtrs, PPtrs);
    Type = A.GetProcess();

    Component *Stored_RPtrs[2];

    Stored_RPtrs[0] = theNULLSpecies;
    Stored_RPtrs[1] = theNULLSpecies;
    Component *Stored_PPtrs[2];

    Stored_PPtrs[0] = theNULLSpecies;
    Stored_PPtrs[1] = theNULLSpecies;
    Process StoredType;

    theReactions.ResetToFront();
    ElementaryRxn ThisRxn;

    while(theReactions)
    {
        ThisRxn = theReactions.Get();
        ThisRxn.GetParticipators(Stored_RPtrs, Stored_PPtrs);
        StoredType = ThisRxn.GetProcess();

        if(RPtrs[0] == Stored_RPtrs[0] &&
            RPtrs[1] == Stored_RPtrs[1] &&
            PPtrs[0] == Stored_PPtrs[0] &&
            PPtrs[1] == Stored_PPtrs[1] && Type == StoredType && A.SurfaceCharge == ThisRxn.SurfaceCharge)
        {

            return theReactions.GetPtr();
        }
        ++theReactions;
    }
    return NULL;
}

int Reactions::FindSurfaceSpecies()
{
    // DW: Made theSpecies private in Species.h, so modified the next line
    int arraysize = theSpecies->GetSize();

    assert(Listing != NULL);
    for(int ii = 0; ii < arraysize; ++ii)
        Listing[ii] = NULL;
        
    Component *Ri[2];

    Ri[0] = theNULLSpecies;
    Ri[1] = theNULLSpecies;
    Component *Pi[2];

    Pi[0] = theNULLSpecies;
    Pi[1] = theNULLSpecies;
    thisReaction = NULL;
    SeqList < Component * >Participators;
    Rxns->ResetToFront();
    //  Find all Participating Species (to get thetastar)
    bool Remove = false;        // DWDWDW

    while(*Rxns)
    {
        Remove = false;         // DWDWDW
        thisReaction = Rxns->GetPtr();
        thisReaction->GetParticipators(Ri, Pi);
        switch (thisReaction->GetProcess())
        {
            case ElecTransfer:
            case Ad1:
            {
                if(Pi[0]->Consider && Ri[0]->Consider)
                {
                    Participators.Add(Pi[0]);
                    Participators.Add(Ri[0]);
                }
                else
                    Remove = true;      // DWDWDW
                break;
            }
            case Des1:
            {
                if(Ri[0]->Consider && Pi[0]->Consider)
                {
                    Participators.Add(Pi[0]);
                    Participators.Add(Ri[0]);
                }
                else
                    Remove = true;      // DWDWDW
                break;
            }
            case React:
            case ReactFromGas:
            {
                if(Ri[0]->Consider && Ri[1]->Consider && Pi[0]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Ri[1]);
                    Participators.Add(Pi[0]);
                    if(Ri[0] == Ri[1] && thisReaction->GetProcess() != ReactFromGas)
                    {
                        thisReaction->StatFactor = 0.5;
                    }
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Diss:
            case DissToGas:
            {
                if(Ri[0]->Consider && Pi[0]->Consider && Pi[1]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Pi[0]);
                    Participators.Add(Pi[1]);
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Ad2:
            {
                if(Pi[0]->Consider && Pi[1]->Consider && Ri[0]->Consider)
                {
                    Participators.Add(Pi[0]);
                    Participators.Add(Pi[1]);
                    Participators.Add(Ri[0]);
                    thisReaction->StatFactor = 0.5;
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Des2:
            {
                if(Ri[0]->Consider && Ri[1]->Consider && Pi[0]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Ri[1]);
                    Participators.Add(Pi[0]);
                    if(Ri[0] == Ri[1])
                        thisReaction->StatFactor = 0.5;
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Des2ToGas:
            {
                if(Ri[0]->Consider && Ri[1]->Consider && Pi[0]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Ri[1]);
                    Participators.Add(Pi[0]);
                    if(Ri[0] == Ri[1])
                        thisReaction->StatFactor = 0.5;
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Ad2ToGas:
            {
                if(Ri[0]->Consider && Pi[1]->Consider && Pi[0]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Pi[1]);
                    Participators.Add(Pi[0]);
                    if(Pi[0] == Pi[1])
                        thisReaction->StatFactor = 0.5;
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Diff:
            {
                if(Ri[0]->Consider)
                {
                    Participators.Add(Ri[0]);
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }
            case Disp:
            case EleyDes:
            case EleyAd:
            case CompAd:
            case CompDes:
            case EleyAdEleyDes:
            {
                if(Ri[0]->Consider &&
                    Ri[1]->Consider && Pi[0]->Consider && Pi[1]->Consider)
                {
                    Participators.Add(Ri[0]);
                    Participators.Add(Ri[1]);
                    Participators.Add(Pi[0]);
                    Participators.Add(Pi[1]);
                    if(thisReaction->GetProcess() != EleyAd &&
                        thisReaction->GetProcess() != EleyAdEleyDes && thisReaction->GetProcess() != CompAd &&
                        (Ri[0] == Ri[1]))
                    {
                        thisReaction->StatFactor = 0.5;
                    }
                }
                else
                {
                    Remove = true;      // DWDWDW
                }
                break;
            }

            default: break;  // DW 8/24/04
        } 
        if(Remove)
        {
            Rxns->Remove();
        }
        else
        {
            ++(*Rxns);
        }
    }
    int RxnCnt = 0;

    Rxns->ResetToFront();
    while(*Rxns)
    {
        Rxns->GetPtr()->_id = RxnCnt;
        RxnCnt++;
        ++(*Rxns);
    }
    NumberOfSpecies = 0;
    Participators.ResetToFront();
    bool Found = false;         // DWDWDW  // DWDWDW

    while(Participators)
    {
        Found = false;          // DWDWDW
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            if(Listing[i] == Participators())
            {
                Found = true;   // DWDWDW
            }
        }
        if(!Found)
        {
            Participators()->_idhook = NumberOfSpecies;
            Listing[NumberOfSpecies] = Participators();
            ++NumberOfSpecies;
        }
        ++Participators;
    }
    if(NumberOfSpecies == 0)
        exit(1);
    Listing[NumberOfSpecies] = theNULLSpecies;
    theNULLSpecies->_idhook = NumberOfSpecies;  // last array member is star species
    assert(NumberOfSpecies < MAX_NUM_SPECIES);
    return NumberOfSpecies;
}

ElementaryRxn *Reactions::PickRandomReaction()
{
    IRandom RanRxn;

    RanRxn.SetInterval(0, Rxns->Size() - 1);
    int RxnNo = RanRxn.Draw();

    Rxns->ResetToFront();
    for(int i = 0; i < RxnNo; ++i)
        ++(*Rxns);
    ElementaryRxn *ChosenReaction = Rxns->GetPtr();

    if(ReactionPossible(ChosenReaction))
        return ChosenReaction;
    else
        return PickRandomReaction();
}

ElementaryRxn *Reactions::ReactionPossible(ElementaryRxn * ChosenReaction)
{
    switch (ChosenReaction->GetProcess())
    {
        case Diff:
        case ElecTransfer:
        {
            if(ChosenReaction->Reactant[0]->NumberOfAtoms > 0)
            {
                return ChosenReaction;
            }
            break;
        }
        case Ad1:
        {
            if (theSimulation->Gas == 1)
            {
              if(ChosenReaction->Reactant[0]->GetPressure() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case Ad2:
        {
            if (theSimulation->Gas == 1)
            {
              if(ChosenReaction->Reactant[0]->GetPressure() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case Des2ToGas:
        {
            if (theSimulation->Gas == 1)
            {
              if(ChosenReaction->Reactant[0]->GetPressure() > 1.e-40 && ChosenReaction->Reactant[1]->NumberOfAtoms > 0)
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40 && ChosenReaction->Reactant[1]->NumberOfAtoms > 0)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case Ad2ToGas:
        {
            if (theSimulation->Gas == 1)
            {
              if(ChosenReaction->Reactant[0]->GetPressure() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case ReactFromGas:
        {
            if (theSimulation->Gas == 1)
            {
              if(ChosenReaction->Reactant[0]->GetPressure() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case Des1:
        case DissToGas:
        case Diss:
        {
            if(ChosenReaction->Reactant[0]->NumberOfAtoms > 0)
            {
                return ChosenReaction;
            }
            break;
        }
        case Des2:
        case Disp:
        case EleyDes:
        case React:
        {
            if(ChosenReaction->Reactant[0]->NumberOfAtoms > 0 &&
                ChosenReaction->Reactant[1]->NumberOfAtoms > 0)
            {
                return ChosenReaction;
            }
            break;
        }
        case CompAd:
        {
            if (theSimulation->Gas == 1)
            {
              if((ChosenReaction->Reactant[0]->GetPressure() > 1.e-40 && ChosenReaction->Reactant[1]->GetPressure() > 1.e-40) )
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if(ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40 && ChosenReaction->Reactant[1]->GetConcentration() > 1.e-40)
              {
                  return ChosenReaction;
              }
            }
            break;
        }
        case CompDes:
        {
            if(ChosenReaction->Reactant[0]->NumberOfAtoms > 0 && ChosenReaction->Reactant[1]->NumberOfAtoms > 0)
            {
                return ChosenReaction;
            }
            break;
        }
        case EleyAd:
        case EleyAdEleyDes:
        {
            if (theSimulation->Gas == 1)
            {
              if((ChosenReaction->Reactant[1]->NumberOfAtoms > 0 && ChosenReaction->Reactant[0]->GetPressure() > 1.e-40) || 
                 (ChosenReaction->Reactant[0]->NumberOfAtoms > 0 && ChosenReaction->Reactant[1]->GetPressure() > 1.e-40) )
              {
                  return ChosenReaction;
              }
            }
            else if (theSimulation->Solvent == 1)
            {
              if((ChosenReaction->Reactant[1]->NumberOfAtoms > 0 && ChosenReaction->Reactant[0]->GetConcentration() > 1.e-40) || 
                  (ChosenReaction->Reactant[0]->NumberOfAtoms > 0 && ChosenReaction->Reactant[1]->GetConcentration() > 1.e-40))
              {
                  return ChosenReaction;
              }
            }
            break;
        }

        default: break;  // DW 8/24/04
    }
    return NULL;
}
