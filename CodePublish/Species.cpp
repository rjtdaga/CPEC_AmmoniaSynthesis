#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstdlib>
using namespace std;

#include "Species.h"
#include "Geometry.h"
#include "Component.h"
#include "Base/Estring.h"

Component *theNULLSpecies = NULL;


// DW 8/23/04: I think that the member variable by the same name hides
// this global variable, so I'm commenting it out.
// Actually, commenting it out resulted in a lot of undefined
// references to 'theSpecies'
Species *theSpecies = NULL;


Component **Listing;            // Listing of these species
int NumberOfSpecies;            // Number of SURFACE species


Species::Species(fstream &fin, fstream &fout) : Facilitator(fin, fout)
{
    inp = fin.tellp();
    //outp = fout.tellp();
    Initialize((CharString) "Species", fin);
    //cout<< "\n  Finished Initialize(). inp = " << inp << " outp = " << outp;
    
    int Num = Read(fin);
//    theSpecies.ResetToFront();
//    while (theSpecies)
//    {
//      Component A = theSpecies();
//      cout << A.GetName() << " " << A.Charge << endl;
//      ++(theSpecies);
//    }
//    exit(1);
    if(Num != theSpecies.Size())
    {
        cout << "\nError... Species Redefinition in input file!" << flush;
        exit(1);
    }
    theInteractions = NULL;
    ActiveSpecies = NULL;
    theGeometries = new Geometry(*this, fin, fout);
    assert(theGeometries);    
}

Species::Species(Species *cop, fstream &fin, fstream &fout) :
    Facilitator(fin, fout)
{

    cop->theSpecies.ResetToFront();
    while(cop->theSpecies)
    {
        theSpecies.Add(cop->theSpecies());
        ++cop->theSpecies;
    }
    Num = cop->theSpecies.Size();
    ActiveSpecies = cop->ActiveSpecies;
    theInteractions = cop->theInteractions;
    theGeometries = cop->theGeometries;
}


int Species::Read(fstream &fin)
{
    fin.seekg(inp);             // Seek species section of file
    CharString Next;
    double fNext;
    int iNext;
    CharString End = "END";
    Component New;
    // three Lists needed after completion
    SeqList < CharString > Molecule;
    SeqList < CharString > Strong;
    SeqList < CharString > Weak;
    
    // define the NULL species -always with an id of 0
    Next = "-";
    New.SetName(Next);
    New.SetOrientation((CharString) "None");
    New.Atom = true; // DWDWDW
    New._id = 0;
    New.OnSurface = true; // DWDWDW       // Surface Vacancy
    
    int NumberFound = 1;        // the NULL
    
    theNULLSpecies = set(New);
    
    while(fin >> Next && Next.lower() != End.lower())
    {
        New.OnSurface = false; // DWDWDW
        New.SetName(Next);
        Molecule.Add(Next);
        fin >> fNext;
        New.SetCore(fNext);
        fin >> fNext;
        New.SetVDW(fNext);
        fin >> fNext;
        New.SetAtomizationEnergy(fNext);
        fin >> fNext;
        New.SetMolWt(fNext);
        fin >> Next;
        New.SetClosedShell(PositiveResponse(Next));
        fin >> Next;
        Strong.Add(Next);
        fin >> Next;
        Weak.Add(Next);
        fin >> Next;
        New.SetOrientation(Next);
        fin >> iNext;
        New.SetElectrons(iNext);
        New._id = NumberFound;
        fin >> iNext;
        New.Charge = iNext;
        set(New);
        ++NumberFound;
    }
    Strong.ResetToFront();
    Weak.ResetToFront();
    Molecule.ResetToFront();
    thisComponent = NULL;
    Component *StrongPtr = NULL, *WeakPtr = NULL;
    while(Molecule)
    {
        // DW: removed the following one line 8/24 
        //assert(has(Strong()) && has(Weak()));

        // DW: added the following ...TmpString lines
        // and made some other modifications 8/24
        CharString moleculeTmpString = Molecule();
        CharString weakTmpString = Weak();
        CharString strongTmpString = Strong();

        // DW: changed Molecule() etc args to ...TmpString  8/24
        thisComponent = has(moleculeTmpString);
        WeakPtr = has(weakTmpString);
        StrongPtr = has(strongTmpString);
        
        // DW: added the following assert lines to replace the one that
        // I removed above (8/24)
        assert(WeakPtr != NULL);
        assert(StrongPtr != NULL);
        
        thisComponent->addFragments(StrongPtr, WeakPtr);
        ++Strong;
        ++Weak;
        ++Molecule;
    }

  
    theSpecies.ResetToFront();
    while(theSpecies)
    {
        thisComponent = theSpecies.GetPtr();
        if(isAtom(thisComponent->GetName()))
        {
            thisComponent->Atom = true; // DWDWDW
        }
        else
            thisComponent->Atom = false; // DWDWDW
        ++theSpecies;
    }

    return NumberFound;
}

Component *Species::add(Component A)
{
    CharString MyName = A.GetName();
    Component *Find = has(MyName);

    if(Find == NULL)
    {
        theSpecies.Add(A);
        return has(MyName);
    }
    return Find;
}

Component *Species::add(CharString A)
{
    Component *Find = has(A);

    if(Find == NULL)
    {
        // Make a new component
        Component New;

        New.SetName(A);
        theSpecies.Add(New);
        return has(A);
    }
    return Find;
}

Component *Species::has(CharString & A)
{
    CharString ListName;

    theSpecies.ResetToFront();
    
    // DW: I think this is the cast of SeqList to int, which returns
    // the number of species in the list to the right of the current elt. 
    while(theSpecies) 
    {
        ListName = theSpecies().GetName();
        if(A == ListName)
        {
            Component *theAddress = theSpecies.GetPtr();

            return theAddress;
        }
        ++theSpecies;
    }
    return NULL;
}

Component *Species::get(CharString & A)
{
    Component *Find = has(A);
    
    return Find;
}

Component *Species::set(Component A)
{
    Component *Find = add(A);
    
    *Find = A;                  // overloaded = operator for component 
    return Find;
}

Eshort Species::Getid(CharString A)
{
    Component *B = has(A);

    if(B)
    {
        return B->_id;
    }
    cout << "\nBad getid for species" << A << endl;
    exit(1);
    return true; // DWDWDW
}


bool Species::isAtom(CharString A)
{
    Component *Question = has(A);

    if(Question)
    {
        if(Question->isAtom())
        {
            return true; // DWDWDW
        }
    }
    return false; // DWDWDW
}
