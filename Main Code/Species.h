#include <fstream>
#include "Base/DynArrays.h"
#include "Base/Constants.h"
#include "Base/Facilitator.h"

class Geometry;
class Interactions;
#include "Component.h"
#include "Base/Estring.h"
#ifndef species_h
#define species_h

#ifndef MAX_NUM_SPECIES 
#define MAX_NUM_SPECIES 25
#endif

using namespace std;
extern Component* theNULLSpecies;
extern Component **Listing;            // Listing of these species
extern int NumberOfSpecies;            // Number of SURFACE species

class Species : public Facilitator
{

public:

    // constructors
    Species(fstream &fin, fstream &fout);
    Species(Species *S, fstream &fin, fstream &fout);

    inline CharString& GetName(Eshort id)
    {
        theSpecies.ResetToFront();
        for(Eshort i = 0; i < id; ++i )
        {
            ++theSpecies;
        }
        return theSpecies().GetName();
    }
    
    CharString getName(Eshort id)
    {
      Component *A = get(id);
      return A->GetName();
    }
    
    inline Component* get(Eshort id)
    {
        if(id <= 0 || id > theSpecies.Size())
        {
            return theNULLSpecies;
        }
        theSpecies.ResetToFront();
        for(Eshort i=0; i < id; ++i)
        {
            ++theSpecies;
        }
        return theSpecies.GetPtr();
    }
    
    Component* get(CharString &A);
    
    Component* add(CharString A);
    
    Eshort Getid(CharString Name);
    bool isAtom(CharString Name); // DWDWDW
    
    // Inspector for presence of A
    Component* has(CharString &A);
    
    static void colorSelectedCallback (double red, 
        double green, double blue, Component* );

    Interactions *theInteractions;
    Component* ActiveSpecies;
    
    Geometry* theGeometries;	
    
    void OutputSpecies(fstream &fout);
    int Read(fstream &fin);
    
    // Base add function for setting up
    // the species
    Component* add(Component A);
    
    // set function for literal copy of A to theSpecies
    Component* set(Component A);
    int GetSize() {return theSpecies.Size();}
    
    int Num;
    
    SeqList<Component> theSpecies;
};



// DW 08/23/04: This has the same name as the member variable:
//   SeqList<Component> theSpecies;
// Commenting it out causes Reactions.o (and other .o's) to fail to be built.
extern Species* theSpecies;


#endif
