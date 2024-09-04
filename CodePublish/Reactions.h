
#include "Base/DynArrays.h"
#include "Base/Facilitator.h"

class ElementaryRxn;
class Component;
#ifndef reactions_h
#define reactions_h

#ifndef MAX_NUM_REACTIONS 
#define MAX_NUM_REACTIONS 100
#endif

class Reactions : public Facilitator
{
public:
    
    // constructors
    Reactions(fstream &fin, fstream &fout); 
    int Read(fstream &fin);
    
    void OutputReactions(fstream &fout);
    void PrintReactions();
    
    ElementaryRxn* add(ElementaryRxn New);

    ElementaryRxn* PickARandomReaction();
	
    ElementaryRxn* has(ElementaryRxn A);
    void SetBondDissociationEnergies();
    int FindSurfaceSpecies();
    ElementaryRxn* PickRandomReaction();
    ElementaryRxn* ReactionPossible(ElementaryRxn* ChosenReaction);    
private:
    ElementaryRxn* ActiveReaction;
    int NumberOfReactions;
    SeqList<ElementaryRxn> theReactions;
};



extern SeqList<ElementaryRxn>* Rxns;
extern Reactions* theReactions;

#endif
