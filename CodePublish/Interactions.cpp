
#include <cstdio>
#include <cassert>

#include "Species.h"
#include "Component.h"
#include "Interactions.h"
#include "Grid.h"

Interactions* theInteractions;

Interactions::Interactions(fstream &fin, fstream &fout) : Facilitator(fin,fout)
{
	Scale = new double;
	*Scale = 1.0;
	Initialize((CharString)"Interactions", fin);
  int No = theSpecies->Num;
	CharString Next;
	CharString End = "FINISH";
	CharString A, B;
  for (int i = 0; i < 10; i++)
  {
    Pair[0][0][i] = NULL;
  }  
	double val;
  fin >> Next;
  int i = 0;
  while(Next.lower() != End.lower())
  {
    double SurfChar;
    fin >> SurfChar;
    Pair[0][0][i] = SurfChar;
    while (fin >> Next && Next.lower() != "surfacecharge" && Next.lower() != End.lower())
    {
      A = Next;
      fin >> Next;
      B = Next;
      fin >> val;
      Pair[theSpecies->Getid(A)][theSpecies->Getid(B)][i] = val;
      Pair[theSpecies->Getid(B)][theSpecies->Getid(A)][i] = val;
    }
    ++i;
	}
  SetScale(1.0);
    
	if(fin >> Next && Next.lower() == "scale")
  {
        fin >> val;
        SetScale(val);
	}
  theSpecies->theInteractions = this;
}


