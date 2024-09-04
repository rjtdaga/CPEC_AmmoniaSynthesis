#include "Base/Matrix.h"
#include "Base/Facilitator.h"
class Species;
#ifndef Interaction_h
#define Interaction_h

class Interactions : public Facilitator
{
public:
    
    Interactions(fstream &fin, fstream &fout); 
    
    double GetPairSC(int i)
    {
      return Pair[0][0][i];
    }
    
    double Get(const Component *A, const Component *B, double SC) 
    {
        int got = 0;
        int i = 0;
        while (got == 0)
        {
          if (Pair[0][0][i] == SC)
          {
            got = 1;
            return Pair[A->_id][B->_id][i];
          }
          ++i;
        }   
    }
    
    double Get(const Component A, const Component B, double SC)
    {
        int got = 0;
        int i = 0;
        while (got == 0)
        {
          if (Pair[0][0][i] == SC)
          {
            got = 1;
            return Pair[A._id][B._id][i];
          }
          ++i;
        }
    }
    
    void  SetScale(double val)
    {
        *Scale = val;
        return;
    }
    
    double GetScale() const
    {
        return *Scale;
    }
    
    double *GetScalePointer()
    {
        return Scale;
    }
    
private:    
    double Pair[10][10][10];
    double *Scale;
};

extern Interactions *theInteractions;

#endif
