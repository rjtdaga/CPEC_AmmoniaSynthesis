#include "ShortRoutines.h"
#include "Base/Random.h"

Random Num(0.0, 1.0);

double FindRandomNumber()
{
    double r = Num.Draw();
    
    if(r == 0.0 || r == 1.0)
    {
        r = Num.Draw();
    }
    return r;
}

template < class Sw > void Swap(Sw & A, Sw & B)
{
    Sw TempSw = A;

    A = B;
    B = TempSw;
    return;
}
