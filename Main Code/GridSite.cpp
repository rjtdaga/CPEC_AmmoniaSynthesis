// GridSite.cpp
// ------------


#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdlib>
using namespace std;


#include "GridSite.h"

GridSite *thisSite;

bool operator!=(GridSite & A, GridSite & B);  // DWDWDW
bool operator==(GridSite & A, GridSite & B);  // DWDWDW


GridSite::GridSite()
{
    BondSize = 0;
    SiteType = 0;
    SetType(theNULLSpecies);
    Location.SetXDistance(0.0);
    Location.SetXDistance(0.0);
    Location.SetXDistance(0.0);
    SurroundSize = 0;
    Orientation = 0;
    GridPos[0] = GridPos[1] = GridPos[2] = 0;
    for(Eshort i=0 ; i<PatternDim ; ++i)
    {
        Surround[i] = NULL;
    }
}


GridSite::GridSite(GridSite & a)
{
    BondSize = a.BondSize;
    
    for(Eshort i = 0; i < BondSize; ++i)
    {
        Bonds[i] = a.Bonds[i];
    }
    
    MetalBonds = a.MetalBonds;
    
    for(Eshort ii = 0; ii < a.SurroundSize; ++ii)
    {
        Surround[ii] = a.Surround[ii];
    }
    
    SurroundSize = a.SurroundSize;
    Location = a.Location;
    GridPos[2] = a.GridPos[2];
    GridPos[1] = a.GridPos[1];
    GridPos[0] = a.GridPos[0];
    Type = a.Type;
    Orientation = 0;
}


GridSite::GridSite(const EPosition < double >&place, Component * T)
{

    Location = place;
    SetType(T);  // It is important to set location first!
    
    for(Eshort i = 0; i < PatternDim; ++i)
    {
        Surround[i] = NULL;
    }
    SurroundSize = 0;
    Orientation = 0;
}


GridSite::GridSite(const double &x, const double &y,
    const double &z, Component * A)
{
    Location.SetXDistance(x);
    Location.SetYDistance(y);
    Location.SetZDistance(z);
    SetType(A);   // It is important to set location first!
    
    for(Eshort i = 0; i < PatternDim; ++i)
    {
        Surround[i] = NULL;
    }
    SurroundSize = 0;
    Orientation = 0;
}

void GridSite::GetGridPos(Eshort CurrPos[])
{
    CurrPos[0] = GridPos[0];
    CurrPos[1] = GridPos[1];
    CurrPos[2] = GridPos[2];
}

void GridSite::SetGridPos(Eshort x, Eshort y, Eshort z)
{
    GridPos[0] = x;
    GridPos[1] = y;
    GridPos[2] = z;
}

void GridSite::SetPosition(const EPosition < double >&Place)
{
    Location = Place;
    return;
}

void GridSite::SetPosition(const double &x,
    const double &y, const double &z)
{
    Location.SetXDistance(x);
    Location.SetYDistance(y);
    Location.SetZDistance(z);
    return;
}

void GridSite::SetType(Component * T)
{
    Type = T;
    return;
}

EPosition < double >GridSite::GetPosition() const
{
    return Location;
}

GridSite & GridSite::operator=(GridSite & A)
{
    if(this != &A)
    {
        Type = A.GetType();
        Location = A.GetPosition();
    }
    return *this;
}

GridSite & GridSite::operator=(Component * T)
{
    if(GetType() != T)
    {
        SetType(T);
    }
    return *this;
}

GridSite & GridSite::operator=(const EPosition < double >&Place)
{
    if(&Location != &Place)
    {
        Location = Place;
    }
    return *this;
}

bool operator!=(GridSite & A, GridSite & B) // DWDWDW
{
    if(A.GetPosition() != B.GetPosition())
    {
        return true; // DWDWDW
    }
    else
    {
        return false; // DWDWDW
    }
}

bool operator==(GridSite & A, GridSite & B)  // DWDWDW
{
    if(A.GetPosition() != B.GetPosition())
    {
        return false; // DWDWDW
    }
    else
    {
        return true; // DWDWDW
    }
}

void GridSite::AddBond(GridSite * New)
{
    Bonds[BondSize] = New;
    ++BondSize;
    return;
}


void GridSite::AddMetalBond(GridSite * New)
{
    MetalBonds.Add(New);
    return;
}


double GridSite::FindDistance(GridSite * Away, double GridLength[])
{
    double Delta[3];
    double MyPos[3];
    double NeighPos[3];
    
    MyPos[0] = Location.GetXDistance();
    MyPos[1] = Location.GetYDistance();
    MyPos[2] = Location.GetZDistance();
    NeighPos[0] = Away->Location.GetXDistance();
    NeighPos[1] = Away->Location.GetYDistance();
    NeighPos[2] = Away->Location.GetZDistance();
    Delta[0] = NeighPos[0] - MyPos[0];
    Delta[1] = NeighPos[1] - MyPos[1];
    Delta[2] = NeighPos[2] - MyPos[2];
    double Cmp = 0.5 * GridLength[0];   // x
    
    if(Delta[0] < -Cmp)
    {
        Delta[0] = GridLength[0] + Delta[0];
    }
    else if(Delta[0] > Cmp)
    {
        Delta[0] = -GridLength[0] + Delta[0];
    }
    Cmp = 0.5 * GridLength[1];  // y
    if(Delta[1] < -Cmp)
    {
        Delta[1] = GridLength[1] + Delta[1];
    }
    else if(Delta[1] > Cmp)
    {
        Delta[1] = -GridLength[1] + Delta[1];
    }
    double Distance =
        Delta[0] * Delta[0] + Delta[1] * Delta[1] + Delta[2] * Delta[2];
    Distance = pow(Distance, 0.5);
    return Distance;
}

void GridSite::FindUnitVector(GridSite * Away, double GridLength[],
    double Delta[])
{
    double MyPos[2];
    double NeighPos[2];
    
    MyPos[0] = Location.GetXDistance();
    MyPos[1] = Location.GetYDistance();
    NeighPos[0] = Away->Location.GetXDistance();
    NeighPos[1] = Away->Location.GetYDistance();
    Delta[0] = NeighPos[0] - MyPos[0];
    Delta[1] = NeighPos[1] - MyPos[1];
    for(Eshort k = 0; k < 2; ++k)
    {
        if(abs(Delta[k]) > GridLength[k] * 0.5)
        {
            if(Delta[k] > 0.)
            {
                Delta[k] = -GridLength[k] + Delta[k];
            }
            else
            {
                Delta[k] = GridLength[k] + Delta[k];
            }
        }
    }
}
