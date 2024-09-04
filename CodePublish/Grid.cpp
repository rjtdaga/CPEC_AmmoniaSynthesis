#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
using namespace std;

#include <math.h>

#include <unistd.h> // needed for sleep()
#include "SystemHeaders.h"
#include "rank.h"
#include "Grid.h"
#include "Geometry.h"
#include "Base/Estring.h"
#include "GridSite.h"
#include "Species.h"
#include "Base/IntRandom.h"
#include "Base/Random.h"
#include "ShortRoutines.h"
#include "Simulation.h"
#include "Base/CommonBlock.h"

Grid *theGrid = NULL;

CharString DistinguishSites(int i, int j);

// Function Prototyping
void GetRealCoordinates(int x, int y, int dx, int dy, int &FinalX,
    int &FinalY);

// CONSTRUCTOR
Grid::Grid(fstream &fin, fstream &fout, int Size) : Facilitator(fin, fout)
{
    // Set the name of coordinates file and set new surface to be zero (true)
    Grid::new_surface = 1;
    Grid::coords_file = coords_file;
    for(int i = 0; i < 7; ++i)
    {
        for(int j = 0; j < PatternDim; ++j)
        {
            SiteDistances[j][i] = -1.;
        }
    }
    AllSite = 1;
    Atop_a = 0;
    Bridge_a = 0;
    Hollow_a = 0;
    ::theGrid = this;
    MetalCoordinates_init = true;
    const char *d1 = "Coordinates";
    strcpy(Coordinates_File, d1);
    strcat(Coordinates_File, rnk);
    strcpy(CoordinatesHan_File, d1);
    strcat(CoordinatesHan_File, "Han");
    strcat(CoordinatesHan_File, rnk);
    strcpy(CoordinatesXYZ_File, d1);
    strcat(CoordinatesXYZ_File, rnk);
    const char *d2 = ".xyz";
    strcat(CoordinatesXYZ_File, d2);
    Cutoff = 2.01;
    NumberOfMetalAtoms = 0;
    OldEDim[1] = OldEDim[0] = 0;
    Alloy_Type = theNULLSpecies;
    Alloy_Composition = 0.;
    AdsorbateLayer = 2;
    Metal = theNULLSpecies;
    SetSurface("100");
    Initialize((CharString) "GRID", fin);
    SubstrateDimension = 0;
    AdsorbateDimension = 0;
    //cout << "Start Read" << endl;
    if(Size == 0)
    {
        Read(fin);
    }
    //cout << "Done Read" << endl;
    Init();
    //cout << "Done Initialization" << endl;
    // Checking whether the surrounding sites are correct or not
//    GridSite *Central = LocateDistinguishSites("Hollow");
//    Eshort a[3];
//    Central->SetType(Listing[0]);
//    Central->GetGridPos(a);
//    cout << "CentralSite: " << FindCoordination(a[0],a[1]) << " " << Central->BondSize << " " << a[0] << " " << a[1] << endl;
//    for (int i = 0; i < Central->SurroundSize; i++)
//    {
//      GridSite *thisSite = Central->Surround[i];
//      Eshort c[3];
//      thisSite->GetGridPos(c);
//      cout << "SurroundSite: " << c[0] << " " << c[1] << endl;
//    }
    
//    exit(1);
    
    RandomNumX.Randomize();
    RandomNumY.Randomize();
    RandomNumSurr.Randomize();
    RandomNumX.SetInterval(0, EDim[0] - 1);
    RandomNumY.SetInterval(0, EDim[1] - 1);
    GridSite *See = PickRandom();
    int GridSiteSurroundSize = See->SurroundSize;
    RandomNumSurr.SetInterval(0, GridSiteSurroundSize);
    MainGrid = true; // DWDWDW
    DefineQuadrants();
}

Grid::Grid(fstream &fin, fstream &fout, int Size, double Cut,
    CharString & Mstr, CharString & Face) : Facilitator(fin, fout)
{

    for(int i = 0; i < 7; ++i)
    {
        for(int j = 0; j < PatternDim; ++j)
        {
            SiteDistances[j][i] = -1.;
        }
    }
    Alloy_Type = theNULLSpecies;
    Alloy_Composition = 0.;
    SubstrateDimension = 0;
    AdsorbateDimension = 0;
    MainGrid = false;  // DWDWDW
    ::theGrid = this;
    Cutoff = Cut;
    NumberOfMetalAtoms = 0;
    EDim[0] = EDim[1] = Size;
    OldEDim[0] = OldEDim[1] = 0;
    AdsorbateLayer = 2;
    Metal = theSpecies->get(Mstr);
    SetSurface(Face);
    Initialize((CharString) "GRID", fin);
    if(Size == 0)
    {
        Read(fin);
    }
    else
    {
        Init();
    }
    RandomNumX.Randomize();
    RandomNumY.Randomize();
    RandomNumSurr.Randomize();
    RandomNumX.SetInterval(0, EDim[0] - 1);
    RandomNumY.SetInterval(0, EDim[1] - 1);
    GridSite *See = PickRandom();
    int GridSiteSurroundSize = See->SurroundSize;

    RandomNumSurr.SetInterval(0, GridSiteSurroundSize);
    DefineQuadrants();
}

double Grid::GetCoverage(Component * PassedComponent)
{
    return PassedComponent->NumberOfAtoms * DeltaCoverage;
}

GridSite *Grid::PickRandom()
{

    int x, y;

    x = RandomNumX.Draw();
    y = RandomNumY.Draw();
    GridSite *SurfaceSite = Surface[AdsorbateLayer][x][y];

    while(SurfaceSite == NULL)
    {
        x = RandomNumX.Draw();
        y = RandomNumY.Draw();
        SurfaceSite = Surface[AdsorbateLayer][x][y];
    }
    
    return SurfaceSite;
}

void Grid::GetSurfaceAtoms(int *Coverage, int &Total)
{

    int x, y;

    thisSite = NULL;
    int i = 0;

    Total = 0;
    assert(Coverage);
    for(i = 0; i < NumberOfSpecies; ++i)
    {
        Coverage[i] = 0;
    }
    for(x = 0; x < EDim[0]; ++x)
    {
        for(y = 0; y < EDim[1]; ++y)
        {
            thisSite = Surface[AdsorbateLayer][x][y];
            if(thisSite != NULL)
            {
                if(thisSite->GetType() != theNULLSpecies)
                {
                    i = 0;
                    while(thisSite->GetType() != Listing[i])
                        ++i;
                    Coverage[i] += 1;
                    Total++;
                }
            }
        }
    }
    return;
}

void Grid::Resize(int Size)
{
    OldEDim[0] = EDim[0];
    OldEDim[1] = EDim[1];
    EDim[1] = Size;
    if(GetSurface() == "111")
        EDim[0] = 2 * Size;
    else if (GetSurface() == "100")
        EDim[0] = Size;
    else if (GetSurface().lower() == "graphene")
        EDim[0] = 2 * Size;
    return;
}

void Grid::ChangeSurface(CharString & Face)
{
    SetSurface(Face);
    Init();
    return;
}

void Grid::Read(fstream &fin)
{
    fin.seekg(inp);
    CharString Next;
    CharString End = "END";
    while(fin >> Next && Next.lower() != End.lower())
    {
        if(Next.lower() == "metal")
        {
            fin >> Next;
            Metal = theSpecies->get(Next);
        }
        else if(Next.lower() == "surface")
        {
            fin >> Next;
            SetSurface(Next);
            SurfaceType = Next;
        }
        else if(Next.lower() == "size")
        {            
            fin >> EDim[0];
            Resize(EDim[0]);          
        }
        else if(Next.lower() == "energies")
        {
            fin >> Next;
            double SurfCharg;
            int j = 0;
            while (Next.lower() != "finish")
            {
              fin >> SurfCharg;
              while (fin >> Next && Next.lower() != "finish" && Next.lower() != "surfacecharge")
              {
                Component *A = theSpecies->get(Next);
                double val;
                for(int i = 1; i < 5; ++i)
                {
                    fin >> val;
                    A->BindingEnergy[i][0][j] = val;
                    if (i == 3)
                    {
                      if (SurfaceType.lower() == "graphene")
                      {
                        A->BindingEnergy[6][0][j] = val;
                      }
                      if (SurfaceType == "100")
                      {
                        A->BindingEnergy[4][0][j] = val;
                      }
                    }
                    
                }
                A->BindingEnergy[5][0][j] = SurfCharg;
              }
              j++;
            }
        }
        else if(Next.lower() == "atop")
        {
          Atop_a = 1;
          AllSite = 0;
        }
        else if(Next.lower() == "bridge")
        {
          Bridge_a = 1;
          AllSite = 0;
        }
        else if(Next.lower() == "hollow")
        {
          Hollow_a = 1;
          AllSite = 0;
        }
        else if(Next.lower() == "cutoff")
        {
            fin >> Cutoff;
        }
        else if(Next.lower() == "alloy")
        {
            fin >> Next;
            double val;

            Alloy_Type = theSpecies->get(Next);
            fin >> val;
            Alloy_Composition = val;
        }
        else if(Next.lower() == "pattern")
        {
            fin >> SubstrateDimension;
            assert(SubstrateDimension <= 12);
            for(int i = 0; i < SubstrateDimension; ++i)
            {
                for(int j = 0; j < SubstrateDimension; ++j)
                {
                    fin >> SubstratePattern[i][j];
                }
            }
        }
        else if(Next.lower() == "adsorbate_pattern")
        {
            fin >> AdsorbateDimension;
            assert(AdsorbateDimension < 24);
            for(int i = 0; i < AdsorbateDimension; ++i)
            {
                for(int j = 0; j < AdsorbateDimension; ++j)
                {
                    fin >> Next;
                    if(Next == "0")
                    {
                        AdsorbatePattern[i][j] = theNULLSpecies->_idhook;
                    }
                    else
                    {
                        AdsorbatePattern[i][j] = theSpecies->get(Next)->_idhook;
                    }
                }
            }
        }
    }
}

int Grid::DePopulate()
{

    int z = AdsorbateLayer, x = 0, y = 0;

    for(x = 0; x < EDim[0]; ++x)
    {
        for(y = 0; y < EDim[1]; ++y)
        {
            if(Surface[z][x][y] != NULL)
            {
                Surface[z][x][y]->SetType(theNULLSpecies);
            }
        }
    }
    return 1;
}

void Grid::SetSurface(CharString Surface)
{

    assert(Surface.lower() == "100" || Surface.lower() == "111" || Surface.lower() == "graphene");
    if(Surface.lower() == "100")
    {
        NumberOfBonds = 9;
        SurfaceType = Surface;
        EDim[0] = EDim[1];
        SurfaceCoordination = 4;        // Sites per unit cell
    }
    else if(Surface.lower() == "111")
    {
        NumberOfBonds = 13;
        SurfaceType = Surface;
        
        // DW Removed 8/24/04 and replaced with the line that does not
        // include a double.
        //EDim[0] = 2. * EDim[1];
        EDim[0] = 2 * EDim[1];
        
        SurfaceCoordination = 6;        // Sites per unit cell
    }
    else if (Surface.lower() == "graphene")
    {
        NumberOfBonds = 7;
        SurfaceType = Surface;
        EDim[0] = 2*EDim[1];
        SurfaceCoordination = 3; 
    }
    return;
}

void Grid::AddMetalAtom(GridSite * PassedSite)
{

    if(SubstrateDimension > 1)
    {
        Eshort CurrentPos[3];

        PassedSite->GetGridPos(CurrentPos);

        CurrentPos[0] = ((CurrentPos[0] / 2) % SubstrateDimension);
        CurrentPos[1] = ((CurrentPos[1] / 2) % SubstrateDimension);
        if(SubstratePattern[CurrentPos[0]][CurrentPos[1]] == 1)
        {
            PassedSite->SetType(Alloy_Type);
        }
        else
        {
            PassedSite->SetType(Metal);
        }
    }
    else
    {
        if(FindRandomNumber() < Alloy_Composition)
        {
            PassedSite->SetType(Alloy_Type);
        }
        else
        {
            PassedSite->SetType(Metal);
        }
    }
    return;
}

void Grid::AddOverLayerAdsorbates()
{
    if(AdsorbateDimension > 1)
    {
        GridSite *PassedSite;
        Component *PassComp;

        for(int i = 0; i < theGrid->EDim[0]; ++i)
        {
            for(int j = 0; j < theGrid->EDim[1]; ++j)
            {
                PassedSite = Surface[AdsorbateLayer][i][j];
                if(PassedSite != NULL)
                {
                    Eshort CurrentPos[3];
                    PassedSite->GetGridPos(CurrentPos);
                    CurrentPos[0] = ((CurrentPos[0]) % AdsorbateDimension);
                    CurrentPos[1] = ((CurrentPos[1]) % AdsorbateDimension);
                    if(AdsorbatePattern[CurrentPos[0]][CurrentPos[1]] != theNULLSpecies->_idhook)
                    {
                        PassComp = Listing[AdsorbatePattern[CurrentPos[0]][CurrentPos[1]]];
                        PassedSite->SetType(PassComp);
                        PassComp->NumberOfAtoms++;
                    }
                }
            }
        }
    }
    
    return;
}

// INITIALIZE GRID
void Grid::SetSiteMatrix()
{
    int x, y, z;

    // Set All to NULL
    for(z = 0; z <= AdsorbateLayer; ++z)
    {
        for(x = 0; x < OldEDim[0]; ++x)
        {
            for(y = 0; y < OldEDim[1]; ++y)
            {
                delete Surface[z][x][y];
            }
        }
    }
    Matrix < GridSite * >CopySurf(EDim[0], EDim[1]);
    for(z = 0; z <= AdsorbateLayer; ++z)
    {
        Surface[z] = CopySurf;
    }
    // Set All to NULL
    for(z = 0; z <= AdsorbateLayer; ++z)
    {
        for(x = 0; x < EDim[0]; ++x)
        {
            for(y = 0; y < EDim[1]; ++y)
            {
                Surface[z][x][y] = NULL;
            }
        }
    }
}

void Grid::PutDownMetalAtoms()
{
    MMdistance = Metal->GetCore() * 2.0;
    //int x, y, z;
    SiteCoordination RotateSite = Atop;
    
    MetalLayer = AdsorbateLayer - 1;
    
    NumberOfMetalAtoms = 0;
    
    if(SurfaceType == "100")
    {
        if(AdsorbateLayer % 2 == 0)
        {
            RotateSite = Hollow;
        }
        
        // SiteArea = MMdistance * MMdistance / 4.0; 
        SiteArea = MMdistance * MMdistance;
        
        GridLength[0] = GridLength[1] = (double (EDim[0]) * MMdistance / 2.0);
        
        for(int z = 0; z <= AdsorbateLayer; ++z)
        {
            for(int x = 0; x < EDim[0]; ++x)
            {
                for(int y = 0; y < EDim[1]; ++y)
                {
                    if((FindCoordination(x, y) == RotateSite) ||
                        z == AdsorbateLayer)
                    {
                        Surface[z][x][y] = new GridSite; // z dimension first
                        assert(Surface[z][x][y] != NULL);
                        
                        GridSite *TempSite = Surface[z][x][y];
                        
                        TempSite->SetPosition(double (x * MMdistance / 2.),
                            double (y * MMdistance / 2.),
                            double (z * 0.7071 * MMdistance));

                        TempSite->SetGridPos(x, y, z);
                        
                        if((FindCoordination(x, y) == RotateSite)
                            && (z != AdsorbateLayer))
                        {
                            // Metal Sites
                            AddMetalAtom(TempSite);
                            
                            if(z == MetalLayer)
                            {
                                ++NumberOfMetalAtoms;
                                
                             }
                        }
                        else
                        {
                            TempSite->SetType(theNULLSpecies);
                        }
                    }
                    else
                    {
                        Surface[z][x][y] = NULL;
                    }
                }
            }
            if(RotateSite == Atop)
            {
                RotateSite = Hollow;
            }
            else
            {
                RotateSite = Atop;
            }
        }
        
    }
    else if(SurfaceType == "111")
    {  
        if(AdsorbateLayer % 2 == 0)
        {
            RotateSite = Hollow_3;
        }
        
        // SiteArea = MMdistance * MMdistance * 0.14433757;  
        SiteArea = MMdistance * MMdistance * 0.8660254; //MM^2*sqrt(3)/2 
        GridLength[0] = double (EDim[0]) * MMdistance / 4.;     // x dir
        GridLength[1] = double (EDim[1]) * MMdistance * 0.4330127; // y dir
        
        bool RowToggle = true;  // DWDWDW  // DWDWDW
        
        assert(EDim[1] % 4 == 0);
        bool TogglePutMetalAtom = true; // DWDWDW  // DWDWDW
        
        for(int z = 0; z <= AdsorbateLayer; ++z)
        {
            for(int x = 0; x < EDim[0]; ++x)
            {
                if(RowToggle)
                {
                    RowToggle = (!RowToggle);
                    TogglePutMetalAtom = (!TogglePutMetalAtom);
                }
                else
                {
                    RowToggle = (!RowToggle);
                }
                
                for(int y = 0; y < EDim[1]; ++y)
                {
                    //  TogglePutMetalAtom = (!TogglePutMetalAtom);
                    SiteCoordination ThisCoord = FindCoordination(x, y);
                    
                    if(ThisCoord == Atop && (ThisCoord == RotateSite ||
                            z == AdsorbateLayer))
                    {
                        // Atop Sites
                        
                        Surface[z][x][y] = new GridSite;
                        assert(Surface[z][x][y] != NULL);

                        GridSite *Pic = Surface[z][x][y];
                        Pic->SetType(theNULLSpecies);
                        Pic->SetGridPos(x, y, z);
                        EPosition <double> c = Pic->GetPosition();
                        Pic->SetPosition(double (x * MMdistance / 4.),
                            double (y * MMdistance * 0.4330127),
                            double (z * 0.8167 * MMdistance));
                        if(z != AdsorbateLayer)
                        {
                            AddMetalAtom(Pic);
                            if(z == MetalLayer){
                                ++NumberOfMetalAtoms;
                                
                                }
                        }
                        
                        
                    }
                    
                    else if(ThisCoord == Hollow_3 &&
                        (ThisCoord == RotateSite || z == AdsorbateLayer))
                    {
                    
                        if(TogglePutMetalAtom || z == AdsorbateLayer)
                        {
                            Surface[z][x][y] = new GridSite;
                            assert(Surface[z][x][y] != NULL);
                            
                            GridSite *Pic = Surface[z][x][y];
                            Pic->SetType(theNULLSpecies);
                            Pic->SetGridPos(x, y, z);
                            if(z != AdsorbateLayer)
                            {
                                AddMetalAtom(Pic);
                            }
                            else
                            {
                                Pic->SetType(theNULLSpecies);
                            }
                            if((x + y) % 4 == 0)
                            {
                                Pic->SetPosition(double (x * MMdistance / 4.0),
                                    double ((y+1) * MMdistance * 0.4330127 -
                                        0.28867513 * MMdistance),
                                    double (z * 0.8167 * MMdistance));
                            }
                            else
                            {
                                Pic->SetPosition(double (x * MMdistance / 4.0),
                                    double ((y+1) * MMdistance * 0.4330127 -
                                        0.57735027 * MMdistance),
                                    double (z * 0.8167 * MMdistance));
                            }
                            TogglePutMetalAtom = (!TogglePutMetalAtom);
                        }
                        else
                        {
                            TogglePutMetalAtom = (!TogglePutMetalAtom);
                        }
                    }
                    else if(ThisCoord == Bridge)
                    { 
                        if(z == AdsorbateLayer)
                        {
                            // For Overlayer Only
                            Surface[z][x][y] = new GridSite;
                            assert(Surface[z][x][y] != NULL);
                            
                            GridSite *Pic = Surface[z][x][y];
                            Pic->SetType(theNULLSpecies);
                            Pic->SetGridPos(x, y, z);
                            Pic->SetPosition(double (x * MMdistance / 4.),
                                double (y * MMdistance * 0.4330127),
                                double (z * 0.8167 * MMdistance));
                        }
                    }
                }
            }
            
            if(RotateSite == Atop)
            {
                RotateSite = Hollow_3;
            }
            else
            {
                RotateSite = Atop;
            }
            
        }
        
    }
    else if(SurfaceType.lower() == "graphene")
    {  
        SiteArea = MMdistance * MMdistance * 1.3; //MM^2*sqrt(3)/2 
        GridLength[0] = double (EDim[0]) * 3 * MMdistance / 8.;     // x dir
        GridLength[1] = double (EDim[1]) * MMdistance * 0.4330127 * 2; // y dir
        assert(EDim[1] % 4 == 0);
        for(int z = MetalLayer; z <= AdsorbateLayer; ++z)
        {
            for(int x = 0; x < EDim[0]; ++x)
            {
                for(int y = 0; y < EDim[1]; ++y)
                {
                    SiteCoordination ThisCoord = FindCoordination(x, y);
                    if(ThisCoord == Atop)
                    {
                        Surface[z][x][y] = new GridSite;
                        assert(Surface[z][x][y] != NULL);
                        GridSite *Pic = Surface[z][x][y];
                        Pic->SetType(theNULLSpecies);
                        Pic->SetGridPos(x, y, z);
                        if (x%8 == 0)
                        {
                          Pic->SetPosition(double (3*(x/8)*MMdistance),
                                double ((y-2)/4*MMdistance*1.732 + MMdistance*0.866),
                                double (z*3.35));
                        }
                        else if (x%8 == 2)
                        {
                          Pic->SetPosition(double (MMdistance/2 + 3*(x-2)/8*MMdistance),
                                double (y/4*MMdistance*1.732),
                                double (z*3.35));
                        }
                        else if (x%8 == 4)
                        {
                          Pic->SetPosition(double (3*MMdistance/2 + 3*(x-4)/8*MMdistance),
                                double (y/4*MMdistance*1.732),
                                double (z*3.35));
                        }
                        else if (x%8 == 6)
                        {
                          Pic->SetPosition(double (MMdistance*2 + 3*(x-6)/8*MMdistance),
                                double ((y-2)/4*MMdistance*1.732 + MMdistance*0.866),
                                double (z*3.35));
                        }
                        if(z != AdsorbateLayer)
                        {
                            AddMetalAtom(Pic);
                            if(z == MetalLayer)
                            {
                                ++NumberOfMetalAtoms;
                            }
                        }
                    }
                    else if(ThisCoord == Hollow_6 && z == AdsorbateLayer)
                    {
                          Surface[z][x][y] = new GridSite;
                          assert(Surface[z][x][y] != NULL);
                          GridSite *Pic = Surface[z][x][y];
                          Pic->SetType(theNULLSpecies);
                          Pic->SetGridPos(x, y, z);
                          if((x-3) % 8 == 0)
                          {
                              Pic->SetPosition(double (3*(x-3)/8*MMdistance + MMdistance),
                                  double ((y-2)/4*MMdistance*1.732 + MMdistance*0.866),
                                  double (z*3.35));
                          }
                          else if ((x+1)%8 == 0)
                          {
                              Pic->SetPosition(double (3*(x-7)/8*MMdistance + 2.5*MMdistance),
                                  double (y/4*MMdistance*1.732),
                                  double (z*3.35));
                          }
                    }
                    else if(ThisCoord == Bridge && z == AdsorbateLayer)
                    { 
                          // For Overlayer Only
                          Surface[z][x][y] = new GridSite;
                          assert(Surface[z][x][y] != NULL);
                          GridSite *Pic = Surface[z][x][y];
                          Pic->SetType(theNULLSpecies);
                          Pic->SetGridPos(x, y, z);
                          if((x-1) % 4 == 0)
                          {
                              Pic->SetPosition(double (3*(x-1)/4*MMdistance/2 + MMdistance/4),
                                  double ((y-1)/4*MMdistance*1.732 + MMdistance*0.866/2),
                                  double (z*3.35));
                          }
                          else if ((x-3)%8 == 0)
                          {
                              Pic->SetPosition(double (3*(x-3)/8*MMdistance + MMdistance),
                                  double (y/4*MMdistance*1.732),
                                  double (z*3.35));
                          }
                          else if ((x-7)%8 == 0)
                          {
                              Pic->SetPosition(double (3*(x-7)/8*MMdistance + 2.5*MMdistance),
                                  double ((y-2)/4*MMdistance*1.732 + MMdistance*0.866),
                                  double (z*3.35));
                          }
                    }
                }
            }
        }
    }
}

void Grid::WriteMetalCoordinates(GridSite *Pic)
{   
    std::ofstream fe(Coordinates_File, ios::app);
    EPosition <double> c = Pic->GetPosition();
    Component *a = Metal;
    fe << a->GetName() << " " << c.GetXDistance() << " " << c.GetYDistance() << " " << c.GetZDistance();
    fe << endl;
    fe.close();
    return;
}

void Grid::WriteCoordinates(int Frame)
{   
    std::ofstream fe(Coordinates_File, ios::app);
    fe << Frame;
    fe << endl;
    int x = 0;
    int y = 0;
    while (x<EDim[0])
    {
          y = 0;
          while (y<EDim[1])
          {
              if (Surface[AdsorbateLayer][x][y]!=NULL)
              {
                GridSite *G_Site = Surface[AdsorbateLayer][x][y];
                EPosition <double> c = G_Site->GetPosition();
                Component *a = G_Site->GetType();
                theSpecies->theGeometries->theCoordinates.ResetToFront();
                int coord_cntr = 0;
                Coordinates* thisCoord;
                while (theSpecies->theGeometries->theCoordinates)
                {
                  thisCoord = theSpecies->theGeometries->theCoordinates.GetPtr();
                  if (a->GetName() == thisCoord->SpecName)
                  {
                    for (int i = 0; i < thisCoord->NumAtoms; i++)
                    {
                      EPosition <double> Coords = thisCoord->Position[i];
                      double AtomCoords[3];
                      AtomCoords[0] = c.GetXDistance() + Coords.GetXDistance();
                      AtomCoords[1] = c.GetYDistance() + Coords.GetYDistance();
                      AtomCoords[2] = c.GetZDistance() + Coords.GetZDistance();
                      fe << thisCoord->AtomName[i] << " " << AtomCoords[0] << " " << AtomCoords[1] << " " << AtomCoords[2];
                      fe << endl;
                    }
                    break;
                  }
                  ++coord_cntr;
                  ++(theSpecies->theGeometries->theCoordinates);
                }
//                if (a->GetName() != "-"){
//                    fe << a->GetName() << " " << c.GetXDistance() << " " << c.GetYDistance() << " " << c.GetZDistance();
//                    fe << endl;
//                }
              }
              ++y;
          }
          ++x;
    }
    fe.close();
    return;
}

void Grid::WriteCoordinatesXYZ(int Frame)
{   
    std::ofstream fe(CoordinatesXYZ_File, ios::app);
    int TotalNumberOfAtoms = GetNumMetalAtoms();
    for(int i = 0; i < NumberOfSpecies; ++i)
    {
      int NumAtom = 0;
      Coordinates* thisCoord;
      theSpecies->theGeometries->theCoordinates.ResetToFront();
      while (theSpecies->theGeometries->theCoordinates)
      {
        thisCoord = theSpecies->theGeometries->theCoordinates.GetPtr();
        if (Listing[i]->GetName() == thisCoord->SpecName)
        {
          NumAtom = thisCoord->NumAtoms;
        }
        ++(theSpecies->theGeometries->theCoordinates);
      }
      TotalNumberOfAtoms += Listing[i]->NumberOfAtoms*NumAtom;
    }
    int x = 0;
    int y = 0;
    fe << TotalNumberOfAtoms << endl;
    fe << Frame << endl;
    int m = 0;
    for (int i = 0; i < EDim[0]; i++)
    {
      for (int j = 0; j < EDim[1]; j++)
      {
        if (Surface[AdsorbateLayer][i][j] != NULL && FindCoordination(i,j) == Atop)
        {
          EPosition <double> c = Surface[AdsorbateLayer-1][i][j]->GetPosition();
          CharString a  = Metal->GetName();
          if (a.lower() == "gr")
            a = "C";
          fe << a << " " << c.GetXDistance() << " " << c.GetYDistance() << " " << c.GetZDistance() << endl;
          m++;
        }
      }
    }
//    while (x<EDim[0]){
//          y = 0;
//          while (y<EDim[1]){
//              if (Surface[AdsorbateLayer][x][y]!=NULL && m < TotalNumberOfAtoms){
//                GridSite *G_Site = Surface[AdsorbateLayer][x][y];
//                EPosition <double> c = G_Site->GetPosition();
//                Component *a = G_Site->GetType();
//                theSpecies->theGeometries->theCoordinates.ResetToFront();
//                int coord_cntr = 0;
//                Coordinates* thisCoord;
//                while (theSpecies->theGeometries->theCoordinates)
//                {
//                  thisCoord = theSpecies->theGeometries->theCoordinates.GetPtr();
//                  if (a->GetName() == thisCoord->SpecName)
//                  {
//                    for (int i = 0; i < thisCoord->NumAtoms; i++)
//                    {
//                      EPosition <double> Coords = thisCoord->Position[i];
//                      double AtomCoords[3];
//                      AtomCoords[0] = c.GetXDistance() + Coords.GetXDistance();
//                      AtomCoords[1] = c.GetYDistance() + Coords.GetYDistance();
//                      AtomCoords[2] = c.GetZDistance() + Coords.GetZDistance();
//                      fe << thisCoord->AtomName[i] << " " << AtomCoords[0] << " " << AtomCoords[1] << " " << AtomCoords[2];
//                      fe << endl;
//                      ++m;
//                    }
//                    break;
//                  }
//                  ++coord_cntr;
//                  ++(theSpecies->theGeometries->theCoordinates);
//                }
//                }
//                ++y;
//          }
//          ++x;
//    }
    while (x<EDim[0])
    {
          y = 0;
          while (y<EDim[1])
          {
//              if (Surface[AdsorbateLayer][x][y]!=NULL)
//                cout << Surface[AdsorbateLayer][x][y]->GetType()->GetName() << endl;
              if (Surface[AdsorbateLayer][x][y]!=NULL && Surface[AdsorbateLayer][x][y]->GetType() != theNULLSpecies)
              {
                GridSite *Site = Surface[AdsorbateLayer][x][y];
                //cout << Surface[AdsorbateLayer][x][y]->GetType()->GetName() << " " << m << " " << Site->NumAtoms << endl;
                for (int i = 0; i < Site->GetType()->NumAtoms; ++i)
                {
                  fe << Site->GetType()->AtomName[i] << " " << Site->Ads_Pos[i][0] << " " << Site->Ads_Pos[i][1] << " " << Site->Ads_Pos[i][2] << endl;
                  //cout << Site->GetType()->AtomName[i] << " " << Site->Ads_Pos[i][0] << " " << Site->Ads_Pos[i][1] << " " << Site->Ads_Pos[i][2] << endl;
                  ++m;
                }
                
              }
              ++y;
          }
          ++x;
    }
    fe.close();
    //cout << m << " " << TotalNumberOfAtoms << endl;
    assert(m == TotalNumberOfAtoms);
    return;
}

void Grid::WriteCoordinatesHan()
{   
    std::ofstream fe(CoordinatesHan_File, ios::app);
    int x = 0;
    int y = 0;
    int m = 0;
    for (int i = 0; i < EDim[0]; i++)
    {
      for (int j = 0; j < EDim[1]; j++)
      {
        if (Surface[AdsorbateLayer][i][j] != NULL && FindCoordination(i,j) == Atop)
        {
          EPosition <double> c = Surface[AdsorbateLayer-1][i][j]->GetPosition();
          fe << Metal->GetName() << endl << "0" << endl << c.GetXDistance() << endl << c.GetYDistance() << endl << c.GetZDistance() << endl << endl;
          m++;
        }
      }
    }
    while (x<EDim[0]){
          y = 0;
          while (y<EDim[1]){
              if (Surface[AdsorbateLayer][x][y]!=NULL){
                GridSite *G_Site = Surface[AdsorbateLayer][x][y];
                EPosition <double> c = G_Site->GetPosition();
                Component *a = G_Site->GetType();
                if (a->GetName() != "-"){
                    fe << a->GetName() << endl << "0" << endl << c.GetXDistance() << endl << c.GetYDistance() << endl << c.GetZDistance();
                    fe << endl << endl;
                    m++;
                }
                }
                ++y;
          }
          ++x;
    }
    fe << "END" << endl;
    fe.close();
    return;
}

// INITIALIZE GRID
void Grid::Init()
{   SetSiteMatrix();
    if(EDim[0] > 0)
    {
        PutDownMetalAtoms();
        ConnectMetalAndOverlayer();
    }
    DeltaCoverage = 1.0 / (double)GetNumMetalAtoms();
    return;
}


void Grid::ConnectMetalAndOverlayer()
{
    int dx, dy;
    int x, y;
    typedef struct
    {
        int x, y;
    } Patt;
    
    // a 100 surface
    
    Patt Pattern[13]; // Surrounding sites
    Patt MetalPattern[9]; // Surrounding metal sites
    
    // Variables used to find true x and y surrounding coordinates
    // difficult due to periodic boundary conditions
    
    int MetalIndex;
    int SiteIndex = NumberOfBonds; // = 7
    int FinalX, FinalY;
    int local;
    
    bool DefinedDist[5] = { false, false, false, false, false }; // DWDWDW x6
    
        
    if(SurfaceType == "100")
    {
        
        // Connecting This Surface
        
        MetalPattern[0].x = 0;
        MetalPattern[0].y = 2;
        MetalPattern[1].x = 2;
        MetalPattern[1].y = 2;
        MetalPattern[2].x = 2;
        MetalPattern[2].y = 0;
        MetalPattern[3].x = 2;
        MetalPattern[3].y = -2;
        MetalPattern[4].x = 0;
        MetalPattern[4].y = -2;
        MetalPattern[5].x = -2;
        MetalPattern[5].y = -2;
        MetalPattern[6].x = -2;
        MetalPattern[6].y = 0;
        MetalPattern[7].x = -2;
        MetalPattern[7].y = 2;
        MetalIndex = 8;
        
        Pattern[0].x = 0;
        Pattern[0].y = 0;
        Pattern[1].x = 0;
        Pattern[1].y = 1;       // create pointers in a clockwise direction 
        Pattern[2].x = 1;
        Pattern[2].y = 1;       // starting at the center, then to N, NE, ...
        Pattern[3].x = 1;
        Pattern[3].y = 0;       // for a 100 surface
        Pattern[4].x = 1;
        Pattern[4].y = -1;
        Pattern[5].x = 0;
        Pattern[5].y = -1;
        Pattern[6].x = -1;
        Pattern[6].y = -1;
        Pattern[7].x = -1;
        Pattern[7].y = 0;
        Pattern[8].x = -1;
        Pattern[8].y = 1;
        
        Patt SNormalPattern[400];
        Patt SNormalPatternRotate90[400];
        Patt In;
        
        // Identify the STANDARD way of searching for surrounding sites
        
        int SurrSize = 4;
        int pcnt = 0;
        
        for(int i = -SurrSize; i <= SurrSize; ++i)
        {
            for(int j = -SurrSize; j <= SurrSize; ++j)
            {
                In.x = i;
                In.y = j;
                SNormalPattern[pcnt] = In;
                ++pcnt;
            }
        }

        // Search grid is rotated by 90 degrees, used for bridging site

        pcnt = 0;
        int scnt = 0;
        int TotalL = (2 * SurrSize + 1);

        for(int ii = -SurrSize; ii <= SurrSize; ++ii)
        {
            for(int j = -SurrSize; j <= SurrSize; ++j)
            {
                scnt = (SurrSize - j) * TotalL + SurrSize + ii;
                SNormalPatternRotate90[pcnt] = SNormalPattern[scnt];
                ++pcnt;
            }
        }

        for(x = 0; x < EDim[0]; ++x)
        {
            for(y = 0; y < EDim[1]; ++y)
            {
                for(local = 0; local < MetalIndex; ++local)
                {
                    dx = MetalPattern[local].x;
                    dy = MetalPattern[local].y;
                    GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                    if(FindCoordination(x, y) == Atop)
                    {           // Connecting M-M bonds on surface

                        Surface[MetalLayer][x][y]->
                            AddMetalBond(Surface[MetalLayer][FinalX]
                                [FinalY]);
                    }
                }

                for(local = 0; local < SiteIndex; ++local)
                {
                    dx = Pattern[local].x;
                    dy = Pattern[local].y;
                    GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);

                    if(FindCoordination(x, y) == Atop)
                    {           // Connecting Metal with
                        // Atop Sites with Overlayer
                        Surface[MetalLayer][x][y]->
                            AddBond(Surface[AdsorbateLayer][FinalX]
                                [FinalY]);
                    }

                    if(FindCoordination(FinalX, FinalY) == Atop)
                    {           // Connecting Overlayer sites
                        if(Surface[AdsorbateLayer][x][y] != NULL)
                        {       // With metal layer
                            Surface[AdsorbateLayer][x][y]->
                                AddBond(Surface[MetalLayer][FinalX]
                                    [FinalY]);
                        }
                    }
                    if(FindCoordination(FinalX, FinalY) == Atop)
                    {
                        if(Surface[AdsorbateLayer][FinalX][FinalY]->
                            BondSize > 1)
                        {
                            cout << "\nError in Atop Definition\n" << flush;
                        }
                    }
                    else if(FindCoordination(FinalX, FinalY) == Bridge)
                    {
                        if(Surface[AdsorbateLayer][FinalX][FinalY]->
                            BondSize > 2)
                        {
                            cout << "\nError in Bridge Definition\n" <<
                                flush;
                        }
                    }
                    else if(FindCoordination(FinalX, FinalY) == Hollow)
                    {
                        if(Surface[AdsorbateLayer][FinalX][FinalY]->
                            BondSize > 4)
                        {
                            cout << "\nError in Hollow Definition\n" <<
                                flush;
                        }
                    }
                }

                GridSite *MySite = Surface[AdsorbateLayer][x][y];

                for(local = 0; local < pcnt; ++local)
                {

                    if(x % 2 == 0 && y % 2 != 0)
                    {           // rotated bridge site
                        dx = SNormalPatternRotate90[local].x;
                        dy = SNormalPatternRotate90[local].y;
                    }
                    else
                    {
                        dx = SNormalPattern[local].x;
                        dy = SNormalPattern[local].y;
                    }
                    GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);

                    GridSite *SurrSite =
                        Surface[AdsorbateLayer][FinalX][FinalY];

                    if(SurrSite != NULL)
                    {
                        if(MySite != NULL)
                        {

                            double delx =
                                double (dx) * theGrid->MMdistance / 2.0;
                            double dely =
                                double (dy) * theGrid->MMdistance / 2.0;

                            double Distance =
                                pow((pow(delx, 2.0) + pow(dely, 2.0)), 0.5);
                            
                            if(Distance <= Cutoff * MMdistance)
                            {
                                MySite->Surround[MySite->SurroundSize] =
                                    SurrSite;
                                assert(MySite->SurroundSize < PatternDim);
                                if(!DefinedDist[FindCoordination(x, y)])
                                {
                                    int ival = FindCoordination(x, y);
                                    
                                    // Saving Distances
                                    SiteDistances[MySite->SurroundSize][ival]=
                                        (float)Distance;
                                }
                                ++MySite->SurroundSize;
                            }
                        }
                    }
                }
                DefinedDist[FindCoordination(x, y)] = true; // DWDWDW
            }
        }
    }
    else if(SurfaceType == "111")
    {
        Patt SitePattern[6][300]; // Matrix containing the search pattern of various site types
        Patt MetalMetal[6]; 
        Patt MetalAdsor[13];
        
        MetalMetal[0].x = 4;
        MetalMetal[0].y = 0;
        MetalMetal[1].x = -2;
        MetalMetal[1].y = -2;
        MetalMetal[2].x = -2;      // Connecting metals and metals for metal type 1
        MetalMetal[2].y = 2;
        MetalMetal[3].x = 2;
        MetalMetal[3].y = -2;
        MetalMetal[4].x = 2;
        MetalMetal[4].y = 2;
        MetalMetal[5].x = -4;      
        MetalMetal[5].y = 0;
        MetalIndex = 6;
        
        MetalAdsor[0].x = 0;
        MetalAdsor[0].y = 0;
        MetalAdsor[1].x = 0;
        MetalAdsor[1].y = 1;       
        MetalAdsor[2].x = 0;
        MetalAdsor[2].y = -1;       // Connecting metals and adsorbates for metal type 1
        MetalAdsor[3].x = 1;
        MetalAdsor[3].y = 1;
        MetalAdsor[4].x = 1;
        MetalAdsor[4].y = -1;
        MetalAdsor[5].x = -1;
        MetalAdsor[5].y = -1;
        MetalAdsor[6].x = -1;
        MetalAdsor[6].y = 1;
        MetalAdsor[7].x = 2;
        MetalAdsor[7].y = 0;
        MetalAdsor[8].x = -2;
        MetalAdsor[8].y = 0;
        MetalAdsor[9].x = 2;
        MetalAdsor[9].y = 1;
        MetalAdsor[10].x = 2;
        MetalAdsor[10].y = -1;
        MetalAdsor[11].x = -2;
        MetalAdsor[11].y = -1;
        MetalAdsor[12].x = -2;
        MetalAdsor[12].y = 1;
        Patt In;
        
        int pcnt; // Number of surrounding sites taken into account
        int SurrSize = 8;
        CharString SiteNames[6];
        SiteNames[0] = "Atop_1";
        SiteNames[1] = "Bridge_1";
        SiteNames[2] = "Bridge_2";    // Types of sites
        SiteNames[3] = "Bridge_3";
        SiteNames[4] = "Hollow_1";
        SiteNames[5] = "Hollow_2";
        
        int Index[6]; // Number of Surrounding sites for each type of site
        for (int i = 0; i < 6; i++)
        {
          Index[i] = 0;
        }
        for (int k = 0; k < 6; k++)
        {
          if (Atop_a == 0 && (k == 0))
            continue;
          if (Bridge_a == 0 && (k == 1 || k == 2 || k == 3))
            continue;
          if (Hollow_a == 0 && (k == 5 || k == 4))
            continue;
            
          GridSite *Central = LocateDistinguishSites(SiteNames[k]);
          Eshort c[3];
          Central->GetGridPos(c);
          pcnt = 0;
          
          for (int i = SurrSize; i >= -SurrSize; i--)
          {
              for (int j = SurrSize; j >= -SurrSize; j--)
              {
                  int FinalX, FinalY;
                  GetRealCoordinates(c[0], c[1], i, j, FinalX, FinalY);
                  if (Surface[c[2]][FinalX][FinalY] == NULL)
                    continue;
                  GridSite *thisSite = Surface[c[2]][FinalX][FinalY];
                  if (thisSite != NULL && Central->FindDistance(thisSite, GridLength) < Cutoff*MMdistance)
                  {
                    SiteDistances[pcnt][k] = Central->FindDistance(thisSite, GridLength);
                    SitePattern[k][pcnt].x = i;
                    SitePattern[k][pcnt].y = j;
                    ++pcnt;
                    ++Index[k];
                  }
              }
          }
        }
        //exit(1);
        //cout << EDim[0] << " " << EDim[1] << endl;
        //theSimulation->Univ = 1;
        for(int x = 0; x < EDim[0]; ++x)
        {
            for(int y = 0; y < EDim[1]; ++y)
            {
                // Connecting Metal-Metal Bonds
                if(DistinguishSites(x, y) == "Atop_1")
                {
                  //cout << "yes" << endl;
                  for(local = 0; local < MetalIndex; ++local)
                  {
                      dx = MetalMetal[local].x;
                      dy = MetalMetal[local].y;
                      GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                      Surface[MetalLayer][x][y]->
                        AddMetalBond(Surface[MetalLayer][FinalX][FinalY]);
                      
                  }
                }
                // Connecting Metal-Adsorbate Bonds
                if(DistinguishSites(x, y) == "Atop_1")
                {
                    for(local = 0; local < SiteIndex; ++local)
                    {   
                        dx = MetalAdsor[local].x;
                        dy = MetalAdsor[local].y;
                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                        if (Surface[AdsorbateLayer][FinalX][FinalY] == NULL)
                          continue;
                        Surface[MetalLayer][x][y]->
                            AddBond(Surface[AdsorbateLayer][FinalX][FinalY]);
                        Surface[AdsorbateLayer][FinalX][FinalY]->
                            AddBond(Surface[MetalLayer][x][y]);
                        
                        if(FindCoordination(FinalX, FinalY) == Atop)
                        {
                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                BondSize > 1)
                            {
                                cout << "\nError in Atop Definition\n" << flush;
                                exit(1);
                            }
                        }
                        else if(FindCoordination(FinalX, FinalY) == Bridge)
                        {
                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                BondSize > 2)
                            {
                                cout << "\nError in Bridge Definition\n" << flush;
                                exit(1);
                            }
                        }
                        else if(FindCoordination(FinalX, FinalY) == Hollow_3)
                        {
                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                BondSize > 3)
                            {
                                cout << "\nError in Hollow Definition\n" << flush;
                                exit(1);
                            }
                        }
                    }
                }
                
                if (Surface[AdsorbateLayer][x][y] == NULL)
                  continue;
                GridSite *MySite = Surface[AdsorbateLayer][x][y];
                if(MySite != NULL && FindCoordination(x, y) != NONE)
                {
                    Patt *CurrentPattern = NULL;
                    CharString WhichSite = DistinguishSites(x, y);
                    int indexval;
                    if(WhichSite == "Atop_1")
                    {
                        CurrentPattern = SitePattern[0];
                        indexval = 0;
                        MySite->SiteType = 0;
                    }
                    else if(WhichSite == "Bridge_1")
                    {
                        CurrentPattern = SitePattern[1];
                        indexval = 2;
                        MySite->SiteType = 0;
                    }
                    if(WhichSite == "Bridge_2")
                    {
                        CurrentPattern = SitePattern[2];
                        indexval = 3;
                        MySite->SiteType = 1;
                    }
                    else if(WhichSite == "Bridge_3")
                    {
                        CurrentPattern = SitePattern[3];
                        indexval = 4;
                        MySite->SiteType = 2;
                    }
                    else if(WhichSite == "Hollow_1")
                    {
                        CurrentPattern = SitePattern[4];
                        indexval = 5;
                        MySite->SiteType = 0;
                    }
                    else if(WhichSite == "Hollow_2")
                    {
                        CurrentPattern = SitePattern[5];
                        indexval = 5;
                        MySite->SiteType = 1;
                    }
                    for(local = 0; local < Index[indexval]; ++local)
                    {
                        int dx, dy;
                        dx = CurrentPattern[local].x;
                        dy = CurrentPattern[local].y;
                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                        GridSite *SurrSite = Surface[AdsorbateLayer][FinalX][FinalY];
                        assert(SurrSite != NULL);
                        MySite->Surround[MySite->SurroundSize] = SurrSite;
                        ++MySite->SurroundSize;
                        assert(MySite->SurroundSize < PatternDim);
                    }
                }
            }
        }
    }
    
    
//    Patt BridgeTypeOne[300];
//        Patt BridgeTypeTwo[300];
//        Patt SNormalPattern[300];
//        Patt SNormalPatternRotate180[300];
//        
//        MetalPattern[0].x = 2;
//        MetalPattern[0].y = 2;
//        MetalPattern[1].x = 4;
//        MetalPattern[1].y = 0;
//        MetalPattern[2].x = 2;
//        MetalPattern[2].y = -2;
//        MetalPattern[3].x = -2;
//        MetalPattern[3].y = -2;
//        MetalPattern[4].x = -4;
//        MetalPattern[4].y = 0;
//        MetalPattern[5].x = -2;
//        MetalPattern[5].y = 2;
//        MetalIndex = 6;
//        
//        Pattern[0].x = 0;
//        Pattern[0].y = 0;
//        Pattern[1].x = 0;
//        Pattern[1].y = 1;       // create pointers in a clockwise direction 
//        Pattern[2].x = 1;
//        Pattern[2].y = 1;       // starting at the center, then to N, NE, ...
//        Pattern[3].x = 2;
//        Pattern[3].y = 1;
//        Pattern[4].x = 2;
//        Pattern[4].y = 0;
//        Pattern[5].x = 2;
//        Pattern[5].y = -1;
//        Pattern[6].x = 1;
//        Pattern[6].y = -1;
//        Pattern[7].x = 0;
//        Pattern[7].y = -1;
//        Pattern[8].x = -1;
//        Pattern[8].y = -1;
//        Pattern[9].x = -2;
//        Pattern[9].y = -1;
//        Pattern[10].x = -2;
//        Pattern[10].y = 0;
//        Pattern[11].x = -2;
//        Pattern[11].y = 1;
//        Pattern[12].x = -1;
//        Pattern[12].y = 1;
//        Patt In;
//        
//        // Identify the STANDARD way of searching for surrounding sites
//        
//        int SurrSize = 8;
//        int pcnt = 0;
//        
//        int i;
//        int j;
//        for(j = SurrSize; j >= -SurrSize; --j)
//        {
//            int i;
//            for(i = -SurrSize; i <= SurrSize; ++i)
//            {
//                In.x = i;
//                In.y = j;
//                SNormalPattern[pcnt] = In;
//                ++pcnt;
//            }
//        }
//        pcnt = 0;
//        for(j = -SurrSize; j <= SurrSize; ++j)
//        {
//            int i;
//            for(i = SurrSize; i >= -SurrSize; --i)
//            {
//                In.x = i;
//                In.y = j;
//                SNormalPatternRotate180[pcnt] = In;
//                ++pcnt;
//            }
//        }
//        
//        Patt Snake[7];
//        
//        Snake[0].x = 0;
//        Snake[0].y = 1;
//        Snake[1].x = 0;
//        Snake[1].y = 1;
//        Snake[2].x = 1;
//        Snake[2].y = 0;
//        Snake[3].x = 1;
//        Snake[3].y = 0;
//        Snake[4].x = 0;
//        Snake[4].y = 1;
//        Snake[5].x = 0;
//        Snake[5].y = 1;
//        Snake[6].x = 1;
//        Snake[6].y = 0;
//        
//        int ChooseSnake = 0;
//        int SnakeCount = 0;
//        int SurroundCounter = 0;
//        int Xstart = -SurrSize * 2;
//        int Ystart = 0;
//        
//        Patt StartPattern[2];
//        
//        StartPattern[0].x = 1;
//        StartPattern[0].y = 0;
//        StartPattern[1].x = 2;
//        StartPattern[1].y = -1;
//
//        int StartOrder[16] =
//            { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 };
//        int StartOrderCnt = 0;
//        
//        for(i = 0; i < 2 * SurrSize + 1; ++i)
//        {
//            ++ChooseSnake;
//            if(ChooseSnake > 4)
//                ChooseSnake = 1;
//            // Define the first coordinate pair
//            if(i == 0)
//            {
//                BridgeTypeOne[SurroundCounter].x = Xstart;
//                BridgeTypeOne[SurroundCounter].y = Ystart;
//            }
//            else
//            {
//                BridgeTypeOne[SurroundCounter].x =
//                    Xstart + StartPattern[StartOrder[StartOrderCnt]].x;
//                BridgeTypeOne[SurroundCounter].y =
//                    Ystart + StartPattern[StartOrder[StartOrderCnt]].y;
//                Xstart = BridgeTypeOne[SurroundCounter].x;
//                Ystart = BridgeTypeOne[SurroundCounter].y;
//                ++StartOrderCnt;
//            }
//            ++SurroundCounter;
//
//            for(int j = 0; j < 2 * SurrSize; j += 4)
//            {
//                switch (ChooseSnake)
//                {               // Resetting Snake
//                    case 1:
//                        SnakeCount = 1;
//                        break;
//                    case 2:
//                        SnakeCount = 2;
//                        break;
//                    case 3:
//                        SnakeCount = 3;
//                        break;
//                    case 4:
//                        SnakeCount = 0;
//                        break;
//                }
//                for(int DoSnake = 0; DoSnake < 4; ++DoSnake)
//                {               // Doing Snake
//                    BridgeTypeOne[SurroundCounter].x =
//                        BridgeTypeOne[SurroundCounter - 1].x +
//                        Snake[SnakeCount].x;
//                    BridgeTypeOne[SurroundCounter].y =
//                        BridgeTypeOne[SurroundCounter - 1].y +
//                        Snake[SnakeCount].y;
//                    ++SnakeCount;
//                    ++SurroundCounter;
//                }
//            }
//        }
//        Snake[0].x = 0;
//        Snake[0].y = -1;
//        Snake[1].x = 0;
//        Snake[1].y = -1;
//        Snake[2].x = 1;
//        Snake[2].y = 0;
//        Snake[3].x = 1;
//        Snake[3].y = 0;
//        Snake[4].x = 0;
//        Snake[4].y = -1;
//        Snake[5].x = 0;
//        Snake[5].y = -1;
//        Snake[6].x = 1;
//        Snake[6].y = 0;
//
//        ChooseSnake = 0;
//        SnakeCount = 0;
//        SurroundCounter = 0;
//        Xstart = SurrSize;
//        Ystart = SurrSize;
//
//        StartPattern[0].x = -1;
//        StartPattern[0].y = 0;
//        StartPattern[1].x = -2;
//        StartPattern[1].y = -1;
//        StartOrderCnt = 0;
//
//        for(i = 0; i < 2 * SurrSize + 1; ++i)
//        {
//            ++ChooseSnake;
//            if(ChooseSnake > 4)
//                ChooseSnake = 1;
//            // Define the first coordinate pair
//            if(i == 0)
//            {
//                BridgeTypeTwo[SurroundCounter].x = Xstart;
//                BridgeTypeTwo[SurroundCounter].y = Ystart;
//            }
//            else
//            {
//                BridgeTypeTwo[SurroundCounter].x =
//                    Xstart + StartPattern[StartOrder[StartOrderCnt]].x;
//                BridgeTypeTwo[SurroundCounter].y =
//                    Ystart + StartPattern[StartOrder[StartOrderCnt]].y;
//                Xstart = BridgeTypeTwo[SurroundCounter].x;
//                Ystart = BridgeTypeTwo[SurroundCounter].y;
//                ++StartOrderCnt;
//            }
//            ++SurroundCounter;
//
//            for(int j = 0; j < 2 * SurrSize; j += 4)
//            {
//                switch (ChooseSnake)
//                {               // Resetting Snake
//                    case 1:
//                        SnakeCount = 1;
//                        break;
//                    case 2:
//                        SnakeCount = 0;
//                        break;
//                    case 3:
//                        SnakeCount = 3;
//                        break;
//                    case 4:
//                        SnakeCount = 2;
//                        break;
//                }
//                for(int DoSnake = 0; DoSnake < 4; ++DoSnake)
//                {               // Doing Snake
//                    BridgeTypeTwo[SurroundCounter].x =
//                        BridgeTypeTwo[SurroundCounter - 1].x +
//                        Snake[SnakeCount].x;
//                    BridgeTypeTwo[SurroundCounter].y =
//                        BridgeTypeTwo[SurroundCounter - 1].y +
//                        Snake[SnakeCount].y;
//                    ++SnakeCount;
//                    ++SurroundCounter;
//                }
//            }
//        }
//
//        /*	
//        int ModJ = int (0.57735027 * SurrSize) + 1;
//        int DiffJ = SurrSize - ModJ; assert(DiffJ >= 0);
//        int Trim = (2*SurrSize+1) * DiffJ;
//        pcnt = 0; 
//        for(i = 0; i < (2*ModJ+1)*(2*SurrSize+1); ++i) {
//
//	    BridgeTypeTwo[pcnt].x = BridgeTypeTwo[Trim].x;
//	    BridgeTypeTwo[pcnt].y = BridgeTypeTwo[Trim].y;
//	    BridgeTypeOne[pcnt].x = BridgeTypeOne[Trim].x;
//	    BridgeTypeOne[pcnt].y = BridgeTypeOne[Trim].y;
//	    SNormalPattern[pcnt].x = SNormalPattern[Trim].x;
//	    SNormalPattern[pcnt].y = SNormalPattern[Trim].y;
//	    SNormalPatternRotate180[pcnt].x = SNormalPatternRotate180[Trim].x;
//	    SNormalPatternRotate180[pcnt].y = SNormalPatternRotate180[Trim].y;
//	    ++pcnt; ++Trim;
//        }
//        */
//
//        // Connecting This Surface
//
//
//        Patt IndexPattern[3][300];
//        int Index[3] = { 0, 0, 0 };
//        
//        for(int coord = 1; coord < 4; ++coord)
//        {
//
//            GridSite *MySite = LocateCentralSite(coord);
//            Eshort CurrPos[3];
//
//            MySite->GetGridPos(CurrPos);
//            x = CurrPos[0];
//            y = CurrPos[1];
//            int MyCoord = FindCoordination(x, y);
//            int RealSiteCntr = 0;
//
//            for(local = 0; local < pcnt; ++local)
//            {
//                int dx, dy;
//
//                dx = SNormalPattern[local].x;
//                dy = SNormalPattern[local].y;
//
//                GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
//                GridSite *SurrSite = NULL;
//
//                if(FinalX < 0 || FinalY < 0 || FinalX > theGrid->EDim[0] - 1
//                    || FinalY > theGrid->EDim[1] - 1)
//                {
//                }
//                else
//                    SurrSite = Surface[AdsorbateLayer][FinalX][FinalY];
//
//                if(SurrSite != NULL)
//                {
//                    double Distance = (double)MySite->FindDistance(SurrSite, GridLength);
//                    if(Distance <= Cutoff * MMdistance)
//                    {
//                        if(RealSiteCntr == PatternDim)
//                        {
//                            exit(1);
//                        }
//                        
//                        // Saving Distances
//                        SiteDistances[RealSiteCntr][MyCoord] = (float)Distance;
//                        ++RealSiteCntr;
//                        
//                        if(coord == 1)
//                        {
//                            IndexPattern[0][Index[0]] =
//                                SNormalPattern[local];
//                            SymmetryPattern[Index[0]][0] =
//                                PatternsNormalAtop[local]; // Symmetry pattern might not be used anywhere
//                            ++Index[0]; // assert(Index[0] < PatternDim);   
//                        }
//                        else if(coord == 2)
//                        {
//                            IndexPattern[1][Index[1]] =
//                                SNormalPattern[local];
//                            BridgeTypeOne[Index[1]] = BridgeTypeOne[local];
//                            BridgeTypeTwo[Index[1]] = BridgeTypeTwo[local];
//                            SymmetryPattern[Index[1]][1] =
//                                PatternsNormalBridge[local];
//                            ++Index[1]; // assert(Index[1] < PatternDim);      
//                        }
//                        else if(coord == 3)
//                        {
//                            IndexPattern[2][Index[2]] =
//                                SNormalPattern[local];
//                            SNormalPatternRotate180[Index[2]] =
//                                SNormalPatternRotate180[local];
//                            SymmetryPattern[Index[2]][2] =
//                                PatternsNormalHollow[local];
//                            ++Index[2];
//                            assert(Index[2] < PatternDim);
//                        }
//                    }
//                }
//            }
//            DefinedDist[FindCoordination(x, y)] = true; // DWDWDW
//        }
//
//        for(x = 0; x < EDim[0]; ++x)
//        {
//            for(y = 0; y < EDim[1]; ++y)
//            {
//                GridSite *MySite = Surface[AdsorbateLayer][x][y];
//
//                for(local = 0; local < MetalIndex; ++local)
//                {
//                    dx = MetalPattern[local].x;
//                    dy = MetalPattern[local].y;
//                    GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
//                    if(FindCoordination(x, y) == Atop)
//                    {           // Connecting M-M bonds on surface
//                        Surface[MetalLayer][x][y]->
//                            AddMetalBond(Surface[MetalLayer][FinalX]
//                                [FinalY]);
//                    }
//                }
//                for(local = 0; local < SiteIndex; ++local)
//                {
//                    dx = Pattern[local].x;
//                    dy = Pattern[local].y;
//
//                    GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
//
//                    if(FindCoordination(x, y) == Atop)
//                    {           // Connecting Metal 
//                        // Atop Sites with Overlayer
//                        Surface[MetalLayer][x][y]->
//                            AddBond(Surface[AdsorbateLayer][FinalX]
//                                [FinalY]);
//                        Surface[AdsorbateLayer][FinalX][FinalY]->
//                            AddBond(Surface[MetalLayer][x][y]);
//                        if(FindCoordination(FinalX, FinalY) == Atop)
//                        {
//                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
//                                BondSize > 1)
//                            {
//                                cout << "\nError in Atop Definition\n" <<
//                                    flush;
//                            }
//                        }
//                        else if(FindCoordination(FinalX, FinalY) == Bridge)
//                        {
//                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
//                                BondSize > 2)
//                            {
//                                cout << "\nError in Bridge Definition\n" <<
//                                    flush;
//                            }
//                        }
//                        else if(FindCoordination(FinalX,
//                                    FinalY) == Hollow_3)
//                        {
//                            if(Surface[AdsorbateLayer][FinalX][FinalY]->
//                                BondSize > 4)
//                            {
//                                cout << "\nError in Hollow Definition\n" <<
//                                    flush;
//                            }
//                        }
//                    }
//                }
//                if(MySite != NULL && FindCoordination(x, y) != NONE)
//                {
//                    Patt *CurrentPattern = NULL;
//                    CharString WhichSite = DistinguishSites(x, y);
//                    int indexval = FindCoordination(x, y) - 1;
//
//                    if(WhichSite == "Bridge+30")
//                    {
//                        CurrentPattern = BridgeTypeOne;
//                    }
//                    else if(WhichSite == "Bridge-30")
//                    {
//                        CurrentPattern = BridgeTypeTwo;
//                    }
//                    else if(WhichSite == "Hollow+180")
//                    {
//                        CurrentPattern = SNormalPatternRotate180;
//                    }
//                    else
//                    {
//                        CurrentPattern = IndexPattern[indexval];
//                    }
//                    for(local = 0; local < Index[indexval]; ++local)
//                    {
//                        int dx, dy;
//                        dx = CurrentPattern[local].x;
//                        dy = CurrentPattern[local].y;
//                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
//                        GridSite *SurrSite = Surface[AdsorbateLayer][FinalX][FinalY];
//                        assert(SurrSite != NULL);
//                        MySite->Surround[MySite->SurroundSize] = SurrSite;
//                        ++MySite->SurroundSize;
//                        assert(MySite->SurroundSize < PatternDim);
//                    }
//                }
//            }
//        }
//    }
    else if(SurfaceType.lower() == "graphene")
    {
        Patt SitePattern[6][300]; // Matrix containing the search pattern of various site types
        Patt MetalMetal[2][3]; 
        Patt MetalAdsor[2][7];
        
        MetalMetal[0][0].x = 2;
        MetalMetal[0][0].y = 0;
        MetalMetal[0][1].x = -2;
        MetalMetal[0][1].y = -2;
        MetalMetal[0][2].x = -2;      // Connecting metals and metals for metal type 1
        MetalMetal[0][2].y = 2;
        
        MetalMetal[1][0].x = -2;
        MetalMetal[1][0].y = 0;
        MetalMetal[1][1].x = 2;
        MetalMetal[1][1].y = -2;      // Connecting metals and metals for metal type 2
        MetalMetal[1][2].x = 2;
        MetalMetal[1][2].y = 2;
        MetalIndex = 3;
        
        MetalAdsor[0][0].x = 0;
        MetalAdsor[0][0].y = 0;
        MetalAdsor[0][1].x = 1;
        MetalAdsor[0][1].y = 2;       
        MetalAdsor[0][2].x = 1;
        MetalAdsor[0][2].y = 0;       // Connecting metals and adsorbates for metal type 1
        MetalAdsor[0][3].x = 1;
        MetalAdsor[0][3].y = -2;
        MetalAdsor[0][4].x = -1;
        MetalAdsor[0][4].y = -1;
        MetalAdsor[0][5].x = -3;
        MetalAdsor[0][5].y = 0;
        MetalAdsor[0][6].x = -1;
        MetalAdsor[0][6].y = 1;
        
        MetalAdsor[1][0].x = 0;
        MetalAdsor[1][0].y = 0;
        MetalAdsor[1][1].x = -1;
        MetalAdsor[1][1].y = 2;        
        MetalAdsor[1][2].x = -1;
        MetalAdsor[1][2].y = 0;       // Connecting metals and adsorbates for metal type 2   
        MetalAdsor[1][3].x = -1;
        MetalAdsor[1][3].y = -2;
        MetalAdsor[1][4].x = 1;
        MetalAdsor[1][4].y = -1;
        MetalAdsor[1][5].x = 3;
        MetalAdsor[1][5].y = 0;
        MetalAdsor[1][6].x = 1;
        MetalAdsor[1][6].y = 1;
        Patt In;
        
        int pcnt; // Number of surrounding sites taken into account
        int SurrSize = 8;
        CharString SiteNames[6];
        SiteNames[0] = "Atop_1";
        SiteNames[1] = "Atop_2";
        SiteNames[2] = "Bridge_1";
        SiteNames[3] = "Bridge_2";    // Types of sites
        SiteNames[4] = "Bridge_3";
        SiteNames[5] = "Hollow";
        int Index[6]; // Number of Surrounding sites for each type of site
        for (int i = 0; i < 6; i++)
        {
          Index[i] = 0;
        }
        for (int k = 0; k < 6; k++)
        {
          GridSite *Central = LocateDistinguishSites(SiteNames[k]);
          Eshort c[3];
          Central->GetGridPos(c);
          pcnt = 0;
          
          for (int i = SurrSize; i >= -SurrSize; i--)
          {
              for (int j = SurrSize; j >= -SurrSize; j--)
              {
                  int FinalX, FinalY;
                  GetRealCoordinates(c[0], c[1], i, j, FinalX, FinalY);
                  GridSite *thisSite = Surface[c[2]][FinalX][FinalY];
                  if (thisSite != NULL && Central->FindDistance(thisSite, GridLength) < Cutoff*MMdistance)
                  {
                    SiteDistances[pcnt][k] = Central->FindDistance(thisSite, GridLength);
                    SitePattern[k][pcnt].x = i;
                    SitePattern[k][pcnt].y = j;
                    ++pcnt;
                    ++Index[k];
                  }
              }
          }
        }
        
        for(x = 0; x < EDim[0]; ++x)
        {
            for(y = 0; y < EDim[1]; ++y)
            {
                GridSite *MySite = Surface[AdsorbateLayer][x][y];
                if (MySite != NULL)
                {
                  // Connecting Metal-Metal Bonds
                  if(DistinguishSites(x, y) == "Atop_1")
                  {
                    for(local = 0; local < MetalIndex; ++local)
                    {
                        dx = MetalMetal[0][local].x;
                        dy = MetalMetal[0][local].y;
                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                        Surface[MetalLayer][x][y]->
                          AddMetalBond(Surface[MetalLayer][FinalX][FinalY]);
                    }
                  }
                  else if (DistinguishSites(x, y) == "Atop_2")
                  {
                    for(local = 0; local < MetalIndex; ++local)
                    {
                        dx = MetalMetal[1][local].x;
                        dy = MetalMetal[1][local].y;
                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                        Surface[MetalLayer][x][y]->
                                AddMetalBond(Surface[MetalLayer][FinalX]
                                    [FinalY]);
                    }
                  }
                  // Connecting Metal-Adsorbate Bonds
                  if(DistinguishSites(x, y) == "Atop_1")
                  {
                      for(local = 0; local < SiteIndex; ++local)
                      {   
                          dx = MetalAdsor[0][local].x;
                          dy = MetalAdsor[0][local].y;
                          GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                          Surface[MetalLayer][x][y]->
                              AddBond(Surface[AdsorbateLayer][FinalX][FinalY]);
                          Surface[AdsorbateLayer][FinalX][FinalY]->
                              AddBond(Surface[MetalLayer][x][y]);
                          
                          if(FindCoordination(FinalX, FinalY) == Atop)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                  BondSize > 1)
                              {
                                  cout << "\nError in Atop Definition\n" << flush;
                                  exit(1);
                              }
                          }
                          else if(FindCoordination(FinalX, FinalY) == Bridge)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                  BondSize > 2)
                              {
                                  cout << "\nError in Bridge Definition\n" << flush;
                                  exit(1);
                              }
                          }
                          else if(FindCoordination(FinalX, FinalY) == Hollow_6)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->
                                  BondSize > 6)
                              {
                                  cout << "\nError in Hollow Definition\n" << flush;
                                  exit(1);
                              }
                          }
                      }
                  }
                  
                  if(DistinguishSites(x, y) == "Atop_2")
                  {
                      for(local = 0; local < SiteIndex; ++local)
                      {   
                          dx = MetalAdsor[1][local].x;
                          dy = MetalAdsor[1][local].y;
                          GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                          Surface[MetalLayer][x][y]->
                              AddBond(Surface[AdsorbateLayer][FinalX][FinalY]);
                          Surface[AdsorbateLayer][FinalX][FinalY]->
                              AddBond(Surface[MetalLayer][x][y]);
                          if(FindCoordination(FinalX, FinalY) == Atop)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->BondSize > 1)
                              {
                                  cout << "\nError in Atop Definition\n" << flush;
                                  exit(1);
                              }
                          }
                          else if(FindCoordination(FinalX, FinalY) == Bridge)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->BondSize > 2)
                              {
                                  cout << "\nError in Bridge Definition\n" << flush;
                                  exit(1);
                              }
                          }
                          else if(FindCoordination(FinalX,
                                      FinalY) == Hollow_6)
                          {
                              if(Surface[AdsorbateLayer][FinalX][FinalY]->BondSize > 6)
                              {
                                  cout << "\nError in Hollow Definition\n" << flush;
                                  exit(1);
                              }
                          }
                      }
                  }
                }
                if(MySite != NULL && FindCoordination(x, y) != NONE)
                {
                    Patt *CurrentPattern = NULL;
                    CharString WhichSite = DistinguishSites(x, y);
                    int indexval;
                    if(WhichSite == "Atop_1")
                    {
                        CurrentPattern = SitePattern[0];
                        indexval = 0;
                    }
                    else if(WhichSite == "Atop_2")
                    {
                        CurrentPattern = SitePattern[1];
                        indexval = 1;
                        
                    }
                    else if(WhichSite == "Bridge_1")
                    {
                        CurrentPattern = SitePattern[2];
                        indexval = 2;
                    }
                    if(WhichSite == "Bridge_2")
                    {
                        CurrentPattern = SitePattern[3];
                        indexval = 3;
                    }
                    else if(WhichSite == "Bridge_3")
                    {
                        CurrentPattern = SitePattern[4];
                        indexval = 4;
                    }
                    else if(WhichSite == "Hollow")
                    {
                        CurrentPattern = SitePattern[5];
                        indexval = 5;
                    }
                    for(local = 0; local < Index[indexval]; ++local)
                    {
                        int dx, dy;
                        dx = CurrentPattern[local].x;
                        dy = CurrentPattern[local].y;
                        GetRealCoordinates(x, y, dx, dy, FinalX, FinalY);
                        GridSite *SurrSite = Surface[AdsorbateLayer][FinalX][FinalY];
                        assert(SurrSite != NULL);
                        MySite->Surround[MySite->SurroundSize] = SurrSite;
                        ++MySite->SurroundSize;
                        assert(MySite->SurroundSize < PatternDim);
                    }
                }
            }
        }
    }
    /* Force metal bond search pattern (array order) to match exactly the 
    * the surrounding site search pattern */
    GridSite *MySite;
    for(x = 0; x < EDim[0]; ++x)
    {
        for(y = 0; y < EDim[1]; ++y)
        {
            MySite = Surface[AdsorbateLayer][x][y];
            if(MySite != NULL)
            {
                GridSite *Metal[2];

                Metal[0] = NULL;
                Metal[1] = NULL;
                Metal[0] = MySite->Bonds[0];
                if(MySite->BondSize != 1)
                    Metal[1] = MySite->Bonds[1];
                else
                {
                    Metal[0]->MetalBonds.ResetToFront();
                    Metal[1] = Metal[0]->MetalBonds.Get();
                }
                double MetalUnitVec[2];
                
                // DW: MetalZPosition was unused, so commenting it out (8/24)
                //double MetalZPosition = Metal[0]->Location.GetZDistance();
                
                Metal[0]->FindUnitVector(Metal[1], theGrid->GridLength,
                    MetalUnitVec);
                double arctangent = atan2(MetalUnitVec[0], MetalUnitVec[1]);
                
                // M_PI is math.h pi to high accuracy
                double Angle = rint(arctangent * (180.0 / M_PI));
                
                // rint rounds to nearest integer, represented as
                // double precision 
                MySite->Angle = Angle;
                Eshort MyPos[3];

                MySite->GetGridPos(MyPos);
                CharString SiteType = DistinguishSites((int)MyPos[0], (int)MyPos[1]);
                if(SiteType == "Atop")
                {
                    MySite->OrientationPos = 0;
                }
                else if(SiteType == "Hollow_2")
                {
                    MySite->OrientationPos = 0;
                }
                else if(SiteType == "Hollow_1")
                {
                    MySite->OrientationPos = 1;
                }
                else if(SiteType == "Bridge_1")
                {
                    MySite->OrientationPos = 0;
                }
                else if(SiteType == "Bridge_2")
                {
                    MySite->OrientationPos = 1;
                }
                else if(SiteType == "Bridge_3")
                {
                    MySite->OrientationPos = 2;
                }
                else if (SiteType == "Atop_1")
                  MySite->OrientationPos = 0;
                else if (SiteType == "Atop_2")
                  MySite->OrientationPos = 1;
                else if (SiteType == "Bridge_1")
                  MySite->OrientationPos = 0;
                else if (SiteType == "Bridge_2")
                  MySite->OrientationPos = 1;
                else if (SiteType == "Bridge_3")
                  MySite->OrientationPos = 2;
                else if (SiteType == "Hollow")
                  MySite->OrientationPos = 0;
            }
        }
    }
    return;
}

void GetRealCoordinates(int x, int y, int dx, int dy, int &FinalX, int &FinalY)
{
    int TransformedX, TransformedY;
    TransformedX = dx + x;
    TransformedY = dy + y;
    
    if(TransformedX < 0)
    {
        FinalX = theGrid->EDim[0] + TransformedX;
    }
    else if(TransformedX > theGrid->EDim[0] - 1)
    {
        FinalX = TransformedX - (theGrid->EDim[0]);
    }
    else
    {
        FinalX = TransformedX;
    }
    if(TransformedY < 0)
    {
        FinalY = theGrid->EDim[1] + TransformedY;
    }
    else if(TransformedY > theGrid->EDim[1] - 1)
    {
        FinalY = TransformedY - (theGrid->EDim[1]);
    }
    else
    {
        FinalY = TransformedY;
    }
    return;
}

CharString DistinguishSites(int i, int j)
{
    if(theGrid->SurfaceType == "100")
    {
        if(i % 2 != 0 && j % 2 != 0)
        {
            return (CharString) "Atop";
        }
        else if(i % 2 == 0 && j % 2 == 0)
        {
            return (CharString) "Hollow";
        }
        else if(i % 2 == 0 && j % 2 != 0)
        {
            return (CharString) "Bridge";
        }
        else if(i % 2 != 0 && j % 2 == 0)
        {
            return (CharString) "Bridge-30";
        }
        else
        {
            return (CharString) "NULL";
        }

    }
    else if(theGrid->SurfaceType == "111")
    {
        if((i + j - 1) % 2 == 0 && j % 2 != 0)
        {
            if((i + j - 1) % 4 == 0)
            {
                return (CharString) "Bridge_1";
            }
            else
            {
                return (CharString) "Atop_1";
            }
        }
        else if((i - 3) % 4 == 0 && j % 4 == 0)
        {
            return (CharString) "Bridge_2";  // Bridge-30 is Bridge_2
        }
        else if((i - 3) % 4 == 0 && (j % 4 != 0 && j % 2 == 0))
        {
            return (CharString) "Bridge_3";  // Bridge+30 is Bridge_3
        }
        else if((i - 1) % 4 == 0 && j % 4 == 0)
        {
            return (CharString) "Bridge_3";
        }
        else if((i - 1) % 4 == 0 && (j % 4 != 0 && j % 2 == 0))
        {
            return (CharString) "Bridge_2";
        }
        else if(i % 2 == 0 && j % 2 == 0)
        {
            if((i + j) % 4 == 0)
                return (CharString) "Hollow_1"; // Hollow+180 is Hollow_1
            else
                return (CharString) "Hollow_2"; // Hollow is Hollow_2
        }
        else
        {
            return (CharString) "NULL";
        }
    }
    else if (theGrid->SurfaceType.lower() == "graphene")
    {
        if (theGrid->FindCoordination(i,j) == Atop && ((i+2)%8 == 0 || (i-2)%8 == 0))
          return (CharString) "Atop_1";
        else if (theGrid->FindCoordination(i,j) == Atop && ((i+4)%8 == 0 || i%8 == 0))
          return (CharString) "Atop_2";
        else if (theGrid->FindCoordination(i,j) == Hollow_6)
          return (CharString) "Hollow";
        else if (theGrid->FindCoordination(i,j) == Bridge && j%2 == 0)
          return (CharString) "Bridge_1";
        else if (theGrid->FindCoordination(i,j) == Bridge && (i+2*j+1)%8 == 0)
          return (CharString) "Bridge_2";
        else if (theGrid->FindCoordination(i,j) == Bridge && (i+2*j+5)%8 == 0)
          return (CharString) "Bridge_3";
        else
        {
          return (CharString) "NULL";   
          cout << "The grid position does not correspond to any site" << endl;
          exit(1);
        }
    }
    
    // DW: Here we would add other cases for other patterns.  As long
    // as the SurfaceType is "100" or "111" we won't make it here.
    
    exit(-2);
    return (CharString) "NULL";
}

// GET SITE COORDINATION 
SiteCoordination Grid::FindCoordination(int i, int j)
{

    if(SurfaceType == "100")
    {
        if(i % 2 != 0 && j % 2 != 0)
        {
            return Atop;
        }
        else if(i % 2 == 0 && j % 2 == 0)
        {
            return Hollow;
        }
        else
        {
            return Bridge;
        }
    }
    else if(SurfaceType == "111")
    {

        if((i + j - 1) % 2 == 0 && j % 2 != 0)
        {
            if((i + j - 1) % 4 == 0)
            {
                return Bridge;
            }
            else
            {
                return Atop;
            }
        }
        else if(i % 2 != 0 && j % 2 == 0)
        {
            return Bridge;
        }
        else if(i % 2 == 0 && j % 2 == 0)
        {
              return Hollow_3;
        }
        else
        {
            return NONE;
        }
    }
    else if(SurfaceType.lower() == "graphene")
    {

        if (j%4 == 2 && ((i+2)%8 == 0 || i%8 == 0))
          return Atop;
        else if (j%4 == 0 && ((i-2)%8 == 0 || (i-4)%8 == 0))
          return Atop;
        else if (j%4 == 0 && (i-3)%8 == 0)
          return Bridge;
        else if (j%4 == 1 && ((i-1)%8 == 0 || (i-5)%8 == 0))
          return Bridge;
        else if (j%4 == 2 && (i+1)%8 == 0)
          return Bridge;
        else if (j%4 == 3 && ((i-1)%8 == 0 || (i-5)%8 == 0))
          return Bridge;
        else if (j%4 == 0 && (i+1)%8 == 0)
          return Hollow_6;
        else if (j%4 == 2 && (i-3)%8 == 0)
          return Hollow_6;
        else
          return NONE;
          
    }
    return NONE;
}

GridSite *Grid::LocateCentralSite(int SiteType)
{                               // Get Coord central site

    assert(EDim[1] > 0);
    int mid[2] = { EDim[0] / 2 - 1, EDim[1] / 2 - 1 };
    int i, j;
    GridSite *FoundSite = NULL;
    GridSite *Check;
    Patt CheckMethod[32];
    int k = 0;
    for (int i = 0; i < 8; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        CheckMethod[k].x = i;
        CheckMethod[k].y = j;
        ++k;
      }
    }
    if(SurfaceType == "100" && SiteType == 3)
        SiteType = 4;
    // Locate a site of the specified type
    for(int k = 0; k < 32 && FoundSite == NULL; ++k)
    {
        i = mid[0] + CheckMethod[k].x;
        j = mid[1] + CheckMethod[k].y;
        Check = Surface[AdsorbateLayer][i][j];
            if (FindCoordination(i,j) == SiteType)
            {
              FoundSite = Check;
              return FoundSite;
            }
            //cout << FindCoordination(i,j) << endl;
//        if(SurfaceType == "100")
//        {
//            if((int)Check->BondSize == SiteType)
//            {
//                FoundSite = Check;
//                return FoundSite;
//            }
//        }
//        else if(SurfaceType == "111")
//        {
//            if (FindCoordination(i,j) == SiteType)
//            {
//              FoundSite = Check;
//              return FoundSite;
//            }
//        }
//        else if(SurfaceType.lower() == "graphene")
//        {
//            if(FindCoordination(i,j) == SiteType)
//            {
//                FoundSite = Check;
//                return FoundSite;
//            }
//        }
    }
    cout << "Check the code to find site, something might be missing" << endl;
    exit(1);
}

GridSite *Grid::LocateDistinguishSites(CharString SiteType)
{   
    assert(EDim[1] > 0);
    int mid[2] = { EDim[0] / 2 - 1, EDim[1] / 2 - 1 };
    int i, j;
    GridSite *FoundSite = NULL;
    GridSite *Check;
    Patt CheckMethod[32];
    int k = 0;

    for (int i = 0; i < 8; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        CheckMethod[k].x = i;
        CheckMethod[k].y = j;
        ++k;
      }
    }
    // Locate a site of the specified type
    for(int k = 0; k < 32 && FoundSite == NULL; ++k)
    {
        i = mid[0] + CheckMethod[k].x;
        j = mid[1] + CheckMethod[k].y;
        Check = Surface[AdsorbateLayer][i][j];
        if (Check!=NULL)
        {
              if(DistinguishSites(i, j) == SiteType)
              {
                  FoundSite = Check;
                  return FoundSite;
              }
              if(DistinguishSites(i, j) == SiteType)
              {
                  FoundSite = Check;
                  return FoundSite;
              }
        }
    }
}

float FindDistance(GridSite * MySite, GridSite * NeighSite)
{

    int Coord = MySite->BondSize;

    // search surrounding for position in array
    for(int i = 0; i < MySite->SurroundSize; ++i)
        if(MySite->Surround[i] == NeighSite)
            return theGrid->SiteDistances[i][Coord];
    return 9999999.;
}

float FindDistance(GridSite * MySite, int index)
{

    int Coord = MySite->BondSize;

    assert(index < MySite->SurroundSize && index >= 0);
    return theGrid->SiteDistances[index][Coord];
}

void Grid::DefineQuadrants()
{

    double eps = 1.e-9;
    GridSite *SurrSite = NULL;
    CharString Surface = GetSurface();
    if(GetSurface() == "111")
    {
        for(int i = 1; i < 4; ++i)
        {
            GridSite *Center = theGrid->LocateCentralSite(i);
            
            EPosition <double> c = Center->GetPosition();
            //cout << c.GetXDistance() << " " << c.GetYDistance() << " " << c.GetZDistance() << endl;
            for(int k = 0; k < Center->SurroundSize; ++k)
            {
                SurrSite = Center->Surround[k];
                if(SurrSite != Center)
                {
                    Eshort Pos[3];

                    Center->GetGridPos(Pos);
                    Eshort NeighPos[3];

                    SurrSite->GetGridPos(NeighPos);
                    double y =
                        Center->Location.GetYDistance() -
                        SurrSite->Location.GetYDistance();
                    double x =
                        Center->Location.GetXDistance() -
                        SurrSite->Location.GetXDistance();
                    double rads = atan2(y, x);
                    double degrees = rads * (180 / (3.1415));

                    if(degrees >= -(180 + eps) && degrees < -(120 + eps))
                    {
                        Quadrant[k][i - 1] = 0;
                    }
                    else if(degrees >= -(120 + eps) &&
                        degrees < -(60 + eps))
                    {
                        Quadrant[k][i - 1] = 1;
                    }
                    else if(degrees >= -(60 + eps) && degrees < -eps)
                    {
                        Quadrant[k][i - 1] = 2;
                    }
                    else if(degrees >= -eps && degrees < (60 - eps))
                    {
                        Quadrant[k][i - 1] = 3;
                    }
                    else if(degrees >= (60 - eps) && degrees < (120 - eps))
                    {
                        Quadrant[k][i - 1] = 4;
                    }
                    else
                        Quadrant[k][i - 1] = 5;
                }
                else
                    Quadrant[k][i - 1] = 0;
            }
        }
    }
    else if(Surface.lower() == "graphene")
    {
        int j;
        for(int i = 1; i < 4; ++i)
        {
            j = i;
            if (j == 3)
            {
              j = 6;
            }
            GridSite *Center = theGrid->LocateCentralSite(j);
            // EPosition <double> c = Center->GetPosition();
            // cout << c.GetXDistance() << " " << c.GetYDistance() << " " << c.GetZDistance() << endl;
            for(int k = 0; k < Center->SurroundSize; ++k)
            {
                SurrSite = Center->Surround[k];
                if(SurrSite != Center)
                {
                    Eshort Pos[3];

                    Center->GetGridPos(Pos);
                    Eshort NeighPos[3];

                    SurrSite->GetGridPos(NeighPos);
                    double y =
                        Center->Location.GetYDistance() -
                        SurrSite->Location.GetYDistance();
                    double x =
                        Center->Location.GetXDistance() -
                        SurrSite->Location.GetXDistance();
                    double rads = atan2(y, x);
                    double degrees = rads * (180 / (3.1415));

                    if(degrees >= -(180 + eps) && degrees < -(120 + eps))
                    {
                        Quadrant[k][i - 1] = 0;
                    }
                    else if(degrees >= -(120 + eps) &&
                        degrees < -(60 + eps))
                    {
                        Quadrant[k][i - 1] = 1;
                    }
                    else if(degrees >= -(60 + eps) && degrees < -eps)
                    {
                        Quadrant[k][i - 1] = 2;
                    }
                    else if(degrees >= -eps && degrees < (60 - eps))
                    {
                        Quadrant[k][i - 1] = 3;
                    }
                    else if(degrees >= (60 - eps) && degrees < (120 - eps))
                    {
                        Quadrant[k][i - 1] = 4;
                    }
                    else
                        Quadrant[k][i - 1] = 5;
                }
                else
                    Quadrant[k][i - 1] = 0;
            }
        }
    }
    else if(GetSurface() == "100")
    {
        for(int i = 1; i < 4; ++i)
        {
            if(i == 3)
                i = 4;
            GridSite *Center = theGrid->LocateCentralSite(i);

            for(int k = 0; k < Center->SurroundSize; ++k)
            {
                SurrSite = Center->Surround[k];
                if(SurrSite != Center)
                {
                    Eshort Pos[3];

                    Center->GetGridPos(Pos);
                    Eshort NeighPos[3];

                    SurrSite->GetGridPos(NeighPos);
                    double y =
                        Center->Location.GetYDistance() -
                        SurrSite->Location.GetYDistance();
                    double x =
                        Center->Location.GetXDistance() -
                        SurrSite->Location.GetXDistance();
                    double rads = atan2(y, x);
                    double degrees = rads * (180 / (3.1415));

                    if(degrees >= -(180 + eps) && degrees < -(90 + eps))
                    {
                        Quadrant[k][i - 1] = 0;
                    }
                    else if(degrees >= -(90 + eps) && degrees < -(eps))
                    {
                        Quadrant[k][i - 1] = 1;
                    }
                    else if(degrees >= -(eps) && degrees < 90 - eps)
                    {
                        Quadrant[k][i - 1] = 2;
                    }
                    else
                        Quadrant[k][i - 1] = 3;
                }
                else
                {
                    Quadrant[k][i - 1] = 0;
                }
            }
        }
    }

}

int SymmetryPattern[PatternDim][3];

int PatternsNormalBridge[289] = {
    56, 0, 57, 0, 58, 0, 59, 0, 60, 0, 59, 0, 58, 0, 57, 0, 56,
    47, 48, 49, 50, 51, 52, 53, 54, 55, 54, 53, 52, 51, 50, 49, 48, 47,
    42, 0, 43, 0, 44, 0, 45, 0, 46, 0, 45, 0, 44, 0, 43, 0, 42,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 40, 39, 38, 37, 36, 35, 34, 33,
    28, 0, 29, 0, 30, 0, 31, 0, 32, 0, 31, 0, 30, 0, 29, 0, 28,
    21, 20, 19, 18, 17, 16, 13, 14, 15, 14, 13, 16, 17, 18, 19, 20, 21,
    27, 0, 26, 0, 8, 0, 10, 0, 4, 0, 10, 0, 8, 0, 26, 0, 27,
    24, 23, 22, 11, 12, 6, 5, 3, 1, 3, 5, 6, 12, 11, 22, 23, 24,
    25, 0, 9, 0, 7, 0, 2, 0, 0, 0, 2, 0, 7, 0, 9, 0, 25,
    24, 23, 22, 11, 12, 6, 5, 3, 1, 3, 5, 6, 12, 11, 22, 23, 24,
    27, 0, 26, 0, 8, 0, 10, 0, 4, 0, 10, 0, 8, 0, 26, 0, 27,
    21, 20, 19, 18, 17, 16, 13, 14, 15, 14, 13, 16, 17, 18, 19, 20, 21,
    28, 0, 29, 0, 30, 0, 31, 0, 32, 0, 31, 0, 30, 0, 29, 0, 28,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 40, 39, 38, 37, 36, 35, 34, 33,
    42, 0, 43, 0, 44, 0, 45, 0, 46, 0, 45, 0, 44, 0, 43, 0, 42,
    47, 48, 49, 50, 51, 52, 53, 54, 55, 54, 53, 52, 51, 50, 49, 48, 47,
    56, 0, 57, 0, 58, 0, 59, 0, 60, 0, 59, 0, 58, 0, 57, 0, 56
};

int PatternsNormalHollow[289] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, 35, 31, 36, 37, 36, 31, 35, -1, -1, -1, -1, -1,
    34, 0, 33, 0, 32, 0, 30, 0, 24, 0, 30, 0, 32, 0, 33, 0, 34,
    29, 28, 27, 26, 19, 15, 23, 18, 17, 18, 23, 15, 19, 26, 27, 28, 29,
    25, 0, 20, 0, 16, 0, 10, 0, 11, 0, 10, 0, 16, 0, 20, 0, 25,
    22, 21, 14, 13, 12, 9, 6, 8, 7, 8, 6, 9, 12, 13, 14, 21, 22,
    20, 0, 13, 0, 5, 0, 4, 0, 2, 0, 4, 0, 5, 0, 13, 0, 20,
    19, 16, 12, 9, 6, 4, 3, 1, 0, 1, 3, 4, 6, 9, 12, 16, 19,
    15, 0, 10, 0, 8, 0, 2, 0, 1, 0, 2, 0, 8, 0, 10, 0, 15,
    23, 18, 17, 11, 7, 8, 6, 4, 3, 4, 6, 8, 7, 11, 17, 18, 23,
    24, 0, 18, 0, 10, 0, 9, 0, 5, 0, 9, 0, 10, 0, 18, 0, 24,
    31, 30, 23, 15, 19, 16, 12, 13, 14, 13, 12, 16, 19, 15, 23, 30, 31,
    35, 0, 32, 0, 26, 0, 20, 0, 21, 0, 20, 0, 26, 0, 32, 0, 35,
    -1, -1, -1, 33, 27, 28, 29, 25, 22, 25, 29, 28, 27, 33, -1, -1, -1,
    -1, 0, -1, 0, 34, 0, 38, 0, 39, 0, 38, 0, 34, 0, -1, 0, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

int PatternsNormalAtop[289] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, 23, 20, -1, -1, -1, -1, -1, -1, -1, 20, 23, -1, -1, -1,
    23, 0, 22, 0, 18, 0, 19, 0, 21, 0, 19, 0, 18, 0, 22, 0, 23,
    20, 18, 16, 17, 16, 15, 12, 14, 13, 14, 12, 15, 16, 17, 16, 18, 20,
    19, 0, 15, 0, 9, 0, 11, 0, 10, 0, 11, 0, 9, 0, 15, 0, 19,
    13, 14, 12, 11, 8, 7, 8, 6, 5, 6, 8, 7, 8, 11, 12, 14, 13,
    14, 0, 10, 0, 6, 0, 3, 0, 4, 0, 3, 0, 6, 0, 10, 0, 14,
    12, 11, 8, 6, 5, 4, 2, 1, 2, 1, 2, 4, 5, 6, 8, 11, 12,
    9, 0, 7, 0, 3, 0, 1, 0, 0, 0, 1, 0, 3, 0, 7, 0, 9,
    12, 11, 8, 6, 5, 4, 2, 1, 2, 1, 2, 4, 5, 6, 8, 11, 12,
    14, 0, 10, 0, 6, 0, 3, 0, 4, 0, 3, 0, 6, 0, 10, 0, 14,
    13, 14, 12, 11, 8, 7, 8, 6, 5, 6, 8, 7, 8, 11, 12, 14, 13,
    19, 0, 15, 0, 9, 0, 11, 0, 10, 0, 11, 0, 9, 0, 15, 0, 19,
    20, 18, 16, 17, 16, 15, 12, 14, 13, 14, 12, 15, 16, 17, 16, 18, 20,
    23, 0, 22, 0, 18, 0, 19, 0, 21, 0, 19, 0, 18, 0, 22, 0, 23,
    -1, -1, -1, 23, 20, -1, -1, -1, -1, -1, -1, -1, 20, 23, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};
