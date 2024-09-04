#include <cmath>
#include <cassert>
#include "Model.h"
#include "Grid.h"
#include "Species.h"
#include "Component.h"
#include "Interactions.h"
#include "GridSite.h"
#include "Geometry.h"
#include "Base/Estring.h"
#include "ShortRoutines.h"
#include "Simulation.h"
#include "Reactions.h"
#include "SystemHeaders.h"


Matrix<double> ThroughSpaceInteraction[89][12][12][4][10];

/* Old Simulation Routines, Not Supported although they may
	be available in the code.
	Cleaning House */

// Function Prototypes for Model

/* return the BOC function for a molecule */
static double BOC_Molecule_Function(const Component * ThisType, const Component * MetalType, const Eshort Coord);




/* return the BOC function for an atom */
static double BOC_Function(const double &E, const Eshort & Coordination);




/* BOC Conservation calculation of thisatom */
static void BondOrder(GridSite * ThisAtom, double[]);




/* BOC */
/* DW 8/24/04: Note that I removed the const part of arg 6 */
static double BondOrderConservation(const double PairEnergy[],
    Eshort & NumberOfPairs, const Eshort & AdsorbateCoordination,
    const Eshort NeighborCoordination[], const Component * ThisType,
    /*const*/ Component * NeighType[], const Component * MetalType);



/* static double MetalAtomCenterEnergy(GridSite * MetalAtom); */

/* DW 8/24/04: Note that I removed the const part of arg 4 */
/* static double MetalAtomCenterConservation(const double PairEnergy[],
    Eshort & NumberOfPairs, const Eshort NeighborCoordination[],
    Component * NeighType[], const Component * MetalType); */


#include "backprop.cpp"

Matrix<double> Molecule_Metal;
Model *theModel;

void SetComponentEnergiesForAlloy();

Model::Model(fstream &fin, fstream &fout) : Facilitator(fin, fout)
{
    ::theModel = this;
    // DW 9/1/04: I think this 'theInteractions' is the global one!
    theInteractions = theSpecies->theInteractions;
    
    theModel = BOC;
    RadialMethod = NoCoreRepulsion;
    Initialize((CharString) "Model", fin);
    
    Dielectric = new double;
    DistDepend = new double;

    *DistDepend = 1.0;
    *Dielectric = 10.0;
    
    int LayerSize[3] = { 9, 4, 1 };
    int NumberOfLayers = 3;
    Read(fin);
    InitializeThroughSpaceMatrix();
    NumberOfOrientations[0] = 12;
    NumberOfOrientations[1] = 2;
    NumberOfOrientations[2] = 3;
    NumberOfOrientations[3] = 4;
    NumberOfOrientations[5] = 6;
    
    for (int i = 0; i < 25; ++i)
    {
      for (int j = 0; j < 10; ++j)
      {
        for (int k = 0; k < 20; ++k)
        {
          for (int l = 0; l < 100; ++l)
          {
            for (int m = 0; m < 3; ++m)
            {
              for (int n = 0; n < 6; ++n)  
              {
                TempCoord[n][i][j][k][l][m] = 0.0; // n: BondSize - 1; i: Component->_idhook; j: OrientationPos; k: Orientation; l:Atom Number; m: X (0) or Y (1) or Z (2)
              }
            }
          }
        }
      }
    }
}

void SetComponentEnergiesForAlloy()
{
    double SC = theSimulation->SurfaceCharge;
    // interaction routine uses sim temperature which is initially
    // not defined
    DistanceType InitialDistanceType = theModel->RadialMethod;
    theModel->RadialMethod = NoCoreRepulsion;
    int ind = theSimulation->FindChargeIndex(SC);
    
    for(int i = 0; i < NumberOfSpecies; ++i)
    {
//        cout << "\n\n\t" << Listing[i]->GetName().Ec_str() << endl;
//        cout << setw(8) << " Site " << setw(8) << " E " << setw(8) <<
//            " BOC_E " << setw(8) << " BOC_E_Alloy " << endl;
//        
        for(int j = 0; j < 5; ++j)
        {
            if (theGrid->GetSurface().lower() == "graphene" && j == 3)
                j = 6;
            if(Listing[i]->BindingEnergy[j][0][ind] > 0.001)
            {                   // if specified binding energy, not shut down
                GridSite *Site = theGrid->LocateCentralSite(j);
                
                if(Site != NULL && theGrid->Alloy_Type != theNULLSpecies)
                {
                    
                    double SavePrimaryMetalBindingEnergy =
                        Listing[i]->BindingEnergy[j][0][ind];
                    
                    bool SaveOnSurface = Listing[i]->OnSurface;
                    Component *OriginalMetal[10];
                    
                    Listing[i]->OnSurface = true; // DWDWDW
                    Listing[i]->BindingEnergy[j][0][ind] = -1.0;
                    int ii = 0;
                    
                    for(ii = 0; ii < Site->BondSize; ++ii)
                    {
                        OriginalMetal[ii] = Site->Bonds[ii]->GetType();
                        Site->Bonds[ii]->SetType(theGrid->Alloy_Type);
                    }
                    Site->SetType(Listing[i]);
                    double BOC_AlloyEnergy =
                        theModel->CalculateBindingEnergy(Site);
                    for(ii = 0; ii < Site->BondSize; ++ii)
                    {
                        Site->Bonds[ii]->SetType(theGrid->Metal);
                    }
                    double BOC_MetalEnergy =
                        theModel->CalculateBindingEnergy(Site);
                    Site->SetType(theNULLSpecies);
                    for(ii = 0; ii < Site->BondSize; ++ii)
                    {
                        Site->Bonds[ii]->SetType(OriginalMetal[ii]);
                    }
                    Listing[i]->BindingEnergy[j][0][ind] =
                        SavePrimaryMetalBindingEnergy;
                    Listing[i]->BindingEnergy[j][1][ind] =
                        SavePrimaryMetalBindingEnergy * BOC_AlloyEnergy / BOC_MetalEnergy;
                    Listing[i]->OnSurface = SaveOnSurface;
//                    cout << setw(8) << j << " " << setw(8) <<
//                        SavePrimaryMetalBindingEnergy << " " << setw(8) <<
//                        BOC_MetalEnergy << " " << setw(8) << BOC_AlloyEnergy
//                        << endl;
                }
            }
        }
    }
    theModel->RadialMethod = InitialDistanceType;
    return;
}

double *InteractionAddress;
void Model::CheckInteractions()
{
    if (theGrid->GetSurface().lower() == "graphene")
    {
      CharString SiteNames[6]; // stored in stack
      //CharString *SiteNames = new CharString[6];
      SiteNames[0] = "Atop_1";
      SiteNames[1] = "Atop_2";
      SiteNames[2] = "Bridge_1";
      SiteNames[3] = "Bridge_2";    // Types of sites
      SiteNames[4] = "Bridge_3";
      SiteNames[5] = "Hollow";
      int Index[6]; // stored in stack // Number of Surrounding sites for each type of site
      for (int i = 0; i < 6; i++)
      {
        Index[i] = 0;
      }
      for (int k = 0; k < 6; k++)
      {
        if (theGrid->Atop_a == 0 && (k == 0 || k == 1))
          continue;
        if (theGrid->Bridge_a == 0 && (k == 2 || k == 3 || k == 4))
          continue;
        if (theGrid->Hollow_a == 0 && k == 5)
          continue;
        GridSite *Central = theGrid->LocateDistinguishSites(SiteNames[k]);
        Component * CentralComp = Central->GetType();
        int tor = Central->Orientation;
        for (int i = 0; i < NumberOfSpecies; ++i)
        {
          if (Listing[i]->OnSurface == true)
          {
            Central->SetType(Listing[i]);
            for (int j = 0; j < NumberOfOrientations[Central->BondSize-1]; ++j)
            {
              Central->Orientation = j;
              GetOrientationPosition(Central);
            }
          }
        }
        Central->Orientation = tor;
        Central->SetType(CentralComp);
      }
    }
    else if (theGrid->GetSurface().lower() == "111")
    {
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
        //cout << theGrid->Atop_a << " " << theGrid->Bridge_a << " " << theGrid->Hollow_a << endl;
        if (theGrid->Atop_a == 0 && (k == 0))
          continue;
        if (theGrid->Bridge_a == 0 && (k == 1 || k == 2 || k == 3))
          continue;
        if (theGrid->Hollow_a == 0 && (k == 4 || k == 5))
          continue;
        //cout << "Hi" << endl;
        GridSite *Central = theGrid->LocateDistinguishSites(SiteNames[k]);
        Component * CentralComp = Central->GetType();
        int tor = Central->Orientation;
        for (int i = 0; i < NumberOfSpecies; ++i)
        {
          if (Listing[i]->OnSurface == true)
          {
            Central->SetType(Listing[i]);
            for (int j = 0; j < NumberOfOrientations[Central->BondSize-1]; ++j)
            {
               Central->Orientation = j;
               GetOrientationPosition(Central);
            }
          }
        }
        Central->Orientation = tor;
        Central->SetType(CentralComp);
      }
    }
    return;
}

void Model::GetOrientationPosition(GridSite *Site)
{
    double SC = theSimulation->SurfaceCharge;
    assert(Site->Type != theNULLSpecies);
    double Origin[3] = { 0. }; // stored in stack
    int Coordination = Site->BondSize; // stored in stack
    int MyType = Site->Type->_idhook; // stored in stack

    EPosition < double > Location = Site->GetPosition(); // stored in stack
    
    Origin[0] = Location.GetXDistance();
    Origin[1] = Location.GetYDistance();
    Origin[2] = Location.GetZDistance();
    
    Eshort SiteID = Listing[MyType]->_id; // stored in stack
    
    theSpecies->theGeometries->theCoordinates.ResetToFront();
    Coordinates SiteCoords; // stored in stack

    bool Finish = false; // DWDWDW  // DWDWDW // stored in stack
    
    while(theSpecies->theGeometries->theCoordinates && !Finish)
    {
          if(SiteID == theSpecies->theGeometries->theCoordinates().SpeciesName
              && Site->BondSize == theSpecies->theGeometries->theCoordinates().MetalCoord 
              && theSpecies->theGeometries->theCoordinates().SurfaceCharge == SC)
          {
              SiteCoords = theSpecies->theGeometries->theCoordinates();
              Finish = true; // DWDWDW
          }
        
        ++theSpecies->theGeometries->theCoordinates;
    }
    
    if(!(Finish))
    {
       //cout << "Site Not found in Geometry. Check again. " << Site->Type->GetName() << " " << Site->BondSize << endl;
       //cout << Site->Type->GetName() << " " << Site->BondSize << endl;
       return;
    }
    else
    {
      //cout << Listing[MyType]->GetName() << " " << Site->BondSize << endl;
    }
    
    double LocPointer[3] = { 0. };   // stored in stack   // Exact location of atom 1 
    double x, y, z; // stored in stack
    
    double SiteA_Angle = Site->Angle; // stored in stack
    double AdsA_Angle = 0.0; // stored in stack
    double A_Net_Angle = 0.0; // stored in stack
    double MetalZPosition = Site->Bonds[0]->Location.GetZDistance(); // stored in stack
    double GridZPosition = Site->Location.GetZDistance(); // stored in stack
    int MyOrientation = Site->Orientation; // stored in stack
    int MyOrientationPos = Site->OrientationPos; // stored in stack
    
    if(Coordination == 2)
    {
        // rotation by 180 degrees for bridge site 
        AdsA_Angle = double (MyOrientation) * 180.0;     
    }
    else if(Coordination == 3)
    {
        AdsA_Angle = double (MyOrientation) * 120.0;
    }
    else if(Coordination == 4)
    {
        AdsA_Angle = double (MyOrientation) * 90.0;
    }
    else if(Coordination == 6)
    {
        AdsA_Angle = double (MyOrientation) * 60.0;
    }
    else
    {
        AdsA_Angle = double (MyOrientation) * 30.0;
    }
    
    A_Net_Angle = (SiteA_Angle + AdsA_Angle);
    while(A_Net_Angle >= 360.0)
    {
        A_Net_Angle = A_Net_Angle - 360.0;
    }
    while(A_Net_Angle < 0.0)
    {
        // place angle between 0 and 360
        A_Net_Angle = 360. + A_Net_Angle;
    }
    
    int no = 0; // stored in stack
    Site->NumAtoms = SiteCoords.NumAtoms;
    /* Looping through atoms of the adsorbate in the center site */
    for(int iNum = 0; iNum < SiteCoords.NumAtoms; ++iNum)
    {
        x = SiteCoords.Position[iNum].GetXDistance();
        y = SiteCoords.Position[iNum].GetYDistance();
        z = SiteCoords.Position[iNum].GetZDistance();
        double newarctangent = atan2(x, y) + (A_Net_Angle) * (M_PI / 180.); // stored in stack
        double measure = pow((x * x + y * y), 0.5); // stored in stack
        LocPointer[0] = Origin[0] + sin(newarctangent) * measure;
        LocPointer[1] = Origin[1] + cos(newarctangent) * measure;
        LocPointer[2] = z + GridZPosition;
        TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][0] = LocPointer[0] - Origin[0];
        TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][1] = LocPointer[1] - Origin[1];
        TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][2] = LocPointer[2] - Origin[2];
        //cout << Site->Type->GetName() << " " << MyOrientation << " " << LocPointer[0] - Origin[0] << " " << LocPointer[1] - Origin[1] << " " << LocPointer[2] - Origin[2] << endl;
    }
    return;
}

void Model::Read(fstream &fin)
{
    fin.seekg(inp);
    CharString Next;
    CharString End = "END";
    
    while(fin >> Next && Next.lower() != End.lower())
    {
        if(Next.lower() == "radial")
        {
            fin >> Next;
            if(Next.lower() == "none")
            {
                RadialMethod = NoCoreRepulsion;
            }
            else if(Next.lower() == "hard")
            {
                RadialMethod = HardSphere;
            }
            else if(Next.lower() == "mmff")
            {
                RadialMethod = AtomAtomMMFF;
            }
            else
              cout << "Specified method does not exist" << endl;

        }
        if(Next.lower() == "boc")
        {
            ChangeModel(BOC);
        }
        if(Next.lower() == "none")
        {
            ChangeModel(NoModel);
        }
        else if(Next.lower() == "dielectric")
        {
            double val;
            
            fin >> val;
            *Dielectric = val;
        }
        else if(Next.lower() == "distance_dependence")
        {
            double val;

            fin >> val;
            *DistDepend = val;
        }
    }
    
}

double TSE;
double Model::CalculateBindingEnergy(GridSite * PassedSite, bool FixedOrientation) // DWDWDW
{
    // This function must calculate the binding energy at PassedSite 
    // This function must set PassedSite->BindingEnergy if site is occupied
    double Val = 0.0;

    if(PassedSite == NULL)
        return 0.0;
    if(PassedSite->GetType() == theNULLSpecies)
        return 0.0;
    TSE = 0.0;

    if(RadialMethod != NoCoreRepulsion)
    {
        // double Test =
        //     PassedSite->Type->BindingEnergy[PassedSite->BondSize][0];
        // if(Test < 0.01) {
        // PassedSite->BindingEnergy = -1.e20;
        // return -1.e20;
        // }
        if(FixedOrientation)
        {
            TSE = CalculateThroughSpaceFixedOrientation(PassedSite);
        }
        else
        {
            TSE = CalculateThroughSpace(PassedSite);
        }
    }
    double theEnergies[2];
    
    switch (theModel)
    {
        case BOC:
        {
            BondOrder(PassedSite, theEnergies);
            Val = theEnergies[1] - TSE;
//            if (theSimulation->SurfaceCharge == 0 && PassedSite->GetType()->GetName() == "NH2")
//            {
//              cout << PassedSite->GetType()->GetName() << " BOC: " << theEnergies[1] << " Interaction: " << TSE << endl;
//            }
            // if(Val <= 0.) Val = -0.01;
            //PassedSite->BindingEnergy = Val;
            return Val;
        }
        default:
        {
            //PassedSite->BindingEnergy = 0.0;
            return -0.01;
        }
    }
}

double Model::CalculateInteractionEnergy(GridSite * PassedSite)
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    if(PassedSite != NULL)
    {
        thisComponent = PassedSite->GetType();
        double TSE = 0.0;
        if(RadialMethod == HardSphere)
        {
            TSE = CalculateThroughSpace(PassedSite);
            double Test = PassedSite->Type->BindingEnergy[PassedSite->BondSize][0][ind];
            if(TSE > Test)
            {
                //PassedSite->BindingEnergy = 0.0;
                return -0.01;   // repulsive
            }
        }
        else if(RadialMethod == AtomAtomMMFF)
        {
            TSE = CalculateThroughSpace(PassedSite);
            return TSE;
        }
        double theEnergies[2];
        double InteractionEnergy;
        CharString Option;

        switch (theModel)
        {
            case BOC:

                BondOrder(PassedSite, theEnergies);
                InteractionEnergy = theEnergies[1] - theEnergies[0] - TSE;
                return InteractionEnergy;

            default:
                return 0;
        }
    }
    else
    {
        return 0.0;
    }
}


static void BondOrder(GridSite * ThisAtom, double theEnergies[])
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    
    Component *ThisType = ThisAtom->GetType();

    if(ThisType == theNULLSpecies)
        return;                 // site not occupied
    Eshort AdsorbateCoordination = ThisAtom->BondSize;

    if(ThisType->OnSurface == false) // DWDWDW
    {
        theEnergies[1] = ThisType->BindingEnergy[AdsorbateCoordination][0][ind];
        theEnergies[0] = theEnergies[1];
        return;
    }
    if(theInteractions->GetScale() < 1.e-8)
    {
        if(theGrid->Alloy_Type == theNULLSpecies)
        {
            theEnergies[1] = ThisType->BindingEnergy[AdsorbateCoordination][0][ind];
            theEnergies[0] = theEnergies[1];
            return;
        }
    }
    theEnergies[0] = 0.0;
    theEnergies[1] = 0.0;
    Eshort NeighborCoordination[12];
    Component *NeighborArray[12];
    double PairEnergy[12];
    GridSite *MetalAtom = NULL;
    for (int i = 0; i < 12; i++)
    {
        NeighborCoordination[i] = 0;
        NeighborArray[i] = NULL;
        PairEnergy[i] = 0.0;
    }
    register GridSite *Neighbor = NULL;

    // SECTION FOR RADIAL CONTRIBUTION METHOD

    double MyQ = ThisType->BindingEnergy[AdsorbateCoordination][0][ind];

    if(MyQ < 1.e-8 && MyQ > -1.e-8)
    {
        return;
    }
    double BindingEnergy = 0.0;
    Eshort NumberOfPairs = 0;

    double ZeroBindingEnergy = 0.0;

    for(int i = 0; i < ThisAtom->BondSize; ++i)
    {
        MetalAtom = ThisAtom->Bonds[i];
        PairEnergy[NumberOfPairs] = theInteractions->Get(ThisType, MetalAtom->Type, SC);
        ++NumberOfPairs;
        ZeroBindingEnergy += BondOrderConservation(PairEnergy, NumberOfPairs, AdsorbateCoordination, NeighborCoordination, ThisType, NeighborArray, MetalAtom->Type);
        for(int j = 0; j < MetalAtom->BondSize; ++j)
        {
            if(MetalAtom->Bonds[j]->Type != theNULLSpecies)
            {
                Neighbor = MetalAtom->Bonds[j]; // i think its faster this way
                if(Neighbor != ThisAtom && Neighbor->Type->GetName() != "H")
                {
                    NeighborCoordination[NumberOfPairs] = Neighbor->BondSize;
                    NeighborArray[NumberOfPairs] = Neighbor->Type;
                    PairEnergy[NumberOfPairs] = theInteractions->Get(Neighbor->Type, MetalAtom->Type, 0);
                    ++NumberOfPairs;
                }
            }
        }
        BindingEnergy += BondOrderConservation(PairEnergy, NumberOfPairs,AdsorbateCoordination, NeighborCoordination, ThisType, NeighborArray, MetalAtom->Type);
        NumberOfPairs = 0;
    }
    
    theEnergies[1] = BindingEnergy;
    theEnergies[0] = ZeroBindingEnergy;
    
    return;
}


static double BondOrderConservation(const double PairEnergy[],
    Eshort & NumberOfPairs, const Eshort & AdsorbateCoordination,
    const Eshort NeighborCoordination[], const Component * ThisType,
    /*const*/ Component * NeighType[], const Component * MetalType)
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    int MetalCoord = 1 + NumberOfPairs;
    double MyQ;
    Eshort AM = 0;

    if(MetalType != theGrid->Alloy_Type)
    {
        AM = 0;
    }
    else
    {
        AM = 1;
    }
    MyQ = ThisType->BindingEnergy[AdsorbateCoordination][AM][ind];
    double NeighborQ[10];
    if(MyQ)// > 1.e-8)
    {                           // user defined, do nothing
        MyQ = MyQ / (double)AdsorbateCoordination;
    }
    else if(ThisType->Atom)
    {
        MyQ = BOC_Function(PairEnergy[0], AdsorbateCoordination);
    }
    else
    {
        MyQ = BOC_Molecule_Function(ThisType, MetalType, AdsorbateCoordination);
    }

    for(Eshort i = 1; i < NumberOfPairs; ++i)
    {
        NeighborQ[i] = NeighType[i]->BindingEnergy[NeighborCoordination[i]][AM][ind];
        if(NeighborQ[i] > 0.00001)
        {                       // user defined
            NeighborQ[i] = NeighborQ[i] / (double)NeighborCoordination[i];
        }
        else if(NeighType[i]->Atom)
        {
            NeighborQ[i] =
                BOC_Function(PairEnergy[i], NeighborCoordination[i]);
        }
        else
        {
            NeighborQ[i] = BOC_Molecule_Function(NeighType[i], MetalType, NeighborCoordination[i]);
        }
    }
    if(NumberOfPairs == 1)
    {
        return MyQ;
    }
    double Correction = theInteractions->GetScale();
    double E = 0.0;
    double Fraction = 0.0;
    double TotalQ = 1/MyQ;

    for(Eshort j = 1; j < NumberOfPairs; ++j)
    {
        TotalQ += 1/NeighborQ[j];
    }
    Fraction = 1 + (1-MetalCoord)*(1/MyQ/TotalQ); //MyQ / (MyQ + TotalQ);
    Correction = theInteractions->GetScale();
    E = MyQ * (2 * Fraction - Fraction * Fraction);
//    if (theSimulation->Univ == 7)
//      cout << ThisType->SpeciesName << " " << MyQ << " " << E << endl;
    return (MyQ * (1.0 - Correction) + Correction * E);
}

static double BOC_Function(const double &E, const Eshort & Coordination)
{

    return (E / Coordination) * (2.0 - 1.0 / (Coordination));

}

static double BOC_Molecule_Function(const Component * ThisType, const Component * MetalType, const Eshort Coord)
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    double DissociationEnergy = ThisType->GetAtomizationEnergy();
    Component *Weak = ThisType->Weak;
    Component *Strong = ThisType->Strong;
    if(ThisType->GetOrientation() == Parallel)
    {
        DissociationEnergy = DissociationEnergy - Weak->GetAtomizationEnergy() - Strong->GetAtomizationEnergy();
        
        double StrongQ = theInteractions->Get(Strong->GetBindingAtom(), MetalType, SC);
        
        // Closed shell means weak binding
        if(!ThisType->GetClosedShell())
        {
            StrongQ = StrongQ * (2.0 - 1.0 / (double)Coord);
        }
        double Energy = (4.5 * StrongQ * StrongQ) / (3.0 * StrongQ + 8.0 * DissociationEnergy);
        
        return Energy / double (Coord);
    }
    else if(ThisType->GetOrientation() == Perpendicular)
    {
        // DissociationEnergy is set to D_{AB}.
        DissociationEnergy = DissociationEnergy - Weak->GetAtomizationEnergy();
        double StrongQ = theInteractions->Get(Strong->GetBindingAtom(), MetalType, SC);
        double Energy = 0.0;
        // NOTE: closed shell means weak binding...
        if(!ThisType->GetClosedShell())
        {
            Energy = StrongQ * (2.0 - 1.0 / double (Coord));
            Energy = Energy * Energy / (Energy + DissociationEnergy);
            return Energy / double (Coord);
        }
        else
        {
            Energy = StrongQ * StrongQ / ((StrongQ / (double)Coord) + DissociationEnergy);
            return Energy / double (Coord);
        }
    }
    else
    {
        return 0.0;
    }
}

double ModelError(CharString & FileName, double val)
{
    double SC = theSimulation->SurfaceCharge;
    int ind = theSimulation->FindChargeIndex(SC);
    ifstream patin(FileName.Ec_str());
    int Coord;
    int SType;
    double Error = 0;
    double Energy;
    double Desired;
    GridSite *Csite;
    while(patin >> Coord)
    {
        Csite = theGrid->LocateCentralSite(Coord);
        for(int i = 0; i < Csite->SurroundSize; ++i)
        {
            patin >> SType;
            Csite->Surround[i]->Type = Listing[SType];
        }
        Energy = theModel->CalculateBindingEnergy(Csite);
        patin >> Desired;
        Error += pow((Energy - Desired) * (Energy - Desired), 0.5);
    }
    return Error;
}

#include "Simplex.h"
// void Model::FindOptimalParameters()
// {

//     double *Values[3];

//     Values[0] = Dielectric;
//     Values[1] = theInteractions->GetScalePointer();
//     Values[2] = DistDepend;
//     double Limits[6];

//     Limits[0] = 1.0;
//     Limits[1] = 0.0;
//     Limits[2] = 1.0;
//     Limits[3] = 10.0;
//     Limits[4] = 1.0;
//     Limits[5] = 2.0;
//     double Noise = 0.0;
//     double ftol = 1.e-10;
//     CharString InputFile = "Patterns.dat";
//     CharString OutputFile = "Simplex.out";

//     // Amoeba(Values, 3, OutputFile, InputFile, ModelError, Limits, Noise, ftol);
    
//     cout << "\nThe dielectric constant is : " << *Dielectric;
//     cout << "\nThe scaling parameter is : " << theInteractions->GetScale();
//     cout << "\nThe DistDepend Parameter is : " << *DistDepend;
//     cout << endl;
//     return;
// }


void Model::InitializeThroughSpaceMatrix()
{ 
    for(int i = 0; i < 89; ++i)
    {
        for(int iO = 0; iO < 12; ++iO)
        {
            for(int jO = 0; jO < 12; ++jO)
            {
                for(int kO = 0; kO < 4; ++kO)
                {
                    for (int x = 0; x < 10; ++x)
                    {
                        ThroughSpaceInteraction[i][iO][jO][kO][x].Resize(NumberOfSpecies, NumberOfSpecies);
                        for(int ii = 0; ii < NumberOfSpecies; ++ii)
                        {
                            for(int jj = 0; jj < NumberOfSpecies; ++jj)
                            {
                                ThroughSpaceInteraction[i][iO][jO][kO][x][ii][jj] = 1.e20;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return;
}

double Model::CalculateThroughSpaceInteraction(GridSite * Site, GridSite * NeighborSite)
{
    double SC = theSimulation->SurfaceCharge;
    double Interaction = 0.0;
    assert(Site->Type != theNULLSpecies);
    assert(NeighborSite->Type != theNULLSpecies);
    double Origin[3] = { 0. };
    double NeighborOrigin[3] = { 0. };
    int Coordination = Site->BondSize;
    int NeighType = NeighborSite->Type->_idhook;
    int MyType = Site->Type->_idhook;

    /* Temporarily remove the local environment */
    Component *CurrentEnvironment[PatternDim];
    int ii;
    for(ii = Site->SurroundSize; ii--;)
    {
        CurrentEnvironment[ii] = Site->Surround[ii]->Type;
        if(Site->Surround[ii] != Site && Site->Surround[ii] != NeighborSite)
        {
            Site->Surround[ii]->Type = theNULLSpecies;
        }
    }
    
    double BOC_InteractionEnergy = 0.0;
    
    // double theEnergies[2];
    // BondOrder(Site, theEnergies);
    // BOC_InteractionEnergy = theEnergies[0] - theEnergies[1];
    
    /* Identify the position (origin) of center site and its neighbor */
    EPosition < double >Location = Site->GetPosition();
    EPosition < double >NeighLocation = NeighborSite->GetPosition();
    
    NeighborOrigin[0] = NeighLocation.GetXDistance();
    NeighborOrigin[1] = NeighLocation.GetYDistance();
    NeighborOrigin[2] = NeighLocation.GetZDistance();
    Origin[0] = Location.GetXDistance();
    Origin[1] = Location.GetYDistance();
    Origin[2] = Location.GetZDistance();
    
    /* Reorient the origin to account for periodic boundary
    * conditions, see GridSite.c */ 
    double Delta[3];
    
    Delta[0] = NeighborOrigin[0] - Origin[0];
    Delta[1] = NeighborOrigin[1] - Origin[1];
    Delta[2] = NeighborOrigin[2] - Origin[2];
    double Cmp = 0.5 * theGrid->GridLength[0];
    
    if(Delta[0] < -Cmp)
    {
        NeighborOrigin[0] = theGrid->GridLength[0] + NeighborOrigin[0];
    }
    else if(Delta[0] > Cmp)
    {
        NeighborOrigin[0] = -theGrid->GridLength[0] + NeighborOrigin[0];
    }
    Cmp = 0.5 * theGrid->GridLength[1];
    if(Delta[1] < -Cmp)
    {
        NeighborOrigin[1] = theGrid->GridLength[1] + NeighborOrigin[1];
    }
    else if(Delta[1] > Cmp)
    {
        NeighborOrigin[1] = -theGrid->GridLength[1] + NeighborOrigin[1];
    }
    
    /* Retrieve Species Coordinates for Center Site */
    
    Eshort SiteID = Listing[MyType]->_id;
    
    theSpecies->theGeometries->theCoordinates.ResetToFront();
    Coordinates SiteCoords;

    bool Finish = false; // DWDWDW  // DWDWDW
    
    while(theSpecies->theGeometries->theCoordinates && !Finish)
    {
        if (theSpecies->theGeometries->theCoordinates().SurfaceCharge == SC)
        {
          if(SiteID == theSpecies->theGeometries->theCoordinates().SpeciesName
              && Site->BondSize == theSpecies->theGeometries->theCoordinates().MetalCoord)
          {
              SiteCoords = theSpecies->theGeometries->theCoordinates();
              Finish = true; // DWDWDW
          }
        }
        ++theSpecies->theGeometries->theCoordinates;
    }
    
    /* Retrieving Species Coordinates for Surrounding Site */
    
    SiteID = Listing[NeighType]->_id;
    theSpecies->theGeometries->theCoordinates.ResetToFront();
    Coordinates NeighCoords;
    bool Finish2 = false; // DWDWDW  // DWDWDW
     
    while(theSpecies->theGeometries->theCoordinates && !Finish2)
    {
        if (theSpecies->theGeometries->theCoordinates().SurfaceCharge == SC)
        {
          if(SiteID == theSpecies->theGeometries->theCoordinates().SpeciesName
              && NeighborSite->BondSize ==
              theSpecies->theGeometries->theCoordinates().MetalCoord)
          {
              NeighCoords = theSpecies->theGeometries->theCoordinates();
              Finish2 = true; // DWDWDW
          }
        }
        ++theSpecies->theGeometries->theCoordinates;
    }
    
    /* If atomic Geometries for these species have not been defined,
    * set interaction to zero and return */ 
    if(!(Finish && Finish2))
    {
        for(ii = Site->SurroundSize ; ii-- ; /*none*/ )
        {
            Site->Surround[ii]->Type = CurrentEnvironment[ii];
        }
        return 1.e10;
    }
    
    Eshort MyMMFF;
    Eshort NeighMMFF;
    double MyCharge, NeighCharge;
    
    double LocPointer[3] = { 0. };      // Exact location of atom 1 
    double NeighLocPointer[3] = { 0. }; // Exact Location of atom 2
    double Distance = 0.0;       // Distance between atom 1 and 2
    double x, y, z;
    
//    double SiteA_Angle = Site->Angle;
//    double AdsA_Angle = 0.0;
//    double A_Net_Angle = 0.0;
//    double MetalZPosition = Site->Bonds[0]->Location.GetZDistance();
//    double GridZPosition = Site->Location.GetZDistance();
    int MyOrientation = Site->Orientation;
    int MyOrientationPos = Site->OrientationPos;
    int NeighOrientation = NeighborSite->Orientation;
    int NeighOrientationPos = NeighborSite->OrientationPos;
    
//    if(Coordination == 2)
//    {
//        // rotation by 180 degrees for bridge site 
//        AdsA_Angle = double (MyOrientation) * 180.0;     
//    }
//    else if(Coordination == 3)
//    {
//        AdsA_Angle = double (MyOrientation) * 120.0;
//    }
//    else if(Coordination == 4)
//    {
//        AdsA_Angle = double (MyOrientation) * 90.0;
//    }
//    else if(Coordination == 6)
//    {
//        AdsA_Angle = double (MyOrientation) * 60.0;
//    }
//    else
//    {
//        AdsA_Angle = double (MyOrientation) * 30.0;
//    }
//    
//    A_Net_Angle = (SiteA_Angle + AdsA_Angle);
//    while(A_Net_Angle >= 360.0)
//    {
//        A_Net_Angle = A_Net_Angle - 360.0;
//    }
//    while(A_Net_Angle < 0.0)
//    {
//        // place angle between 0 and 360
//        A_Net_Angle = 360. + A_Net_Angle;
//    }
    
//    double NeighMetalZPosition = NeighborSite->Bonds[0]->Location.GetZDistance();
    
//    double SiteB_Angle = NeighborSite->Angle;
//    double AdsB_Angle = 0.0;
//    double B_Net_Angle = 0.0;
//    int NeighCoordination = NeighborSite->BondSize;
//    int NeighOrientation = NeighborSite->Orientation;
    
//    if(NeighCoordination == 2)
//    {
//        // rotation by 180 degrees for bridge site
//        AdsB_Angle = (double)NeighOrientation *180.0;
//    }
//    else if(NeighCoordination == 3)
//    {
//        // rotation by 120 degrees for hollow site
//        AdsB_Angle = (double)NeighOrientation *120.0;
//    }
//    else if(NeighCoordination == 4)
//    {
//        // rotation by 90 degrees for hollow site for 100 surface
//        AdsB_Angle = (double)NeighOrientation *90.0;
//    }
//    else if(NeighCoordination == 6)
//    {
//        // rotation by 90 degrees for hollow site for 100 surface
//        AdsB_Angle = (double)NeighOrientation *60.0;
//    }
//    else
//    {
//        AdsB_Angle = (double)NeighOrientation *30.0;
//    }
    
//    B_Net_Angle = (SiteB_Angle + AdsB_Angle);
//    while(B_Net_Angle >= 360.)
//    {
//        B_Net_Angle = B_Net_Angle - 360.0;
//    }
//    while(B_Net_Angle < 0.)
//    {
//        // place angle between 0 and 360
//        B_Net_Angle = 360. + B_Net_Angle;
//    }
    Site->NumAtoms = SiteCoords.NumAtoms;
    /* Looping through atoms of the adsorbate in the center site */
    for(int iNum = 0; iNum < SiteCoords.NumAtoms; ++iNum)
    {
//        x = SiteCoords.Position[iNum].GetXDistance();
//        y = SiteCoords.Position[iNum].GetYDistance();
//        z = SiteCoords.Position[iNum].GetZDistance();
        MyMMFF = SiteCoords.MMFF94_Type[iNum];  // Merck model atom type
        MyCharge = SiteCoords.Charge[iNum];     // Charge on the atom 
//        double newarctangent = atan2(x, y) + (A_Net_Angle) * (M_PI / 180.);
//        double measure = pow((x * x + y * y), 0.5);
//        LocPointer[0] = Origin[0] + sin(newarctangent) * measure;
//        LocPointer[1] = Origin[1] + cos(newarctangent) * measure;
//        //LocPointer[2] = z + MetalZPosition;
//        LocPointer[2] = z + GridZPosition;
//        if (TempCoordLearn == 1)
//        {
//          cout << "Hi" << endl;
//          TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][0] = LocPointer[0] - Origin[0];
//          TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][1] = LocPointer[1] - Origin[1];
//          TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][2] = LocPointer[2] - Origin[2];
//        }
//        if (theSimulation->Univ == 1)
//          cout << SiteCoords.AtomName[iNum] << endl;
//        Site->Ads_Atom[iNum] = SiteCoords.AtomName[iNum];
        /* Looping through the atoms of the adsorbate in the
        * neighboring site */ 
        for(int jNum = 0; jNum < NeighCoords.NumAtoms; ++jNum)
        {
//            x = NeighCoords.Position[jNum].GetXDistance();
//            y = NeighCoords.Position[jNum].GetYDistance();
//            z = NeighCoords.Position[jNum].GetZDistance();
            NeighMMFF = NeighCoords.MMFF94_Type[jNum];
            NeighCharge = NeighCoords.Charge[jNum];

//            double newarctangent =
//                atan2(x, y) + (B_Net_Angle) * (M_PI / 180.);
//
//            double measure = pow((x * x + y * y), 0.5);

            for (int i = 0; i < 3; ++i)
            {
              NeighLocPointer[i] = NeighborOrigin[i] + TempCoord[NeighborSite->BondSize-1][NeighborSite->GetType()->_idhook][NeighOrientationPos][NeighOrientation][jNum][i];
              LocPointer[i] = Origin[i] + TempCoord[Site->BondSize-1][Site->GetType()->_idhook][MyOrientationPos][MyOrientation][iNum][i];
            }
            
//            NeighLocPointer[0] =
//                NeighborOrigin[0] + sin(newarctangent) * measure;
//            NeighLocPointer[1] =
//                NeighborOrigin[1] + cos(newarctangent) * measure;
//            //NeighLocPointer[2] = z + NeighMetalZPosition;
//            NeighLocPointer[2] = z + GridZPosition;
            
            x = LocPointer[0] - NeighLocPointer[0];
            y = LocPointer[1] - NeighLocPointer[1];
            z = LocPointer[2] - NeighLocPointer[2];
            Distance = pow((x * x + y * y + z * z), 0.5);
            double dum;
            dum = GetMMFFEnergy(MyMMFF, NeighMMFF, MyCharge, NeighCharge,
                    Distance);
//            if (theSimulation->Univ == 6)
//              cout << dum << endl;
            Interaction += dum;
//            if (theSimulation->Univ == 1)
//            {
//              cout << SiteCoords.AtomName[iNum] << " " << NeighCoords.AtomName[jNum] << " " << x << " " << y << " " << z << " " << Distance << endl; 
//            }
            
        }
    }
//    if (theSimulation->Univ == 1)
//      cout << Interaction << endl;
    /* Reinstate the local environment */
    for(ii = Site->SurroundSize; ii--;)
    {
        Site->Surround[ii]->Type = CurrentEnvironment[ii];
    }

    Interaction = BOC_InteractionEnergy + Interaction;
    if(Interaction > 1.e18)
    {
        Interaction = 1.e18;    // been here, done that    
    }
    
    return Interaction;
    
    
//    double SC = theSimulation->SurfaceCharge;
//    double Interaction = 0.0;
//    assert(Site->Type != theNULLSpecies);
//    assert(NeighborSite->Type != theNULLSpecies);
//
//    double Origin[3] = { 0. };
//    double NeighborOrigin[3] = { 0. };
//    int Coordination = Site->BondSize;
//    int NeighType = NeighborSite->Type->_idhook;
//    int MyType = Site->Type->_idhook;
//
//    /* Temporarily remove the local environment */
//    Component *CurrentEnvironment[90];
//    int ii;
//    for(ii = Site->SurroundSize; ii--;)
//    {
//        CurrentEnvironment[ii] = Site->Surround[ii]->Type;
//        if(Site->Surround[ii] != Site && Site->Surround[ii] != NeighborSite)
//        {
//            Site->Surround[ii]->Type = theNULLSpecies;
//        }
//    }
//    
//    double BOC_InteractionEnergy = 0.0;
//    
//    // double theEnergies[2];
//    // BondOrder(Site, theEnergies);
//    // BOC_InteractionEnergy = theEnergies[0] - theEnergies[1];
//    
//    /* Identify the position (origin) of center site and its neighbor */
//    EPosition < double >Location = Site->GetPosition();
//    EPosition < double >NeighLocation = NeighborSite->GetPosition();
//    
//    NeighborOrigin[0] = NeighLocation.GetXDistance();
//    NeighborOrigin[1] = NeighLocation.GetYDistance();
//    NeighborOrigin[2] = NeighLocation.GetZDistance();
//    Origin[0] = Location.GetXDistance();
//    Origin[1] = Location.GetYDistance();
//    Origin[2] = Location.GetZDistance();
//    
//    /* Reorient the origin to account for periodic boundary
//    * conditions, see GridSite.c */ 
//    double Delta[3];
//    
//    Delta[0] = NeighborOrigin[0] - Origin[0];
//    Delta[1] = NeighborOrigin[1] - Origin[1];
//    Delta[2] = NeighborOrigin[2] - Origin[2];
//    double Cmp = 0.5 * theGrid->GridLength[0];
//    
//    if(Delta[0] < -Cmp)
//    {
//        NeighborOrigin[0] = theGrid->GridLength[0] + NeighborOrigin[0];
//    }
//    else if(Delta[0] > Cmp)
//    {
//        NeighborOrigin[0] = -theGrid->GridLength[0] + NeighborOrigin[0];
//    }
//    Cmp = 0.5 * theGrid->GridLength[1];
//    if(Delta[1] < -Cmp)
//    {
//        NeighborOrigin[1] = theGrid->GridLength[1] + NeighborOrigin[1];
//    }
//    else if(Delta[1] > Cmp)
//    {
//        NeighborOrigin[1] = -theGrid->GridLength[1] + NeighborOrigin[1];
//    }
//    
//    /* Retrieve Species Coordinates for Center Site */
//    
//    Eshort SiteID = Listing[MyType]->_id;
//    
//    theSpecies->theGeometries->theCoordinates.ResetToFront();
//    Coordinates SiteCoords;
//
//    bool Finish = false; // DWDWDW  // DWDWDW
//    
//    while(theSpecies->theGeometries->theCoordinates && !Finish)
//    {
//        if (theSpecies->theGeometries->theCoordinates().SurfaceCharge == SC)
//        {
//          if(SiteID == theSpecies->theGeometries->theCoordinates().SpeciesName
//              && Site->BondSize == theSpecies->theGeometries->theCoordinates().MetalCoord)
//          {
//              SiteCoords = theSpecies->theGeometries->theCoordinates();
//              Finish = true; // DWDWDW
//          }
//        }
//        ++theSpecies->theGeometries->theCoordinates;
//    }
//    
//    /* Retrieving Species Coordinates for Surrounding Site */
//    
//    SiteID = Listing[NeighType]->_id;
//    theSpecies->theGeometries->theCoordinates.ResetToFront();
//    Coordinates NeighCoords;
//    bool Finish2 = false; // DWDWDW  // DWDWDW
//     
//    while(theSpecies->theGeometries->theCoordinates && !Finish2)
//    {
//        if (theSpecies->theGeometries->theCoordinates().SurfaceCharge == SC)
//        {
//          if(SiteID == theSpecies->theGeometries->theCoordinates().SpeciesName
//              && NeighborSite->BondSize ==
//              theSpecies->theGeometries->theCoordinates().MetalCoord)
//          {
//              NeighCoords = theSpecies->theGeometries->theCoordinates();
//              Finish2 = true; // DWDWDW
//          }
//        }
//        ++theSpecies->theGeometries->theCoordinates;
//    }
//    
//    /* If atomic Geometries for these species have not been defined,
//    * set interaction to zero and return */ 
//    if(!(Finish && Finish2))
//    {
//        for(ii = Site->SurroundSize ; ii-- ; /*none*/ )
//        {
//            Site->Surround[ii]->Type = CurrentEnvironment[ii];
//        }
//        return 1.e10;
//    }
//    
//    Eshort MyMMFF;
//    Eshort NeighMMFF;
//    double MyCharge, NeighCharge;
//    
//    double LocPointer[3] = { 0. };      // Exact location of atom 1 
//    double NeighLocPointer[3] = { 0. }; // Exact Location of atom 2
//    double Distance = 0.0;       // Distance between atom 1 and 2
//    double x, y, z;
//    
//    double SiteA_Angle = Site->Angle;
//    double AdsA_Angle = 0.0;
//    double A_Net_Angle = 0.0;
//    double MetalZPosition = Site->Bonds[0]->Location.GetZDistance();
//    int MyOrientation = Site->Orientation;
//    
//    if(Coordination == 2)
//    {
//        // rotation by 180 degrees for bridge site 
//        AdsA_Angle = double (MyOrientation) * 180.0;     
//    }
//    else if(Coordination == 3)
//    {
//        AdsA_Angle = double (MyOrientation) * 120.0;
//    }
//    else if(Coordination == 4)
//    {
//        AdsA_Angle = double (MyOrientation) * 90.0;
//    }
//    else if(Coordination == 6)
//    {
//        AdsA_Angle = double (MyOrientation) * 60.0;
//    }
//    else
//    {
//        AdsA_Angle = double (MyOrientation) * 30.0;
//    }
//    
//    A_Net_Angle = (SiteA_Angle + AdsA_Angle);
//    while(A_Net_Angle >= 360.0)
//    {
//        A_Net_Angle = A_Net_Angle - 360.0;
//    }
//    while(A_Net_Angle < 0.0)
//    {
//        // place angle between 0 and 360
//        A_Net_Angle = 360. + A_Net_Angle;
//    }
//    
//    double NeighMetalZPosition = NeighborSite->Bonds[0]->Location.GetZDistance();
//    double SiteB_Angle = NeighborSite->Angle;
//    double AdsB_Angle = 0.0;
//    double B_Net_Angle = 0.0;
//    int NeighCoordination = NeighborSite->BondSize;
//    int NeighOrientation = NeighborSite->Orientation;
//    
//    if(NeighCoordination == 2)
//    {
//        // rotation by 180 degrees for bridge site
//        AdsB_Angle = (double)NeighOrientation *180.0;
//    }
//    else if(NeighCoordination == 3)
//    {
//        // rotation by 120 degrees for hollow site
//        AdsB_Angle = (double)NeighOrientation *120.0;
//    }
//    else if(NeighCoordination == 4)
//    {
//        // rotation by 90 degrees for hollow site for 100 surface
//        AdsB_Angle = (double)NeighOrientation *90.0;
//    }
//    else if(NeighCoordination == 6)
//    {
//        // rotation by 90 degrees for hollow site for 100 surface
//        AdsB_Angle = (double)NeighOrientation *60.0;
//    }
//    else
//    {
//        AdsB_Angle = (double)NeighOrientation *30.0;
//    }
//    
//    B_Net_Angle = (SiteB_Angle + AdsB_Angle);
//    while(B_Net_Angle >= 360.)
//    {
//        B_Net_Angle = B_Net_Angle - 360.0;
//    }
//    while(B_Net_Angle < 0.)
//    {
//        // place angle between 0 and 360
//        B_Net_Angle = 360. + B_Net_Angle;
//    }
//    
//    /* Looping through atoms of the adsorbate in the center site */
//    for(int iNum = 0; iNum < SiteCoords.NumAtoms; ++iNum)
//    {
//        x = SiteCoords.Position[iNum].GetXDistance();
//        y = SiteCoords.Position[iNum].GetYDistance();
//        z = SiteCoords.Position[iNum].GetZDistance();
//        MyMMFF = SiteCoords.MMFF94_Type[iNum];  // Merck model atom type
//        MyCharge = SiteCoords.Charge[iNum];     // Charge on the atom 
//        double newarctangent = atan2(x, y) + (A_Net_Angle) * (M_PI / 180.);
//        double measure = pow((x * x + y * y), 0.5);
//
//        LocPointer[0] = Origin[0] + sin(newarctangent) * measure;
//        LocPointer[1] = Origin[1] + cos(newarctangent) * measure;
//        LocPointer[2] = z + MetalZPosition;
//        
//        /* Looping through the atoms of the adsorbate in the
//        * neighboring site */ 
//        for(int jNum = 0; jNum < NeighCoords.NumAtoms; ++jNum)
//        {
//            x = NeighCoords.Position[jNum].GetXDistance();
//            y = NeighCoords.Position[jNum].GetYDistance();
//            z = NeighCoords.Position[jNum].GetZDistance();
//            NeighMMFF = NeighCoords.MMFF94_Type[jNum];
//            NeighCharge = NeighCoords.Charge[jNum];
//
//            double newarctangent =
//                atan2(x, y) + (B_Net_Angle) * (M_PI / 180.);
//
//            double measure = pow((x * x + y * y), 0.5);
//
//            NeighLocPointer[0] =
//                NeighborOrigin[0] + sin(newarctangent) * measure;
//            NeighLocPointer[1] =
//                NeighborOrigin[1] + cos(newarctangent) * measure;
//            NeighLocPointer[2] = z + NeighMetalZPosition;
//            
//            x = LocPointer[0] - NeighLocPointer[0];
//            y = LocPointer[1] - NeighLocPointer[1];
//            z = LocPointer[2] - NeighLocPointer[2];
//            Distance = pow((x * x + y * y + z * z), 0.5);
//
//            Interaction +=
//                GetMMFFEnergy(MyMMFF, NeighMMFF, MyCharge, NeighCharge,
//                    Distance);
//        }
//    }
//    /* Reinstate the local environment */
//    for(ii = Site->SurroundSize; ii--;)
//    {
//        Site->Surround[ii]->Type = CurrentEnvironment[ii];
//    }
//
//    Interaction = BOC_InteractionEnergy + Interaction;
//    if(Interaction > 1.e18)
//    {
//        Interaction = 1.e18;    // been here, done that    
//    }
//    
//    return Interaction;
}


double Model::CalculateThroughSpaceFixedOrientation(GridSite * PassedSite)
{
    Component *MyType = PassedSite->GetType();
    double SC = theSimulation->SurfaceCharge;
    
    int ind = theSimulation->FindChargeIndex(SC);
    if(MyType == theNULLSpecies)
    {
        return 0;
    }
    if(RadialMethod == NoCoreRepulsion)
        return 0.0;
    double Energy = 0.0;
    int Myid = MyType->_idhook;
    register GridSite *const *Surround = PassedSite->Surround;
    int MyCoord = PassedSite->BondSize;
    double BaseEnergy = MyType->BindingEnergy[MyCoord][0][ind];
    GridSite *Neighbor;
    int MyCoordIndex = MyCoord - 1;

    // assert(MyCoordIndex != 3);
    InteractionAddress = NULL;
    double InteractionEnergy;
    assert(theNULLSpecies->_id == 0);
    for(int i = PassedSite->SurroundSize; i--;)
    {
        if(Surround[i]->Type->_id)
        {
            // this might be faster, but assumes _id is 0 for theNULLSpecies
            if(Surround[i] != PassedSite)
            {
                Neighbor = Surround[i];
                int Neighid = Neighbor->Type->_idhook;
                int NeOr = Neighbor->Orientation;
                int NeighCoordination = Surround[i]->BondSize;
                int NeighIndex;
                
                if(NeighCoordination == 2)
                {
                    NeighIndex = 2 * Neighbor->OrientationPos + NeOr;
                    // 2 possible orientations for bridging adsorbate
                    // 3 possible orientations for bridging metal site
                }
                else if(NeighCoordination == 3)
                {
                    NeighIndex = 3 * Neighbor->OrientationPos + NeOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(NeighCoordination == 4)
                {
                    NeighIndex = 4 * Neighbor->OrientationPos + NeOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(NeighCoordination == 6)
                {
                    NeighIndex = 6 * Neighbor->OrientationPos + NeOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else
                {
                    NeighIndex = NeOr;  // rotation by 30.  degrees for atop 
                    // 12 possible orientations for atop adsorbate
                    // 1 possible orientations for atop metal site
                }
                
                int MyOr = PassedSite->Orientation;
                int MyIndex;
                
                if(MyCoord == 2)
                {
                    MyIndex = 2 * PassedSite->OrientationPos + MyOr;
                    // 2 possible orientations for bridging adsorbate
                    // 3 possible orientations for bridging metal site
                }
                else if(MyCoord == 3)
                {
                    MyIndex = 3 * PassedSite->OrientationPos + MyOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(MyCoord == 4)
                {
                    MyIndex = 4 * PassedSite->OrientationPos + MyOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(MyCoord == 6)
                {
                    MyIndex = 6 * PassedSite->OrientationPos + MyOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else
                {
                    MyIndex = MyOr;
                    // 12 possible orientations for bridging adsorbate
                    // 1 possible orientations for bridging metal site
                }
                InteractionEnergy = ThroughSpaceInteraction[i][NeighIndex][MyIndex]
                    [MyCoordIndex][ind][Myid][Neighid];
                InteractionAddress =
                    &ThroughSpaceInteraction[i][NeighIndex][MyIndex]
                    [MyCoordIndex][ind][Myid][Neighid];
                if(InteractionEnergy > 1.e19)
                {
                    InteractionEnergy = CalculateThroughSpaceInteraction(PassedSite, Neighbor);
                    ThroughSpaceInteraction[i][NeighIndex][MyIndex][MyCoordIndex][ind][Myid][Neighid] = InteractionEnergy;
                }
                if(RadialMethod == HardSphere)
                {
                    if(InteractionEnergy > BaseEnergy)
                        return 1.e20;
                }
                else
                {
                    Energy += InteractionEnergy;                    
                }
            }
        }
    }
    EPosition < double >Location = PassedSite->GetPosition();
      double Origin[3];
      Origin[0] = Location.GetXDistance();
      Origin[1] = Location.GetYDistance();
      Origin[2] = Location.GetZDistance();
      for (int no_at = 0; no_at < PassedSite->GetType()->NumAtoms; ++no_at)
      {
        for (int no = 0; no < 3; ++no)
        {
          PassedSite->Ads_Pos[no_at][no] = Origin[no] + TempCoord[PassedSite->BondSize-1][PassedSite->GetType()->_idhook][PassedSite->OrientationPos][PassedSite->Orientation][no_at][no];
        }
      }
    return Energy;
}

double Model::CalculateThroughSpace(GridSite * PassedSite)
{
    EPosition < double >Location = PassedSite->GetPosition();
    double Origin[3];
    Origin[0] = Location.GetXDistance();
    Origin[1] = Location.GetYDistance();
    Origin[2] = Location.GetZDistance();
    double SC = theSimulation->SurfaceCharge;
    Component *MyType = PassedSite->GetType();
    // DW 8/24/04: NeighSite unused, so commented out
    //register GridSite *NeighSite;

    if(MyType == theNULLSpecies)
    {
        return 0;
    }
    if(RadialMethod == NoCoreRepulsion)
    {
        return 0.0;
    }
    
    // DW 9/6/04: Eric only initialize 0-9,11 (skipped 10).
    // This occasionally gave weird results.
    // Specifically, the debug version from VC++ initializes
    // doubles to a very large negative number.
    double Energy[12];
    for(int jjj=0 ; jjj<12 ; jjj++)
    {
        Energy[jjj] = 0.0;
    }
    
    int Myid = MyType->_idhook;
    register GridSite *const *Surround = PassedSite->Surround;
    int MyCoord = PassedSite->BondSize;
    int OldOrient = PassedSite->Orientation;
    int ind = theSimulation->FindChargeIndex(SC);
    double BaseEnergy = MyType->BindingEnergy[MyCoord][0][ind];
    GridSite *Neighbor;
    double Low = 1.e5;
    int MyCoordIndex = MyCoord - 1;

    //      assert(MyCoordIndex != 3);
    double InteractionEnergy;

    // DW 8/24/04: LocalNullSpecies not used, so commented out
    //Component *LocalNULLSpecies = theNULLSpecies;
    int NumOr = NumberOfOrientations[MyCoordIndex];

    //int Eval;       
    //int ispeed = PassedSite[i]->Type->_id;
    //int i = 0;
    //while(!Eval && PassedSite[i]) Eval += PassedSite[i]->Type->_id; 
    //for(ispeed=i ; ispeed < Lasti+10; ++i)
    //    Eval[ispeed] += PassedSite[ispeed]->Type->_id;
    // if(RadialMethod == AtomAtomMMFF) {
    
    assert(theNULLSpecies->_id == 0);
    
    for(int i = PassedSite->SurroundSize; i--;)
    {
        if(Surround[i]->Type->_id)
        {
            // this might be faster, but assumes _id is 0 for theNULLSpecies
            if(Surround[i] != PassedSite)
            {
                Neighbor = Surround[i];
                int Neighid = Neighbor->Type->_idhook;
                int NeOr = Neighbor->Orientation;
                int NeighCoordination = Surround[i]->BondSize;
                int NeighIndex;

                if(NeighCoordination == 2)
                {
                    NeighIndex = 2 * Neighbor->OrientationPos + NeOr;
                    // 2 possible orientations for bridging adsorbate
                    // 3 possible orientations for bridging metal site
                }
                else if(NeighCoordination == 3)
                {
                    NeighIndex = 3 * Neighbor->OrientationPos + NeOr;
                    // 3 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(NeighCoordination == 4)
                {
                    NeighIndex = 4 * Neighbor->OrientationPos + NeOr;
                    // 4 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else if(NeighCoordination == 6)
                {
                    NeighIndex = 6 * Neighbor->OrientationPos + NeOr;
                    // 4 possible orientations for hollow adsorbate
                    // 2 possible orientations for hollow metal site
                }
                else
                {
                    NeighIndex = NeOr;  // rotation by 30.  degrees for atop 
                    // 12 possible orientations for bridging adsorbate
                    // 1 possible orientations for bridging metal site
                }
                
                int k;
                for(k = 0; k < NumOr; k++)
                {
                    int MyIndex;
                    Component *thisComp = PassedSite->GetType();
                    
                    if(MyCoord == 2)
                    {
                        MyIndex = 2 * PassedSite->OrientationPos + k;
                        // 2 possible orientations for bridging adsorbate
                        // 3 possible orientations for bridging metal site
                    }
                    else if(MyCoord == 3)
                    {
                        MyIndex = 3 * PassedSite->OrientationPos + k;
                        // 3 possible orientations for hollow adsorbate
                        // 2 possible orientations for hollow metal site
                    }
                    else if(MyCoord == 4)
                    {
                        MyIndex = 4 * PassedSite->OrientationPos + k;
                        // 4 possible orientations for hollow adsorbate
                        // 2 possible orientations for hollow metal site
                    }
                    else if(MyCoord == 6)
                    {
                        MyIndex = 6 * PassedSite->OrientationPos + k;
                    }
                    else
                    {
                        MyIndex = k;
                    }
                    
                    PassedSite->Orientation = k;
                    InteractionEnergy = ThroughSpaceInteraction[i][NeighIndex][MyIndex][MyCoordIndex][ind][Myid][Neighid];
                    if(InteractionEnergy > 1.e19)
                    {
                        // haven't been here, not done that
                        InteractionEnergy = CalculateThroughSpaceInteraction(PassedSite, Neighbor);
                        ThroughSpaceInteraction[i][NeighIndex][MyIndex][MyCoordIndex][ind][Myid][Neighid] = InteractionEnergy;
                    }
                    if(RadialMethod == HardSphere)
                    {
                        if(InteractionEnergy > BaseEnergy)
                        {
                            Energy[k] += 1.e20;
                        }
                        else
                        {
                            Energy[k] += InteractionEnergy;
                        }
                    }
                    else
                    {
                        Energy[k] += InteractionEnergy;
                    }

                }
            }
        }
    }
    // }
    int Orien;
    bool Unchanged = 1;
    for(int k = NumOr; k--;)
    {
        if(Energy[k] < Low)
        {
            Low = Energy[k];
            Orien = k;
            Unchanged = 0;
        }
    }
    if (Unchanged == 1)
    {
      Orien = OldOrient;
    }
      PassedSite->Orientation = Orien;
      for (int no_at = 0; no_at < PassedSite->GetType()->NumAtoms; ++no_at)
      {
        for (int no = 0; no < 3; ++no)
        {
          PassedSite->Ads_Pos[no_at][no] = Origin[no] + TempCoord[PassedSite->BondSize-1][PassedSite->GetType()->_idhook][PassedSite->OrientationPos][Orien][no_at][no];
        }
      }
    if(RadialMethod == HardSphere)
    {
        if(Low < 1.e5)
            return 0.0;
        else
            return 1.e20;
    }
    return Low;
}

double Model::GetMMFFEnergy(Eshort & MyMMFF, Eshort & NeighMMFF,
    double &MyCharge, double &NeighCharge, double &Distance)
{
    MMFF.ResetToFront();
    MMFF_Parameters Params;
    MMFF_Parameters MyParams;
    MMFF_Parameters NeighParams;

    bool FoundParameter[2] = { false, false }; // DWDWDW x3

    while(MMFF)
    {
        Params = MMFF.Get();
        if(Params.id == MyMMFF)
        {
            FoundParameter[0] = true; // DWDWDW
            MyParams = Params;
        }
        if(Params.id == NeighMMFF)
        {
            FoundParameter[1] = true; // DWDWDW
            NeighParams = Params;
        }
        ++MMFF;
    }

    if(!FoundParameter[0])
    {
        Eshort PType;
        char LineRead[256];
        ifstream Pin("MMFF94.dat");

        for(int i = 21; i--;)
            Pin.getline(LineRead, 255);
        while(Pin >> PType && PType != MyMMFF)
            Pin.getline(LineRead, 255);
        assert(PType == MyMMFF);
        MyParams.id = PType;
        Pin >> MyParams.alpha;
        Pin >> MyParams.N;
        Pin >> MyParams.A;
        Pin >> MyParams.G;
        Pin >> MyParams.DA;
        MyParams.AoverNsqrt = pow((MyParams.alpha / MyParams.N), 0.5);
        MyParams.A_alpha_factor = MyParams.A * pow(MyParams.alpha, 0.25);
        MMFF.Add(MyParams);
        Pin.close();
    }
    if(!FoundParameter[1])
    {
        Eshort PType;
        char LineRead[256];
        ifstream Pin("MMFF94.dat");

        for(int i = 21; i--;)
            Pin.getline(LineRead, 255);
        while(Pin >> PType && PType != NeighMMFF)
            Pin.getline(LineRead, 255);
        assert(PType == NeighMMFF);
        NeighParams.id = PType;
        Pin >> NeighParams.alpha;
        Pin >> NeighParams.N;
        Pin >> NeighParams.A;
        Pin >> NeighParams.G;
        Pin >> NeighParams.DA;
        NeighParams.AoverNsqrt =
            pow((NeighParams.alpha / NeighParams.N), 0.5);
        NeighParams.A_alpha_factor =
            NeighParams.A * pow(NeighParams.alpha, 0.25);
        MMFF.Add(NeighParams);
        Pin.close();
    }

    double Rstar;
    bool DA_pair = false; // DWDWDW  // DWDWDW

    if(MyMMFF == NeighMMFF)
    {
        Rstar = MyParams.A_alpha_factor;
    }
    else
    {
        double MyRstar = MyParams.A_alpha_factor;
        double NeighRstar = NeighParams.A_alpha_factor;

        if((MyParams.DA == 'D' && NeighParams.DA == 'A') ||
            (MyParams.DA == 'A' && NeighParams.DA == 'D'))
        {
            DA_pair = true; // DWDWDW
            Rstar = 0.8 * 0.5 * (MyRstar + NeighRstar);
        }
        else
        {
            double gamma = (MyRstar - NeighRstar) / (MyRstar + NeighRstar);

            Rstar = 0.5 * (MyRstar + NeighRstar)
                * (1. + 0.2 * (1. - exp(-12. * gamma * gamma)));
        }
    }
    double Eps;
    double SimpRstar = Rstar * Rstar * Rstar * Rstar * Rstar * Rstar;

    if(DA_pair)
    {
        double simplifyRstar = 1.25 * Rstar;

        simplifyRstar =
            simplifyRstar * simplifyRstar * simplifyRstar * simplifyRstar *
            simplifyRstar * simplifyRstar;
        Eps =
            90.58 * MyParams.G * NeighParams.G * MyParams.alpha *
            NeighParams.alpha / ((MyParams.AoverNsqrt +
                                     NeighParams.AoverNsqrt) * simplifyRstar);
    }
    else
    {
        Eps =
            181.16 * MyParams.G * NeighParams.G * MyParams.alpha *
            NeighParams.alpha / ((MyParams.AoverNsqrt +
                                     NeighParams.AoverNsqrt) * SimpRstar);
    }

    double Factor1 = 1.07 * Rstar / (Distance + 0.07 * Rstar);
    double Rstar7 = SimpRstar * Rstar;
    double Distance7 =
        Distance * Distance * Distance * Distance * Distance * Distance *
        Distance;
    double Factor1_7 =
        Factor1 * Factor1 * Factor1 * Factor1 * Factor1 * Factor1 * Factor1;
    double Factor2 = 1.12 * Rstar7 / (Distance7 + 0.12 * Rstar7);
    double EVDW = Eps * Factor1_7 * (Factor2 - 2.)*4.184;
    double aaaa = Distance + 0.05;
    double MMFF_Energy = 332.0716 * MyCharge * NeighCharge / (*Dielectric * pow(aaaa, *DistDepend))*4.184; // Dielectric constant = 10
    EVDW += MMFF_Energy;
    //cout << MyCharge << " NeighCharge: " << NeighCharge << " " << MMFF_Energy << endl;
    return (double)EVDW;
}
