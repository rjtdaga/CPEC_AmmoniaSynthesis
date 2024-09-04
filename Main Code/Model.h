#ifndef Model_H
#define Model_H


#include "Base/Facilitator.h"

class Model;
class Attributes;
class Component;
class Interactions;
class Grid;
class GridSite;
class Species;

#include "Base/Matrix.h"
#include "Base/DynArrays.h"

enum ModelType{BOC, NoModel};
enum DistanceType{HardSphere, NoCoreRepulsion, AtomAtomMMFF};

extern Matrix<Eshort> BlockingPattern[89][3];

typedef struct
{
	Eshort id;
	double alpha;
	double N;
	double A;
	double G;
	char DA;
	double AoverNsqrt;
	double A_alpha_factor; 
} MMFF_Parameters;

class Model : public Facilitator
{
public:
    
    void WritePattern();
    // void FindOptimalParameters();
    
    // DWDWDW // DWDWDW
    double CalculateBindingEnergy(GridSite* Site, bool FixedOrientation=false);
    double TempCoord[6][25][10][20][100][3];
    void GetOrientationPosition(GridSite *Central);
    
    double CalculateInteractionEnergy(GridSite* Site);
    
    Model(fstream &fin, fstream &fout);
    
    void ChangeModel(ModelType A) {theModel = A;}
    ModelType GetModel() {return theModel;}
    
    void Read(fstream &fin);
    
    
    double * DistDepend;
    void InitializeThroughSpaceMatrix();
    
        
    double GetMMFFEnergy(Eshort &MyMMFF, Eshort &NeighMMFF, double&,
        double&, double &Distance);
    
    double CalculateThroughSpace(GridSite* PassedSite);
    double CalculateThroughSpaceFixedOrientation(GridSite* PassedSite);
    double CalculateThroughSpaceInteraction(GridSite* Site, GridSite* NeighborSite);
    
    void CheckInteractions();
    
    
    
    // Sadly, these must be public for a while.
    int NumberOfOrientations[6];
    DistanceType RadialMethod;
    
private:
    ModelType theModel;
    SeqList<MMFF_Parameters> MMFF;
    double *Dielectric;
};

extern Model* theModel;
#endif
