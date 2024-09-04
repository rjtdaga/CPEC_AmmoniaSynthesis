#include "Base/Constants.h"
#include "Base/Estring.h"
#include "Component.h"
#ifndef elementaryRxn_h
#define elementaryRxn_h

enum Process{NoProc, Ad1, Ad2, Des1, Des2, Diff, Diss, React, Disp, EleyAd, EleyDes, EleyAdEleyDes, DissToGas, ReactFromGas, ElecTransfer, NonSurface_e, Des2ToGas, Ad2ToGas, CompDes, CompAd};

#define Exact 1
#define Approximate 0

class ElementaryRxn
{

public:

    ElementaryRxn();
    void Toggle();
    double ActivationEnergy;
    double UserSetActivationEnergy;

    bool isReverse;
    double AvActivationEnergy;
    int RxnCounter;
    double deltaH;
    double SurfaceCharge;
    
    bool Importance; // DWDWDW	

    void SetPreexponential(double  P) 
    {
        Preexp = P;
    }

    void SetUpdPreexponential(double  P) 
    {
        Upd_Preexp = P;
    }
    
    double GetPreexponential() 
    {
        return Preexp;
    }
    
    double GetUpdPreexponential()
    {
      return Upd_Preexp;
    }
    
    void SetActivationCorrection(double E)
    {
        EactCorrection = E;
//        if(ReverseReaction != NULL)
//        {
//            ReverseReaction->EactCorrection = E; 
//        }
    }

    
    double GetActivationCorrection()
    {
        return EactCorrection;
    }

    void SetName(CharString A) {Name = A; return;}
    CharString& GetName() {return Name;}
	
    Process GetProcess() {return ReactionType;}
    Process GetReverseProcess();
    CharString GetProcessName();

    void SetProcess(Process Type)
    {
        ReactionType = Type;
    }

    
    void SetReverseReaction(ElementaryRxn* A);

    ElementaryRxn* GetReverseReaction() {return ReverseReaction;}

    void  GetParticipators(Component* Reactant[], Component* Product[]);
    void  SetParticipators(Component* Reactant[], Component* Product[]);
    
    void Print();
    
    double BondDissociationEnergy;
    int _id;
    int _idhook;

    double MacroscopicBarrier;	
    double AverageSimulationRate;
    double Rate;
    double Total_Rate;
    bool UserEquil;
    bool isExec; //To check whether the reaction is executed enough times or not
    bool isSlow; // To check whether the reaction is slow or fast
    bool isTotExec;
    bool isQE; //To check whether the reaction is quasi-equilibriated or not
    bool Checked_dummy;
    int Forw_Rev_Diff;
    int First;
    int Last;
    
    bool Active; // DWDWDW
    CharString Name; // Reaction Name
    ElementaryRxn* ReverseReaction; // Pointer to Reverse Reaction
    Component* Reactant[2]; // Reactant species
    Component* Product[2]; // Product species
    Process    ReactionType; // Type of reaction
    int NumberOfScenarios;

    int Possible;
    double StatFactor;
    Eshort RxnMatrix[90][4];
    Eshort RxnMatrixProd0[90][4];
    Eshort RxnMatrixProd1[90][4];
    double Scale;
    double CycleAveragedRate;
    double CycleRateLog[1000];
private:
    double Preexp; // Preexponential for this reaction
    double EactCorrection; // Correction to the activation energy for this reaction
    double Upd_Preexp;
    double MacroscopicPrefactor;
     
};
	  
extern ElementaryRxn* thisReaction;
	    
#endif






