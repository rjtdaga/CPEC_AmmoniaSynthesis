#include "ElementaryRxn.h"
#include "Component.h"

ElementaryRxn *thisReaction; // freq used, easy hook

ElementaryRxn::ElementaryRxn()
{
	for(int i =0; i < 2; ++i)
    {
        Reactant[i] = theNULLSpecies;
        Product[i]  = theNULLSpecies;
	}
	ReverseReaction = NULL;
	EactCorrection = 0.;
	Preexp         = 0.;
	ReactionType= NoProc;
	Active = true; // DWDWDW
	ActivationEnergy = 0.;
	UserSetActivationEnergy = -1.;
	Importance = false; // DWDWDW
	MacroscopicBarrier = -1.;
	MacroscopicPrefactor = -1.;
	AverageSimulationRate = 0.;
	BondDissociationEnergy = 0.;
	StatFactor = 1.;
  RxnCounter = 0;
  isSlow = 1;
  isExec = 0;
  isQE = 0;
  isTotExec = 0;
  CycleAveragedRate = 0;
  for (int i = 0; i < sizeof(CycleRateLog)/sizeof(double); i++)
  {
    CycleRateLog[i] = 0;
  }
  Scale = 1;
  Checked_dummy = 0;
  Forw_Rev_Diff = 0;
  First = 0;
  Last = 0;
  deltaH = 0.;
  AvActivationEnergy = 0.;
}

void ElementaryRxn::Toggle()
{
	if(Active) 
        Active = false; // DWDWDW
	else
        Active = true; // DWDWDW
	return;
}
	
void ElementaryRxn::GetParticipators(Component* R[], Component* P[])
{
	R[0] = Reactant[0];
	P[0] = Product[0];
	R[1] = Reactant[1];
	P[1] = Product[1];
	return;
}

void ElementaryRxn::SetParticipators(Component* R[], Component* P[]) {
	Reactant[0] = R[0];
	Reactant[1] = R[1];
	Product[0]  = P[0];
	Product[1]  = P[1];
	return;
}

CharString ElementaryRxn::GetProcessName()
{
	switch(ReactionType)
    {
        case Ad2ToGas:case Ad1:case Ad2: return "Adsorption";
        case Des2ToGas:case Des1:case Des2: return "Desorption";
        case ElecTransfer: return "Electron Transfer";
        case NonSurface_e: return "Non surface reaction";
        case Diff: return "Diffusion";
        case React: return "Reaction";
        case Diss: return "Dissociation";
        case Disp: return "Disproportionation";
        case EleyAd: return "EleyAdsorption";
        case EleyDes: return "EleyDesorption";
        case EleyAdEleyDes: return "EleyAdsorptionToEleyDesorption";
        case ReactFromGas: return "ReactFromGas";
        case DissToGas: return "DissToGas";
        case CompAd: return "Complete Adsorption";
        case CompDes: return "Complete Desorption";

        default: break;  // DW 8/24/04 
	}
	return "";
}

Process ElementaryRxn::GetReverseProcess()
{
	switch(ReactionType)
    {
        case Ad1: return Des1;
        case ElecTransfer: return ElecTransfer;
        case NonSurface_e: return NonSurface_e;
        case Ad2: return Des2;
        case Des1: return Ad1;
        case Des2: return Ad2;
        case Ad2ToGas: return Des2ToGas;
        case Des2ToGas: return Ad2ToGas;
        case React: return Diss;
        case Diss: return React;
        case DissToGas: return ReactFromGas;
        case ReactFromGas: return DissToGas;
        case EleyAdEleyDes: return EleyAdEleyDes;
        case Disp: return Disp;
        case EleyAd: return EleyDes;
        case EleyDes: return EleyAd;
        case CompDes: return CompAd;
        case CompAd: return CompDes;

        default: break; // DW 8/24/04 
	}
	return NoProc;
}	

void ElementaryRxn::SetReverseReaction(ElementaryRxn* A) {
	ReverseReaction = A;
	return;
}


void ElementaryRxn::Print()
{
    cout << " Name=" << Name
        << " id=" << _id 
        << " idhook=" << _idhook
        << " MacroBarrier=" << MacroscopicBarrier
        << " MacroPrefactor= "<< MacroscopicPrefactor
        << " AvgSimRate=" << AverageSimulationRate
        << " Rate=" << Rate
        << " Active?=" << Active
        << " Preexp=" << Preexp
        << " ActEnergy=" << ActivationEnergy
        << " EactCorrection=" << EactCorrection
        << " NumberOfScenarios=" << NumberOfScenarios
        << " Possible=" << Possible
        << " StatFactor=" << StatFactor;
    cout << endl << flush;
}

        
