
#include <iomanip>
#include <cassert>
#include <cmath>
using namespace std;

#include "Base/Constants.h"
#include "Grid.h"
#include "Component.h"


// Constructor
//RD: Redefinition error. Already defined in component.h file

int Component::CrossReference()
{

/*	if(_idhook < 0) {
	  cout << "\nBad Cross-reference" << endl;
	  cout << "Make sure a diffusion reaction is defined for" << GetName().Ec_str() << endl;
          return 0;
	} */
    
	return _idhook;
}

Component::Component() {

	VDWRadius = MolWt = CoreRadius = Atomize = 0.0;
	// ClosedShell = true; // DWDWDW
	ClosedShell = false; // RD Initialized as false avoid free radical being counted as non radical.
	Weak = Strong = NULL;
    
    //_id = NULL; // DW: Note: _id is an Eshort, which is typedef'd to short
    _id = 0; // DW: Changed to this 8/23/04
    
    _idhook = -1;
	Atom = true; // DWDWDW
	Pressure = 0.0;
  Concentration = 0.0;
  NumAtoms = 0;
	NumberOfAtoms = 0;
	TotalDesorbed = 0.0;
  Charge = 0;
  
    int i;
  for(i = 0; i < 7; ++i)
  {
      for (int j = 0; j < 10; ++j)
      {
        for (int k = 0; k < 2; ++k)
        {
          BindingEnergy[i][k][j] = 0.;
        }
      }
	}
  for (int i = 0; i < 10; ++i)
  {
    SimEnergy[i] = 0.;
    StableSiteEnergy[i] = 0.;
  }
	SimBindingEnergyCount = 0.0;
	SimBindingEnergy = 0;
	SimVacancyEnergy = 0;

    int j;
    for(j = 0; j < 7; ++j)
    {
        SimVacancyEnergyCount[j] = 0.0;
        VacancyCount[j] = 0.0;
	}
	OnSurface = false; // DWDWDW
  EquilCov = NULL;
	NumberCount = 0;
	DiffusionReaction = NULL;
	Consider = false; // DWDWDW
  for (int j = 0; j < 10; j++)
  {
      for(int k = 0; k < MAX_NUM_SPECIES; ++k)
      {
            for(i = 0; i < 7; ++i)
            {
                FavorableDist[k][i][j] = 1.0e10;
            }
    	}  
  }   
}

double Component::CalculateAdsorptionRate(double Temp)
{
    /*
    double Rate = Pressure * 6022. * 133.28947 * theGrid->SiteArea / 
    pow(2. * 3141.5 * MolWt * 8.314 * Temp , 0.5); WRONG 

    Rate = Pressure{bar} * (100000 {Pa} / 1 {bar} ) * (1 {m^2} / 10^20 {Angstrom^2}) * SiteArea {Angstrom^2}
    / square root(2 * pi * MW {g/mol} * (1 {mol} / 6.022*10^23 {atoms}) *  
    * (1 {kg} / 1000 {g}) * 8.314 {J/mol/K} * (1 {mol} / 6.022*10^23 {atoms})
    * Temperature {K})
    */
    // Finding the total coverage
    if (theSimulation->Ads_Arrhenius == 1)
    {
      if (theSimulation->Gas == true)
      {
        return Pressure;
      }
      if (theSimulation->Solvent == true)
      {
        return Concentration;
      }
    }
    else if (theSimulation->Ads_CollisionTheory == 1)
    {
      if (theSimulation->Gas == true)
      {
        double Rate = Pressure * 1e-20 * theGrid->SiteArea / pow(1.44048382555e-49 * MolWt * Temp , 0.5);
        return Rate;
      }
      if (theSimulation->Solvent == true)
      {
        double Rate = Concentration * 1e-20 * theGrid->SiteArea / pow(1.44048382555e-49 * MolWt * Temp , 0.5);
        return Rate;
      }
    }
    else
    {
      cout << "No rate constant calculation type selected" << endl;
      exit(1);
    }
}



void Component::SetOrientation(CharString A)
{
	if(A.lower() == "parallel")
    {
        theOrientation = Parallel;
	}
	else if(A.lower() == "perpendicular")
    {
        theOrientation = Perpendicular;
	}
	else
    {
        theOrientation = NoOrientation;
	}
}


/*
Component::~Component() {
	delete _id;
	delete Weak;
	delete Strong;
}
*/
void 		Component::SetName(CharString Name) {
	SpeciesName = Name;
	return;
}	
	  
void		Component::SetMolWt(const double val) {
	MolWt = val;
	return;
} 
		  
double*		Component::SetCore(const double val) {
	CoreRadius = val;
	return &CoreRadius;
}

double*		Component::SetVDW(const double val)  {
	VDWRadius = val;
	return &VDWRadius;
}

double*		Component::SetAtomizationEnergy(const double val) {
	Atomize = val;
	return &Atomize;
}
	
void 		Component::addFragments(Component *StrongE,
				        Component *WeakE) {
/*	Strong.Add(StrongE); // add is also a simplelist command
	Weak.Add(WeakE);
	assert(Strong.Size() == Weak.Size());
*/	Weak= WeakE;
	Strong= StrongE; 
	return; 
}		
	
Component& Component::operator=(Component A) {

          SetCore (A.GetCore());
          SetVDW  (A.GetVDW());
          SetAtomizationEnergy( A.GetAtomizationEnergy());
          SetMolWt(A.GetMolWt());
          SetClosedShell(A.GetClosedShell());
	  SetName (A.GetName());
	  Strong = A.Strong;
	  Weak   = A.Weak;
	  _idhook = A._idhook;
/*
	  A.Strong.ResetToFront();
	  A.Weak.ResetToFront();
	  assert(Strong.Size() == Weak.Size());
	  while(A.Strong) {
	    Strong.Add(A.Strong());
	    Weak.Add(A.Weak());
	    ++A.Strong;
	  }  */ 
	  _id = A._id;  
	  return *this;

}	  
/*	
int Component::isAtom(CharString A) {

	// A is the NULL Species string
	if(Strong && Weak) {
	  if( (Strong->GetName() == A) || (Weak->GetName() == A) )
	    return true; // DWDWDW
	}
	return False;
} 
*/

int Component::isAtom() {

	// A is the NULL Species string
	if(Strong && Weak) {
        if( (Strong->_id == 0) || (Weak->_id == 0) )
            return true; // DWDWDW
	}
	return false; // DWDWDW
} 

void Component::SetColor(double r, double g, double b) {

	MyColor.red 	= r;
	MyColor.green 	= g;
	MyColor.blue 	= b;
} 

void Component::OutputColor(double red, double green, double blue) {

    CharString SpecName = GetName();
    CharString Name;
    int Change = true; // DWDWDW
	ifstream thecolors("Colors");
    fstream colorin("Colors",ios::in|ios::out);
	if(thecolors)
    {
        double r=0, g=0, b=0;
        streampos pos;
        while(colorin >> Name && Name != "")
        {
            pos= colorin.tellp();
            colorin >> r >> g >> b;
            if(Name == SpecName) {
                colorin.seekg(pos);
                Change = false; // DWDWDW
                colorin << "\t" << setw(10) << red << "\t" 
                    << setw(10) <<  green << "\t" << setw(10) << blue << endl;
            }
        }
	}
    if(Change) {
        colorin << SpecName << "\t" << setw(10) << red << "\t" 
            << setw(10) << green << "\t" << setw(10) << blue << endl;
    }
}


Component* Component::GetBindingAtom() {

	Component* TempA = this;
        while(!TempA->isAtom()) {
	  TempA = TempA->Strong;
	}
	return TempA;
}








