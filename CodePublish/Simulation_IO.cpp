void Simulation::WriteCoverage()
{
    ofstream covout(Coverage_File, ios::app);
    
    if(Coverage_Init)
    {
        covout << "# ***Surface Coverages as a Function of Time***\n"
            << endl;
        covout << setiosflags(ios::scientific);
        covout << "#" << setw(11) << "time" << " ";
        covout << setw(11) << "temp (K)" << " ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            thisComponent = (Listing[i]);
            covout << setw(7) << thisComponent->GetName().Ec_str() << " ";
            
        }
        covout << endl;
        covout << endl << endl;
        covout << setiosflags(ios::scientific);
        Coverage_Init = false; // DWDWDW
    }
    else
    {
        covout << setw(11) << setprecision(7) << Time << " ";
        covout << setw(11) << setprecision(7) << GetTemperature() << " ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            thisComponent = (Listing[i]);
            covout << setw(12) << setprecision(8) << theGrid->GetCoverage(thisComponent) << " ";
        }
        covout << endl;
    }
    covout.close();
}

void Simulation::WriteEmptySites()
{   
    std::ofstream fe(EmptySites_File, ios::app);
    fe << Time << setw(10);
    int Layer = theGrid->AdsorbateLayer;
    int EmptyCount[3];
    for (int i = 0; i < theGrid->EDim[0]; i++)
    {
      for (int j = 0; j < theGrid->EDim[1]; j++)
      {
        GridSite *thisSite = theGrid->Surface[Layer][i][j];
        if (thisSite != NULL && thisSite->GetType() == theNULLSpecies)
        {
          if (theGrid->FindCoordination(i,j) == Atop)
            ++EmptyCount[0];
          else if (theGrid->FindCoordination(i,j) == Bridge)
            ++EmptyCount[1];
          else if (theGrid->FindCoordination(i,j) == Hollow)
            ++EmptyCount[2];
        }
      }
    }
    fe << EmptyCount[0] << setw(10) << EmptyCount[1] << setw(10) << EmptyCount[2] << endl;
    fe.close();
    return;
}

void Simulation::WriteDesorbed()
{

    ofstream desout(Desorbed_File, ios::app);

    if(Desorbed_Init)
    {

        desout << "# *** Desorbed as a Function of Time***\n" << endl;
        desout << setiosflags(ios::scientific);
        desout << "#" << setw(15) << "time" << " ";
        desout << setw(15) << "temp (K)" << " " << setw(15) << "Charge ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            thisComponent = Listing[i];
            desout << setw(7) << thisComponent->GetName().Ec_str() << " ";
        }
        desout << endl;

        Desorbed_Init = false; // DWDWDW
    }
    else
    {
        desout << setw(16) << setprecision(12) << Time << " ";
        desout << setw(16) << setprecision(12) << GetTemperature() << " " << setw(16) << SurfaceCharge << " ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            thisComponent = (Listing[i]);
            desout << setw(12) << setprecision(8) << thisComponent->
                TotalDesorbed << " ";
        }
        desout << endl;
    }
    desout.close();
}

void Simulation::WriteReactions()
{
    // DW 8/24/04 Had to tweak these ofstream constructions.
    ofstream cfout(Reaction_File, /*ios::scientific | */ ios::app);
    cfout << setiosflags(ios::scientific);
    const char *d1 = "EventsByRate";
    char EventByrate[100];
    strcpy(EventByrate, d1);
    strcat(EventByrate, rnk);
    ofstream pout(EventByrate, /*ios::scientific | */ ios::app);
    pout << setiosflags(ios::scientific);
    
    if(Reaction_Init)
    {
        cfout << "\n\n#*** Reaction Rates ***\n" << endl;
        pout << "\n\n#*** Reaction Rates ***\n" << endl;
        cfout << "#";
        pout << "#";
        cfout << setw(25) << "Temperature (K)"
            << setw(25) << setprecision(5) << theSimulation->Temperature;
        pout << setw(25) << "Temperature (K)"
            << setw(25) << setprecision(5) << theSimulation->Temperature;
        cfout << "\n\n";
        pout << "\n\n";
        cfout << "# Reaction Name";
        cfout << setw(25) << "# Forward Rate (1/s)" << endl;
        cfout << endl << endl;
        pout << "# Reaction Name";
        pout << setw(25) << "# Forward Rate (1/s)" << endl;
        pout << endl << endl;
        thisReaction = NULL;
        for(int i = 0; i < NumberOfSimReactions; ++i)
        {
            thisReaction = SimReactions[i];
            cfout << "#  " << i + 1 << " ";
            pout << "#  " << i + 1 << " ";
            CharString RxnName = thisReaction->GetName();

            cfout << setw(40) << RxnName.Ec_str() << endl;
            pout << setw(40) << RxnName.Ec_str() << endl;
            ++(*Rxns);
        }
        cfout << endl;
        pout << endl;
        Reaction_Init = false; // DWDWDW
    }
    else
    {
        cfout << setw(13) << setprecision(9) << Time << " ";
        cfout << setw(13) << setprecision(9) << GetTemperature() << " ";
        pout << setw(13) << setprecision(9) << Time << " ";
        pout << setw(13) << setprecision(9) << GetTemperature() << " ";
        Rxns->ResetToFront();
        double TotalRate = 0.;

        for(int i = 0; i < NumberOfSimReactions; ++i)
        {
            thisReaction = SimReactions[i];
            TotalRate += thisReaction->AverageSimulationRate;
            cfout << setw(7) << setprecision(3) << thisReaction->
                AverageSimulationRate << " ";
        }
        for(int i2 = 0; i2 < NumberOfSimReactions; ++i2)
        {
            thisReaction = SimReactions[i2];
            double tempvar = EventsByRate[i2];

            if(TotalRate > 0.0)
                EventsByRate[i2] +=
                    VTS_Chance * thisReaction->AverageSimulationRate /
                    TotalRate;
            if(isnan(EventsByRate[i2]))
                EventsByRate[i2] = tempvar;
            pout << setw(7) << setprecision(3) << EventsByRate[i2] << " ";
        }
        cfout << endl;
        pout << endl;
    }
    cfout.close();
    pout.close();
}

void Simulation::WriteBarrier()
{

    // DW 8/24/04 Tweaked construction
    ofstream cfout(Barrier_File, /*ios::scientific | */ ios::app);
    cfout << setiosflags(ios::scientific);
    
    if(Barrier_Init)
    {
        cfout << "\n\n#*** Barriers ***\n" << endl;
        cfout << "#";
        cfout << setw(25) << "Temperature (K)"
            << setw(25) << setprecision(5) << theSimulation->Temperature;
        cfout << "\n\n";
        cfout << "# Reaction Name";
        cfout << setw(25) << "# Forward Barrier (kcal/mol)" << endl;
        cfout << endl << endl;
        thisReaction = NULL;
        for(int i = 0; i < NumberOfSimReactions; ++i)
        {
            thisReaction = SimReactions[i];
            cfout << "#  " << i + 1 << " ";
            CharString RxnName = thisReaction->GetName();

            cfout << setw(40) << RxnName.Ec_str() << endl;
            ++(*Rxns);
        }
        cfout << endl;
        Barrier_Init = false; // DWDWDW
    }
    else
    {
        cfout << setw(13) << setprecision(9) << Time << " ";
        cfout << setw(13) << setprecision(9) << GetTemperature() << " ";
        for(int i = 0; i < NumberOfSimReactions; ++i)
        {
            thisReaction = SimReactions[i];
//            if(thisReaction->MacroscopicBarrier < 1.e5 &&
//                thisReaction->MacroscopicBarrier > 0.)
//            {
//                cfout << setw(7) << setprecision(3) << thisReaction->
//                    MacroscopicBarrier << " ";
//            }
//            else
//            {
//                cfout << setw(7) << setprecision(3) << "x" << " ";
//            }
	    cfout << setw(7) << thisReaction->AvActivationEnergy << " ";
         cfout << setw(7) << thisReaction->GetPreexponential() << " ";
        }
        cfout << endl;
    }
    cfout.close();
}

void Simulation::WriteBindingEnergies()
{

    thisComponent = NULL;
    ofstream fe(BindingEnergy_File, ios::app);

    fe << setiosflags(ios::scientific);
    if(BindingEnergy_Init)
    {
        fe << "# *** Binding Energy as a function of Simulation Time ***\n";
        fe << "#" << setw(13) << "time" << " ";
        fe << setw(13) << "temp (K)" << " ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            thisComponent = (Listing[i]);
            fe << setw(7) << thisComponent->GetName().Ec_str() << " ";
        }

        BindingEnergy_Init = false; // DWDWDW
    }
    else
    {
        fe << setw(13) << setprecision(9) << Time << " ";
        fe << setw(13) << setprecision(9) << GetTemperature() << " ";
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            fe << setw(7) << setprecision(3) << Listing[i]->
                SimBindingEnergy << " ";
        }
    }
    fe << endl;
    fe.close();
    return;
}
