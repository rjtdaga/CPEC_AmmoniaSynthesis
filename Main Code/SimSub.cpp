void Simulation::Read(fstream &fin)
{
    fin.seekg(inp);
    CharString Next;
    CharString End = "END";
    double value;
    
    // DW 8/24/04: SpecVal is not used, so commented it out
    //int SpecVal = -1;

    while(fin >> Next && Next.lower() != End.lower())
    {
        if(Next.lower() == "pressure")
        {
            fin >> Next >> value;
            theSpecies->get(Next)->SetPressure(value);
        }
        if(Next.lower() == "concentration")
        {
            fin >> Next >> value;
            theSpecies->get(Next)->SetConcentration(value);
        }
        if(Next.lower() == "pressure_study")
        {
            RemoveGas = false; // DWDWDW
        }
        if(Next.lower() == "concentration_study")
        {
            ConstConc = false; // DWDWDW
        }
        if(Next.lower() == "solvent")
        {
            Solvent = true;
        }
        if(Next.lower() == "gas")
        {
            Gas = true;
        }
        if(Next.lower() == "ads_arrhenius")
        {
            Ads_Arrhenius = true;
            Ads_CollisionTheory = false;
        }
        if(Next.lower() == "ads_collisiontheory")
        {
            Ads_CollisionTheory = true;
            Ads_Arrhenius = false;
        }
        if(Next.lower() == "temperature")
        {
            fin >> value;
            Temperature = value;
        }
        if(Next.lower() == "coverage")
        {
            
            fin >> Next >> value;
            theSpecies->get(Next)->NumberOfAtoms = 
                (int)(value / theGrid->DeltaCoverage);
        }
        if(Next.lower() == "equilibrium_coverage")
        {
            
            fin >> Next >> value;
            theSpecies->get(Next)->EquilCov =
                (int)(value / theGrid->DeltaCoverage);
        }
        
        if(Next.lower() == "equilibration_factor")
        {
            fin >> value;
            equilibration_factor = value;
        }
        if(Next.lower() == "time")
        {
            fin >> value;
            StopTime = value;
        }
        if(Next.lower() == "delay")
        {
            fin >> value;
            DelayTime = (double)value;
        }
        if(Next.lower() == "beta")
        {
            fin >> value;
            Beta = value;
        }
        if(Next.lower() == "time_step")
        {
            fin >> value;
            running_time_step = value;
            input_time_step = value;
        }
        if(Next.lower() == "debug")
        {
            Debug = true; // DWDWDW
        }
        if(Next.lower() == "populate")
        {
            Population = true; // DWDWDW
        }
        if(Next.lower() == "eliminate")
        {
            /*
            * fin >> Next;
            * Component* SpecToElim = theSpecies->Get(Next);
            * SpecVal = SpecToElim->_idhook;
            */
        }
//        if(Next.lower() == "restart")
//        {
//            theGrid->Restart();
//        }
        if(Next.lower() == "userset")
        {
            UserSet = true;
        }
        if(Next.lower() == "nointeraction")
        {
            NoInteraction = true;
        }
        if(Next.lower() == "numbercharge")
        {
            fin >> NumberChargeType;
        }
        if(Next.lower() == "chargepattern")
        {
            for (int i = 0; i < NumberChargeType; i++)
            {
              fin >> value;
              ChargePattern[i] = value;
            }
        }
        if(Next.lower() == "timeperiod")
        {
            for (int i = 0; i < NumberChargeType; i++)
            {
              fin >> value;
              TimePeriod[i] = value;
            }
        }
        if(Next.lower() == "unknownbe_boc")
        {
            UnknownBE_BOC = 1;;
        }
        //---------------------------------------------
        //Reading Dybeck Parameters
        if(Next.lower() == "rescalingcyclesize")
        {
            fin >> CycleTime;
        }
        if(Next.lower() == "minexecsuperbasin")
        {
            fin >> MinExecSuperbasin;
        }
        if(Next.lower() == "maxrxnqueuesize")
        {
            fin >> MaxRxnQueueSize;
        }
        if(Next.lower() == "QEdiff")
        {
            fin >> QEMaxDiff;
        }
        if(Next.lower() == "k_buff")
        {
            fin >> k_buff;
        }
        //---------------------------------------------
    }
    
}


void Simulation::Read_Analysis(fstream &fin)
{
  CharString Next;
  CharString Param;
  CharString End = "End";
  while(fin >> Next && Next.lower() != End.lower())
  {
    if(Next.lower() == "analyze")
    {
        fin >> Param;
    }
    if (Next.lower() == "period")
    {
      if (Param == "rescalingcyclesize")
      {
        int Period;
        fin >> Period;
        CycleTime = CycleTime + rank*Period;
      }
      else if (Param.lower() == "minexecsuperbasin")
      {
        int Period;
        fin >> Period;
        MinExecSuperbasin = MinExecSuperbasin + rank*Period;
      }
      else if (Param.lower() == "maxrxnqueuesize")
      {
        int Period;
        fin >> Period;
        MaxRxnQueueSize = MaxRxnQueueSize + rank*Period;
      }
      else if (Param.lower() == "qediff")
      {
        int Period;
        fin >> Period;
        QEMaxDiff = QEMaxDiff + rank*Period;
      }
      else if (Param.lower() == "k_buff")
      {
        double Period;
        fin >> Period;
        k_buff = k_buff + rank*Period;
      }
      else if (Param.lower() == "reactionactenergy")
      {
        int RxnNum;
        fin >> RxnNum;
        int forw;
        fin >> forw;
        ElementaryRxn *thisRxn;
        if (forw == 0)
        {
          thisRxn = SimReactions[(RxnNum-1)*2];
        }
        else
        {
          thisRxn = SimReactions[(RxnNum-1)*2+1];
        }
        double Period;
        fin >> Period;
        thisRxn->UserSetActivationEnergy = thisRxn->UserSetActivationEnergy + rank*Period;
      }
      else if (Param.lower() == "preexpconstant")
      {
        int RxnNum;
        fin >> RxnNum;
        int forw;
        fin >> forw;
        ElementaryRxn *thisRxn;
        thisRxn = SimReactions[(RxnNum-1)*2+forw];
        double Period;
        fin >> Period;
        thisRxn->SetPreexponential(thisRxn->GetPreexponential() + rank*Period);
        thisRxn->SetUpdPreexponential(thisRxn->GetUpdPreexponential() + rank*Period);
      }
      else if (Param.lower() == "pressure")
      {
        CharString SpecName;
        fin >> SpecName;
        Component* thisComp;
        for (int i = 0; i < NumberOfSpecies; ++i)
        {
          Component* dumm = Listing[i];
          if (dumm->GetName() == SpecName)
            thisComp = dumm;
        }
        int Period;
        fin >> Period;
        thisComp->SetPressure(thisComp->GetPressure() + rank*Period);
      }
      else if (Param.lower() == "coverage")
      {
        CharString SpecName;
        fin >> SpecName;
        Component* thisComp;
        for (int i = 0; i < NumberOfSpecies; ++i)
        {
          Component* dumm = Listing[i];
          if (dumm->GetName() == SpecName)
          {
            thisComp = dumm;
            break;
          }
        }
        int Period;
        fin >> Period;
        double value = thisComp->NumberOfAtoms*theGrid->DeltaCoverage;
        thisComp->NumberOfAtoms = (int)((value + rank*Period) / theGrid->DeltaCoverage);
      }
      else if (Param.lower() == "temperature")
      {
        int Period;
        fin >> Period;
        Temperature = Temperature + rank*Period;
      }
      else if (Param.lower() == "scale")
      {
        double Period;
        fin >> Period;
        theInteractions->SetScale(theInteractions->GetScale() + rank*Period);
      }
      else if (Param.lower() == "timeperiod")
      {
        double Period;
        int ChargeIndex;
        fin >> Period;
        fin >> ChargeIndex;
        TimePeriod[ChargeIndex-1] = TimePeriod[ChargeIndex-1] + Period*rank;
        cout << ChargePattern[ChargeIndex-1] << " " << TimePeriod[ChargeIndex-1] << endl;
      }
      else
      {
        cout << "The analysis variable does not exist!!" << endl;
        exit(1);
      }
    }
  }
}

double Simulation::StableSiteEnergy(Component * PassedComponent)
{
    double charge = theSimulation->SurfaceCharge;
    double HighestBE = 0.;
    double BE = 0.;
    CharString S = theGrid->GetSurface();
    if(S != "100" && S != "111" && S.lower() != "graphene")
    {
        cout << "Undefined Surface in Sim" << endl;
        exit(1);
    }
    // Finding the index corresponding to provided charge
    int index = theSimulation->FindChargeIndex(charge);
    for(int i = 1; i < 4; ++i)
    {
        if (S.lower() == "graphene" && i == 3)
            i = 6;
        if (S == "100" && i == 3)
            i = 4;
        BE = thisComponent->BindingEnergy[i][0][index];
        if(BE > HighestBE)
        {
            HighestBE = BE;
        }
    }
    return HighestBE;
}

void Simulation::ConsistentBarrierModification(fstream &fout)
{  
    CharString SurfaceType = theGrid->GetSurface();
    double SurfCharge = SurfaceCharge;
    for (int ind = 0; ind < NumberChargeType; ind++)
    {
      SurfaceCharge = ChargePattern[ind];
      theSpecies->theSpecies.ResetToFront();
      while(theSpecies->theSpecies)
      {
          thisComponent = theSpecies->theSpecies.GetPtr();
          ++theSpecies->theSpecies;
          thisComponent->StableSiteEnergy[ind] = StableSiteEnergy(thisComponent);
          thisComponent->SimVacancyEnergy = thisComponent->StableSiteEnergy[ind];
          thisComponent->SimBindingEnergy = thisComponent->StableSiteEnergy[ind];
          thisComponent->SimEnergy[ind] = thisComponent->StableSiteEnergy[ind];
          //cout << thisComponent->GetName() << " ind: " << ind << " BE: " << thisComponent->StableSiteEnergy[ind] << endl;
      }
    }
    SurfaceCharge = SurfCharge;
    theReactions->SetBondDissociationEnergies();
//    Rxns->ResetToFront();
//    while(*Rxns)
//    {
//      ElementaryRxn * thisRxn = Rxns->GetPtr();
//      cout << thisRxn->GetName() << " DeltaH: " << thisRxn->deltaH << " Bond Diss Energy: " << thisRxn->BondDissociationEnergy << endl;
//      ++(*Rxns);
//    }
//    exit(1);
    for (int ind = 0; ind < NumberChargeType; ind++)
    {
      SurfaceCharge = ChargePattern[ind];
      theSpecies->theSpecies.ResetToFront();
      thisComponent = theNULLSpecies;
      thisSite = NULL;
      if (NoInteraction == false) 
      {
        Component *Reactants[2];
        
        Reactants[0] = theNULLSpecies;
        Reactants[1] = theNULLSpecies;
        Component *Products[2];
        
        Products[0] = theNULLSpecies;
        Products[1] = theNULLSpecies;
        thisReaction = NULL;
        Rxns->ResetToFront();
        while(*Rxns)
        {
            thisReaction = Rxns->GetPtr();
            ++(*Rxns);
            if(thisReaction->SurfaceCharge == SurfaceCharge)
            {
             double Eact = thisReaction->UserSetActivationEnergy;
              
              if(Eact > -0.1)
              {                       // if User - Specified Activation Energy
                  thisReaction->ActivationEnergy = Eact;
                  thisReaction->GetParticipators(Reactants, Products);
                  double BE_R[2], BE_P[2];
                  
                  BE_R[0] = Reactants[0]->StableSiteEnergy[ind];
                  BE_R[1] = Reactants[1]->StableSiteEnergy[ind];
                  BE_P[0] = Products[0]->StableSiteEnergy[ind];
                  BE_P[1] = Products[1]->StableSiteEnergy[ind];
                  double Resist, D;
                  double FwRxn, Correct, delH;
                  D = thisReaction->BondDissociationEnergy;
                  switch (thisReaction->GetProcess())
                  {
                      // Set Thermodynamics
                      // BE is Binding Energy, D is Dissociation Energy,
                      // Resist is resistance term 
                      // Eabg is gas phas adsorption barrier
                      case NonSurface_e:
                      case Ad1:
                      {
                          thisReaction->SetActivationCorrection(Eact);
                          break;
                      }
                      case ElecTransfer:
                      {
                          Resist = (BE_R[0] * BE_P[0]) / (BE_R[0] + BE_P[0]); 
                          delH = BE_R[0] - BE_P[0];
                          Correct = Eact - 0.5 * (delH + Resist);
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Des1:
                      {
                          Correct = Eact - BE_R[0];
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case React:
                      {
                          Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                          delH = BE_R[0] + BE_R[1] - BE_P[0] + D;
                          FwRxn = 0.5 * (delH + Resist);
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case ReactFromGas:
                      {
                          delH = BE_R[1] - BE_P[0] + D;
                          FwRxn = delH;
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Diss:
                      {
                          Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                          delH = BE_R[0] + D - BE_P[0] - BE_P[1];
                          FwRxn = 0.5 * (delH + Resist);
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case DissToGas:
                      {
                          delH = BE_R[0] + D - BE_P[1];
                          FwRxn = delH;
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Disp:
                      {
                          delH = BE_R[0] + BE_R[1] + D - BE_P[0] - BE_P[1];
                          if(D >= 0.)
                          {
                              Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                          }
                          else
                          {
                              Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                          }
                          FwRxn = 0.5 * (delH + Resist);
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case EleyAd:
                      case EleyDes:
                      case EleyAdEleyDes:
                      {
                          delH = BE_R[0] + BE_R[1] + D - BE_P[0] - BE_P[1];
                          if(D >= 0.)
                          {
                              Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                          }
                          else
                          {
                              Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                          }
                          if(thisReaction->GetProcess() == EleyAd)
                          {
                              delH = BE_R[1] + D - BE_P[0] - BE_P[1];
                              FwRxn = 0.5 * (delH + Resist);
                          }
                          else if(thisReaction->GetProcess() == EleyDes)
                          {
                              delH = BE_R[1] + BE_R[0] + D - BE_P[1];
                              FwRxn = 0.5 * (delH + Resist);
                          }
                          else if(thisReaction->GetProcess() == EleyAdEleyDes)
                          {
                              delH = BE_R[1] + D - BE_P[1];
                              FwRxn = 0.5 * (delH + Resist);
                          }
                          Correct = Eact - FwRxn;
      
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Ad2:
                      {
                          Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                          delH = D - BE_P[0] - BE_P[1];
//                          if (thisReaction->GetName() == "N2 = N * + N *")
//                            cout << thisReaction->GetName() << " SimSub " << D << " " << BE_P[0] << " " << BE_P[1] << endl;
                          FwRxn = 0.5 * (delH + Resist);
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Des2:
                      {
                          Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                          delH = BE_R[0] + BE_R[1] + D;
                          FwRxn = 0.5 * (delH + Resist);
                          Correct = Eact - FwRxn;
                          thisReaction->SetActivationCorrection(Correct);
                          break;
                      }
                      case Diff:
                      {
                          // a multiplicative speedup
                          if(SurfaceType == "111")
                          {
                              Correct = Eact - BE_R[0] * (1.6666667 - 1.5);
                          }
                          else
                          {
                              Correct = Eact - BE_R[0] * (1.75 - 1.5);
                          }
                          thisReaction->SetActivationCorrection(Correct);
                      }
      
                      default: break;
                  }
              }
              //cout << ind << " " << thisReaction->GetName() << " Charge: " << SurfaceCharge << " Correction: " << thisReaction->GetActivationCorrection() << endl;
            }
        }
        
        if (UserSet == false)
        {
          Rxns->ResetToFront();
          while(*Rxns)
          {
              thisReaction = Rxns->GetPtr();
              ++(*Rxns);
              if(thisReaction->SurfaceCharge == SurfaceCharge)
              {
                theSpecies->theSpecies.ResetToFront();
                double Eact = thisReaction->GetActivationCorrection();
                thisReaction->GetParticipators(Reactants, Products);
                double BE_R[2], BE_P[2];
                BE_R[0] = Reactants[0]->StableSiteEnergy[ind];
                BE_R[1] = Reactants[1]->StableSiteEnergy[ind];
                BE_P[0] = Products[0]->StableSiteEnergy[ind];
                BE_P[1] = Products[1]->StableSiteEnergy[ind];
                double delH, RvRxn, FwRxn, Correct;
                double Resist, D;
                D = thisReaction->BondDissociationEnergy;
                Correct = Eact;
                switch (thisReaction->GetProcess())
                {
                    // Set Thermodynamics
                    // BE is Binding Energy, D is Dissociation Energy, Resist
                    // is resistance term 
                    case NonSurface_e:
                    case Ad1:
                    {
                        if(Correct < -0.1)
                            Correct = 0.;
                        thisReaction->ActivationEnergy = Correct;
                        break;
                    }
                    case ElecTransfer:
                    {
                        Resist = (BE_R[0] * BE_P[0]) / (BE_R[0] + BE_P[0]); 
                        delH = BE_R[0] - BE_P[0];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case Des1:
                    {
                        if(Correct < -0.1)
                            Correct = 0.;
                        thisReaction->ActivationEnergy = BE_R[0] + Correct;
                        break;
                    }
                    case React:
                    {
                        Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        delH = BE_R[0] + BE_R[1] - BE_P[0] + D;
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case ReactFromGas:
                    {
                        delH = BE_R[1] - BE_P[0] + D;
                        FwRxn = delH + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case DissToGas:
                    {
                        delH = BE_R[0] + D - BE_P[1];
                        FwRxn = delH + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case Diss:
                    {
                        Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        delH = BE_R[0] + D - BE_P[0] - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case Disp:
                    {
                        if(D >= 0.)
                        {               // Resist must be defined for positive D
                            Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        }
                        else
                        {
                            Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        }
        
                        delH = BE_R[0] + BE_R[1] + D - BE_P[0] - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                        {               // Adj. BOC negative Barriers 
                            RvRxn -= FwRxn;
                            FwRxn = 0.;
                        }
                        if(RvRxn < 0.)
                        {
                            FwRxn -= RvRxn;
                            RvRxn = 0.;
                        }
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case EleyAdEleyDes:
                    {
                        if(D >= 0.)
                        {               // Resist must be defined for positive D
                            Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        }
                        else
                        {
                            Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        }
        
                        delH = BE_R[1] + D - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                        {               // Adj. BOC negative Barriers 
                            RvRxn -= FwRxn;
                            FwRxn = 0.;
                        }
                        if(RvRxn < 0.)
                        {
                            FwRxn -= RvRxn;
                            RvRxn = 0.;
                        }
                        thisReaction->ActivationEnergy = FwRxn;
        
                        break;
                    }
                    case EleyAd:
                    {
                        if(D >= 0.)
                        {               // Resist must be defined for positive D
                            Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        }
                        else
                        {
                            Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        }
        
                        delH = BE_R[1] + D - BE_P[0] - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                        {               // Adj. BOC negative Barriers 
                            RvRxn -= FwRxn;
                            FwRxn = 0.;
                        }
                        if(RvRxn < 0.)
                        {
                            FwRxn -= RvRxn;
                            RvRxn = 0.;
                        }
                        thisReaction->ActivationEnergy = FwRxn;
        
                        break;
                    }
                    case EleyDes:
                    {
                        if(D >= 0.)
                        {               // Resist must be defined for positive D
                            Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        }
                        else
                        {
                            Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        }
                        delH = BE_R[0] + BE_R[1] + D - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                        {               // Adj. BOC negative Barriers 
                            RvRxn -= FwRxn;
                            FwRxn = 0.;
                        }
                        if(RvRxn < 0.)
                        {
                            FwRxn -= RvRxn;
                            RvRxn = 0.;
                        }
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case Ad2:
                    {
                        Resist = (BE_P[0] * BE_P[1]) / (BE_P[0] + BE_P[1]);
                        delH = D - BE_P[0] - BE_P[1];
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers           
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
                    case Des2:
                    {
                        Resist = (BE_R[0] * BE_R[1]) / (BE_R[0] + BE_R[1]);
                        delH = BE_R[0] + BE_R[1] + D;
                        FwRxn = 0.5 * (delH + Resist) + Correct;
                        RvRxn = FwRxn - delH;
                        if(FwRxn < 0.)
                            FwRxn = 0.; // Adj. BOC negative Barriers   
                        if(RvRxn < 0.)
                            FwRxn -= RvRxn;
                        thisReaction->ActivationEnergy = FwRxn;
                        break;
                    }
        
                    default: break;  // DW 8/24/04
                }
                //cout << thisReaction->GetName() << " Charge: " << thisReaction->SurfaceCharge << " delH: " << delH << " Barrier: " << thisReaction->ActivationEnergy << endl; 
              }
          }
        }
      }
  }
  SurfaceCharge = SurfCharge;
  //exit(1);
}


int Simulation::SuperFastDiffuseEquilibration(double bias)
{
    // Unimportant, resets ghost atoms
    double SC = SurfaceCharge;
    int TotalNumberOfAtoms = 0;
    int ind = FindChargeIndex(SC);
    
    for(int i = 0; i < NumberOfSpecies; ++i)
        TotalNumberOfAtoms += Listing[i]->NumberOfAtoms;
    // Define the Adsorbate List
    int Layer = theGrid->AdsorbateLayer;
    GridSite *PickSite = NULL;
    int index = 0;

    for(int y = 0; y < theGrid->EDim[1]; ++y)
    {
        for(int x = 0; x < theGrid->EDim[0]; ++x)
        {
            PickSite = theGrid->Surface[Layer][x][y];
            if(PickSite != NULL && PickSite->GetType() != theNULLSpecies)
            {
                Adsorbates[index] = PickSite;
                ++index;
            }
        }
    }
    const int ArrayMax = 2 * PatternDim;
    double InitialEnvironmentEnergy[ArrayMax];
    GridSite *Spectators[ArrayMax];
    int NumberOfSpectators;

    // Return if No adsorbates in list
    if(index == 0)
    {
        return 0;
    }
    // Precalculate the binding energies of these adsorbates
    for(int i2 = 0; i2 < index; ++i2)
    {
        PickSite = Adsorbates[i2];
        theModel->CalculateBindingEnergy(PickSite);
    }

    IRandom PickNumber;

    PickNumber.Randomize();
    IRandom PickNumberSurr;

    PickNumberSurr.Randomize();
    PickNumber.SetInterval(0, index - 1);
    int BailOut = (int)(bias * (double)TotalNumberOfAtoms * 89. * (double)equilibration_factor);
    int Events = 0;
    int ChosenAdsorbateNumber = 0;
    GridSite *SurroundingSite = NULL;
    GridSite *IntermediateSite = NULL;
    GridSite *IntermediateSiteNeighbor = NULL;
    Component *A;
    double InitialSpectatorBindingEnergy = 0.;
    double FinalSpectatorBindingEnergy = 0.;
    double theRT = 1. / (R_GasConst * Temperature);
    double Prob = 0;
    double dist;
    double TestEnergy;
    int NumSurrPick;
    bool Diffuse = 1;
    double BE = 0.;
    // DW 8/24/04: cutoff is only used in commented out code, so I
    // commented out cutoff.
    // double cutoff = 1.2 * theGrid->MMdistance;
    
    // int itest = 0;
    while(BailOut--)
    {
        ChosenAdsorbateNumber = PickNumber.Draw();
        PickSite = Adsorbates[ChosenAdsorbateNumber];
        PickNumberSurr.SetInterval(0, PickSite->SurroundSize - 1);
        dist = 1.e6;
        while(dist > theGrid->Cutoff*theGrid->MMdistance) 
        {
          NumSurrPick = PickNumberSurr.Draw();
          SurroundingSite = PickSite->Surround[NumSurrPick];
          dist = FindDistance(PickSite, NumSurrPick);
        }
        
        A = PickSite->GetType();
        assert(A != theNULLSpecies);
        if((SurroundingSite->Type == theNULLSpecies) && (A->BindingEnergy[SurroundingSite->BondSize][0][ind] > 0.0))
        {
            BE = theModel->CalculateBindingEnergy(PickSite, true);
            PickSite->SetType(theNULLSpecies);
            SurroundingSite->SetType(A);
            
            double SurrE = theModel->CalculateBindingEnergy(SurroundingSite);
            //double delH = PickSite->BindingEnergy - SurrE;
            double delH = BE - SurrE;
            
//            InitialSpectatorBindingEnergy = 0.;
//            FinalSpectatorBindingEnergy = 0.;
//            // Identify Spectators
//            IntermediateSite = PickSite;
//            NumberOfSpectators = 0;
//            for(int i = IntermediateSite->SurroundSize; i--; /*none*/ )
//            {
//                IntermediateSiteNeighbor = IntermediateSite->Surround[i];
//                if(IntermediateSiteNeighbor->Type->_id)
//                {
//                    if(IntermediateSiteNeighbor != PickSite &&
//                        IntermediateSiteNeighbor != SurroundingSite)
//                    {
//                        // System Before Movement
//                        Spectators[NumberOfSpectators] = IntermediateSiteNeighbor;
//                        ++NumberOfSpectators;
//                    }
//                }
//            }
//            IntermediateSite = SurroundingSite;
//            for(int i3 = IntermediateSite->SurroundSize; i3--;  )
//            {
//                IntermediateSiteNeighbor = IntermediateSite->Surround[i3];
//                if(IntermediateSiteNeighbor->Type->_id)
//                {
//                    if(IntermediateSiteNeighbor != PickSite && IntermediateSiteNeighbor != SurroundingSite)
//                    {
//                        // System Before Movement
//                        bool Continue = true; // DWDWDW // DWDWDW
//                        for(int jj = 0 ; jj < NumberOfSpectators && Continue ; ++jj)
//                        {
//                            if(IntermediateSiteNeighbor == Spectators[jj])
//                            {
//                                Continue = false; // DWDWDW
//                            }
//                        }
//                        if(Continue)
//                        {
//                            Spectators[NumberOfSpectators] = IntermediateSiteNeighbor;
//                            ++NumberOfSpectators;
//                        }
//                    }
//                }
//            }
//            IntermediateSite = SurroundingSite;
//            for(int i4 = NumberOfSpectators; i4--;  )
//            {
//                IntermediateSiteNeighbor = Spectators[i4];
//                
//                InitialEnvironmentEnergy[i4] = theModel->CalculateBindingEnergy(IntermediateSiteNeighbor, true);
//                
//                InitialSpectatorBindingEnergy += InitialEnvironmentEnergy[i4];
//                
//                TestEnergy = InitialEnvironmentEnergy[i4]; // DWDWDW
//                
//                if (TestEnergy < 0)
//                {
//                  Diffuse = 0;
//                  break;
//                }
//                
//                FinalSpectatorBindingEnergy += TestEnergy;
//            }
//            
            //double Keq = exp((-delH + FinalSpectatorBindingEnergy - InitialSpectatorBindingEnergy) * theRT);
            double Keq = exp(-delH * theRT);
            Keq = Keq/(1+Keq);
            
            if(isnan(Keq))
            {
                Keq = double_limit;
            }
            if (Diffuse == 1)
            {
              if(Keq >= 1.0)
              {
                  // Move List with Adsorbate
                  Adsorbates[ChosenAdsorbateNumber] = SurroundingSite;    
                  Events++;
              }
              else
              {
                  Prob = FindRandomNumber();
                  if(Prob < Keq)
                  {
                      // Move List with Adsorbate
                      Adsorbates[ChosenAdsorbateNumber] = SurroundingSite;
                      Events++;
                  }
                  else
                  {
                      // Movement not possible
                      PickSite->SetType(A);
                      SurroundingSite->SetType(theNULLSpecies);
//                      for(int i = NumberOfSpectators; i--;)
//                      {
//                          Spectators[i]->BindingEnergy = InitialEnvironmentEnergy[i];
//                      }
                  }
              }
            }
            else
            {
              PickSite->SetType(A);
              SurroundingSite->SetType(theNULLSpecies);
//              for(int i = NumberOfSpectators; i--;)
//              {
//                  Spectators[i]->BindingEnergy = InitialEnvironmentEnergy[i];
//              }
            }
        }
        else if(SurroundingSite == PickSite)
        {
            int MyCoord = PickSite->BondSize - 1;
            int NumOr = theModel->NumberOfOrientations[MyCoord];
            
            PickNumberSurr.SetInterval(0, NumOr - 1);
            int MyOr = PickSite->Orientation;
            int NewOr = PickSite->Orientation;
            double InitialEnergy = theModel->CalculateBindingEnergy(PickSite, true); 
            
            while(NewOr == MyOr)
            {
                NewOr = PickNumberSurr.Draw();
            }
            
            PickSite->Orientation = NewOr;
            
            double FinalEnergy =
                theModel->CalculateBindingEnergy(PickSite, true); // DWDWDW
            
            // 2.0 times because neighbors feel the same interaction
            // energy change as the adsorbate under consideration 
            double delH = 2.0 * (InitialEnergy - FinalEnergy);
            
            double Keq = exp(-delH * theRT);
            
            if(isnan(Keq))
            {
                Keq = double_limit;
            }
            
            if(Keq >= 1.0)
            {
                // do nothing, orientation accepted
            }
            else
            {
                Prob = FindRandomNumber();
                if(Prob < Keq)
                {
                    // do nothing, orientation accepted
                }
                else
                {
                    PickSite->Orientation = MyOr;
                    //PickSite->BindingEnergy = InitialEnergy;
                }
            }
        }
    }
    return 1;
}

int Simulation::Populate(int *Difference, fstream &fout)
{
    // Calculate the Number Of Atoms to add or subtract Surface Population
    int *SurfaceAtoms = new int[NumberOfSpecies];
    
    assert(SurfaceAtoms);
    int Total;
    int ind = FindChargeIndex(SurfaceCharge);
    theGrid->GetSurfaceAtoms(SurfaceAtoms, Total);
    
    thisComponent = theNULLSpecies;
    
    for(int i = 0; i < NumberOfSpecies; ++i)
    {
        thisComponent = Listing[i];
        Difference[i] = Difference[i] - SurfaceAtoms[i];
        if(SurfaceAtoms[i] == 0 && Difference[i] <= 0)
            Difference[i] = 0;
    }
    GridSite *AnySite = NULL;
    
    bool Converged = false; // DWDWDW // DWDWDW
    
    double T = 100.0;
    IRandom SpeciesType;
    
    bool Continue = false; // DWDWDW // DWDWDW
    
    for(int i2 = 0; i2 < NumberOfSpecies; ++i2)
    {
        if(Difference[i2] != 0)
        {
            Continue = true; // DWDWDW
        }
    }
    if(Continue == false) // DWDWDW
    {
        Converged = true; // DWDWDW
    }
    double DeltaE, BE;
    Random Ran;
    
    Ran.SetInterval(0., 1.);
    double RanNum;
    double Probability;
    int ToLong = false; // DWDWDW
    int OLDT = 0;
    int MaxCycles = 1000 * theGrid->EDim[1] * theGrid->EDim[1];
    int Cycles = MaxCycles;
    
    thisSite = NULL;
    thisReaction = NULL;
    int niter = 0;
    while(!Converged && Cycles--)
    {
        niter++;
        //if(Cycles % 9999 == 0)
           // SuperFastDiffuseEquilibration(0.1);
        Continue = false; // DWDWDW
        
        for(int i = 0; i < NumberOfSpecies; ++i)
        {
            if(Difference[i] != 0)
            {
                Continue = true; // DWDWDW
            }
        }
        if(Continue == false)
        {
            Converged = true;
        }
        AnySite = theGrid->PickRandom();
        
        // DW 8/24/04 commented out i, moved inside if block
        //i = 0;
        if(AnySite->GetType() != theNULLSpecies)
        {
            int i = 0;
            for(int j = 0; j < NumberOfSpecies; ++j)
            {
                if(Difference[j] < 0)
                {
                    i = j;
                }
            }
            if(Difference[i] < 0 && AnySite->GetType() == Listing[i])
            {
                AnySite->SetType(theNULLSpecies);
                Difference[i] += 1;
            }
        }
        else
        {
            //RD: Rather than drawing from a list of all species, 
            //draw only from the list of species whose Difference > 0
            int List[NumberOfSpecies];
            int NSp = 0;
            for (int Sp = 0; Sp < NumberOfSpecies; Sp++){
                if (Difference[Sp]>0){
                    List[NSp] = Sp;
                    ++NSp;
                }
            }
	          double x = drand48()*NSp;
            int i = x;
            i = List[i];
            bool Okdokey = false; // DWDWDW  // DWDWDW
            
            for(int j = 0; j < NumberOfSpecies; ++j)
            {
                if(Difference[j] > 0)
                {
                    Okdokey = true; // DWDWDW
                    break;
                }
            }
            if(Okdokey)
            {
                Component *PopulateComponent = Listing[i];
                double InitEnergy = theModel->CalculateBindingEnergy(AnySite);
                for (int j = 0; j < AnySite->SurroundSize; j++)
                {
                  InitEnergy += theModel->CalculateBindingEnergy(AnySite->Surround[j]);
                }
                AnySite->SetType(PopulateComponent);
                double FinalEnergy = theModel->CalculateBindingEnergy(AnySite);
                for (int j = 0; j < AnySite->SurroundSize; j++)
                {
                  FinalEnergy += theModel->CalculateBindingEnergy(AnySite->Surround[j]);
                }
                BE = FinalEnergy - InitEnergy;
                double dummy = theModel->CalculateBindingEnergy(AnySite);
                AnySite->SetType(theNULLSpecies);
                if(dummy > 1.0e-4)
                {
                    DeltaE = BE; // - 0.1*PopulateComponent->StableSiteEnergy[ind];
                    RanNum = Ran.Draw();
                    Probability = exp(DeltaE / (R_GasConst * T));
                    T += 100. / (double)theGrid->GetNumMetalAtoms();
                    if(T > 1500.0)
                    {
                        T = 1500.0;
                    }
                    if(Debug && (int)T != OLDT)
                    {
                        OLDT = (int)T;
                        //cout << "\nInit : Temperature = " << T;
                    }
                    if(RanNum < Probability)
                    {
                        if(Debug)
                            //cout << "\n\t\tNEW EVENT" << endl;
                        Cycles = MaxCycles;
                        AnySite->SetType(PopulateComponent);
                        Difference[i] -= 1;
                        //SuperFastDiffuseEquilibration(0.025);
                        EPosition <double> a;
                        a = AnySite->GetPosition();
                        //cout << a.GetXDistance() << " " << a.GetYDistance() << " " << a.GetZDistance() << endl;
                    }
                    
                }
                
            }
            
        }
        Continue = false; // DWDWDW
        for(int i2 = 0; i2 < NumberOfSpecies; ++i2)
        {
            if(Difference[i2] != 0)
            {
                Continue = true; // DWDWDW
            }
        }
        if(Continue == false) // DWDWDW
        {
            Converged = true; // DWDWDW
        }
    }
    
    for(int i3 = 0; i3 < NumberOfSpecies; ++i3)
    {
        SurfaceAtoms[i3] = 0;
    }
    
    theGrid->GetSurfaceAtoms(SurfaceAtoms, Total);
    for(int i4 = 0; i4 < NumberOfSpecies; ++i4)
    {
        Listing[i4]->NumberOfAtoms = SurfaceAtoms[i4];
    }
    // Output the coverage of species from the data stored in theGrid
    int SurfAtom[NumberOfSpecies];
    int Tot;
    theGrid->GetSurfaceAtoms(SurfAtom, Tot);
    for(int i4 = 0; i4 < NumberOfSpecies; ++i4)
    {
        Listing[i4]->NumberOfAtoms = SurfaceAtoms[i4]; 
    }

    fout << endl;
    delete[]SurfaceAtoms;
    return ToLong;
}

void Simulation::SimulationSummary()
{
    int ind = theSimulation->FindChargeIndex(SurfaceCharge);
    
    if(Debug)
    {
        cout << "\n\n************ SIMULATION SUMMARY ************\n\n\n\n"
            << endl;
        cout.setf(ios::scientific);
        cout << setw(18) << "FAST REACTIONS";
        cout << setw(18) << "RATE";
        cout << setw(18) << "~ACTIVATION~";
        cout << endl;
        for(int i = 0; i < NumberOfSimReactions; ++i)
        {
            thisReaction = SimReactions[i];
            cout << "\n" << thisReaction->GetName();
            cout << setw(18) << setprecision(4) << thisReaction->Rate;
            double Val;
            double RT = R_GasConst * Temperature;

            switch (thisReaction->GetProcess())
            {
                case Ad1:
                case Ad2:
                case EleyAd:
                case EleyAdEleyDes:
                    Val = -RT * log(thisReaction->MacroscopicBarrier
                        / (thisReaction->NumberOfScenarios *
                            thisReaction->GetPreexponential() *
                            thisReaction->Reactant[0]->
                            CalculateAdsorptionRate(GetTemperature())));
                    break;
                default:
                    Val = -RT * log(thisReaction->MacroscopicBarrier
                        / (thisReaction->NumberOfScenarios *
                            thisReaction->GetPreexponential()));
                    break;
            }
            cout << setw(18) << Val << flush;
        }
        cout << "\n\n\n";
        cout << setw(Rxns->GetPtr()->GetName().
            length()) << "SLOW REACTIONS";
        cout << setw(18) << "BULK RATE";
        cout << setw(18) << "ADJ. RATE";
        cout << setw(18) << "TIME AVG. RATE";
        cout << setw(18) << "~ACTIVATION~";
        cout << endl;
        cout << setw(26) << "\nBINDING ENERGIES\n";
        for(int i2 = 0; i2 < NumberOfSpecies; ++i2)
        {
            thisComponent = Listing[i2];
            if(thisComponent->SimBindingEnergy > 0.0)
            {
                cout << setw(18) << thisComponent->GetName() << "\t";
                cout << setw(12) << thisComponent->SimBindingEnergy << "\n";
            }
            else if(theGrid->GetCoverage(thisComponent) > 1.e-36)
            {
                cout << setw(18) << thisComponent->GetName() << "\t";
                cout << setw(12) << thisComponent->SimBindingEnergy;
                cout << " ERROR -> NEGATIVE BINDING ENERGY!\n";
            }
        }
        cout << endl;
        cout << "Elapsed time is: " << Time << endl;
    }
    return;
}
