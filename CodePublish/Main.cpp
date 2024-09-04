///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//         This example code is from the book:
//
//           Object-Oriented Programming with C++ and OSF/Motif
//         by
//           Douglas Young
//           Prentice Hall, 1992
//           ISBN 0-13-630252-1	
//
//         Copyright 1991 by Prentice Hall
//         All Rights Reserved
//
//  Permission to use, copy, modify, and distribute this software for 
//  any purpose except publication and without fee is hereby granted, provided 
//  that the above copyright notice appear in all copies of the software.
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
using namespace std;



//////////////////////////////////////////////////////////
// Main.C: Generic main program used by all applications
//////////////////////////////////////////////////////////
#ifdef GUI
  #include "Motif/Application.h"
  #include "GridWindow.h"
#endif


// We can implement main() in the library because the 
// framework completely encapsulates all Xt boilerplate 
// and all central flow of control. 
#include "rank.h"
#include "Species.h"
#include "Reactions.h"
#include "Interactions.h"
#include "Grid.h"
#include "Model.h"
#include "ScenarioBuilder.h"
#include "Simulation.h"
#include "Base/Constants.h"

//Including cpp files
#include "Species.cpp"
#include "Reactions.cpp"
#include "Interactions.cpp"
#include "Grid.cpp"
#include "Model.cpp"
#include "ScenarioBuilder.cpp"
#include "Simulation.cpp"
#include "backprop.cpp"
#include "Component.cpp"
#include "ElementaryRxn.cpp"
#include "Geometry.cpp"
#include "GridSite.cpp"
#include "layer.cpp"
#include "Scenario.cpp"
#include "SimSub.cpp"
#include "ShortRoutines.cpp"
#include "Simulation_IO.cpp"
#include "Base/Color.cpp"
#include "Base/Estring.cpp"
#include "Base/Facilitator.cpp"
#include "Base/IntRandom.cpp"
#include "Base/Matrix.cpp"
#include "Base/Position.cpp"
#include "Base/Random.cpp"

Component* thisComponent;
const char *rnk;
int rank;
unsigned int rank2;
int main(int argc, char *argv[])
{ 	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int sim_nmb = 167;
    rank2 = sim_nmb*12 + rank;
    MPI_Finalize();
    stringstream strs;
    strs << rank;
    string temp_str = strs.str();
    rnk = (char*) temp_str.c_str();
    fstream fin;
    fstream fout;
    char output_file[100];
    strcpy(output_file, argv[2]);
    strcat(output_file, rnk);
    fin.open(argv[1], ios::in);
    fout.open(output_file, fstream::out);
    theSpecies = new Species(fin, fout);
    cout << "Species done" << endl;
    theReactions = new Reactions(fin, fout);
    cout << "Reactions done" << endl;
    theInteractions = new Interactions(fin, fout);
    cout << "Interactions done" << endl;
    theGrid = new Grid(fin, fout);
    cout << "Grid done" << endl;
    theModel = new Model(fin, fout);
    cout << "Model done" << endl;
    theScenarioBuilder = new ScenarioBuilder(fin, fout);
    cout << "Builder done" << endl;
    theSimulation = new Simulation(fin, fout);
    cout << "Simulation done" << endl;
    theSimulation->Simulate(fout);
    return 0;
}
