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

//////////////////////////////////////////////////////////
// "Main" program
//////////////////////////////////////////////////////////

#ifdef GUI

#include "Motif/Application.h"
#include "ReactionWindow.h"
#include "SpeciesWindow.h"
#include "InteractionWindow.h"
#include "GridWindow.h"
#include "ModelWindow.h"
#include "ScenarioWindow.h"
#include "SimulationWindow.h"

Application myApp("CIRE");

SpeciesWindow *theSpeciesWindow = new SpeciesWindow( "SpeciesWindow" );
 
ReactionWindow *theReactionWindow = new ReactionWindow( "ReactionWindow" );

InteractionWindow *theInteractionWindow =
    new InteractionWindow( "InteractionWindow" );

GridWindow *theGridWindow = new GridWindow( "GridWindow" );

ModelWindow *theModelWindow = new ModelWindow( "ModelWindow" );

ScenarioWindow *theScenarioWindow = new ScenarioWindow( "ScenarioWindow" );

SimulationWindow *theSimulationWindow =
    new SimulationWindow( "SimulationWindow" );

#endif

// void Redisplay() {
// #ifdef GUI
// 	theGridWindow->Redisplay();
// #endif
// 	return;
// }
