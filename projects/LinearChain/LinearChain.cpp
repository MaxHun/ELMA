/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Hauke Rabbel)
    ooo                        | 
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

/* *********************************************************************
* The aim of this program is to simulate a linear chain in a periodic box 
* and to produce a .bfm-File that can be used for further analysis.
* *********************************************************************/

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

// New with respect to Example4:
#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char* argv[]){
	

    //first set up the random number generator
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();
    
    //now set up the system
    //here, we use one additional feature called FeatureMoleculesIO
    //strictly speaking this is not a feature in the sense that it 
    //changes the simulation conditions. What it provides is the
    //basic functionalities for writing BFM files. So, in most cases
    //you are going to use this feature in simulations.
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO,FeatureBox,FeatureBondset<>) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;

    mySystem.setBoxX(64);
    mySystem.setBoxY(64);
    mySystem.setBoxZ(64);
    
    mySystem.setPeriodicX(true);
    mySystem.setPeriodicY(true);
    mySystem.setPeriodicZ(true);
    
    mySystem.modifyBondset().addBFMclassicBondset();

    mySystem.synchronize(mySystem);
    
    
    //and add a polymer chain of length N
    uint32_t N=10;
    
    mySystem.modifyMolecules().resize(N);
    
    //add the polymer chain
    mySystem.modifyMolecules()[0].setX(0);
    mySystem.modifyMolecules()[0].setY(0);
    mySystem.modifyMolecules()[0].setZ(0);
    
    for(uint32_t i=1;i<N;i++){
	    mySystem.modifyMolecules()[i].setX(i*2);
	    mySystem.modifyMolecules()[i].setY(i*2);	
	    mySystem.modifyMolecules()[i].setZ(i);
	    
	    mySystem.modifyMolecules().connect(i-1,i);
    }
    
    
    mySystem.synchronize(mySystem);
    
    /* ****************************************************************
      * Now we can set up the task manager with the desired tasks.
      * We want to do two things:
      * 1) Simulate the system. For this we use the updater 
      * UpdaterSimpleSimulator. The source code of this updater
      * can be found in src/updater/UpdaterSimpleSimulator.h
      * 2) Write the trajectory to a file. For this we use the analyzer 
      * AnalyzerWriteBfmFile. The source code of this analyzer can be
      * found in src/analyzer/AnalyzerWriteBfmFile.h
      * ***************************************************************/
    
    //create the task manager
    TaskManager taskmanager;
    
    //add the simulator
    //the syntax says: the simulator should simulate mySystem
    //every time it is executed, it simulates for 10000 steps
    //The 1 as second argument to addUpdater says that the 
    //simulation is to be called in every circle.
    taskmanager.addUpdater(new 
    UpdaterSimpleSimulator<MyIngredients,MoveLocalSc>(mySystem,50000),1);
    
    //add the file output, the trajectory file name will be 
    // "LinearChain_<ChainLength>.bfm"
    std::ostringstream filename;
    filename << "LinearChain_" << N << ".bfm"; 
    
    
    taskmanager.addAnalyzer(new 
    AnalyzerWriteBfmFile<MyIngredients>(filename.str(),mySystem),1);
    
    
    
    taskmanager.initialize();
    taskmanager.run(2000);
    taskmanager.cleanup();
    
    /* **************************************************************
      * The program should produce among others a file called
      * polymer.bfm. This contains the trajectory. You can visualize
      * the trajectory using the provided LemonadeViewerFLTK from
      * the projects folder.
      * *************************************************************/
    
    return 0;
}
