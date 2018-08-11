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
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <UpdaterAddCosolvent.h>
#include <LeMonADE/feature/FeatureAttributes.h> 
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/feature/FeatureLattice.h>




#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include "catchorg/clara/clara.hpp"

int main(int argc, char* argv[]){


try{
	uint32_t N  = 100;
	uint32_t max_mcs=100;
	uint32_t save_interval=100;
        double eps =0.0;
    bool showHelp = false;

    auto parser
    = clara::Opt( N, "polymer length" )
        ["-n"]["--len"]
        ("Length of the polymer")
        .required()
    | clara::Opt( eps, "interaction energy" )
        ["-e"]["--eps"]
        ("ineraction energy between the monomers in 1/1000*kT")
        .required()
    | clara::Opt( [&max_mcs](int const m)
        {
         if (m < 0)
         {
            return clara::ParserResult::runtimeError("Simulation time must be greater than 0");
         }
         else
         {
            max_mcs = m;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "max MCS(=100)" )
        ["-m"]["--max-mcs"]
        ("(required) specifies the total Monte-Carlo steps to simulate.")
        .required()
    | clara::Opt( [&save_interval](int const s)
        {
         if (s < 0)
         {
            return clara::ParserResult::runtimeError("Save intervall must be greater than 0");
         }
         else
         {
            save_interval = s;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
      }, "save MCS(=100)")
        ["-s"]["--save-mcs"]
        ("(required) Save after every <integer> Monte-Carlo steps to the output file." )
        .required()
    | clara::Help( showHelp );

    auto result = parser.parse( clara::Args( argc, argv ) );


    if( !result ) {
    std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
    exit(1);
    }
    else if(showHelp == true)
    {
        std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck" << std::endl
                  << "maximum number of connections per monomer is 6" << std::endl
                  << "Features used: FeatureBondset, FeatureExcludedVolumeSc<FeatureLattice<bool> >" << std::endl
		          << "Updaters used: ReadFullBFMFile, SimpleSimulator" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;

        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {




    //first set up the random number generator
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    //now set up the system
    //here, we use one additional feature called FeatureMoleculesIO
    //strictly speaking this is not a feature in the sense that it
    //changes the simulation conditions. What it provides is the
    //basic functionalities for writing BFM files. So, in most cases
    //you are going to use this feature in simulations.
//    typedef LOKI_TYPELIST_3(FeatureMoleculesIO,
//                            FeatureAttributes, FeatureExcludedVolumeSc<>) Features; 
    typedef LOKI_TYPELIST_5(FeatureMoleculesIO,FeatureBox,FeatureBondset<>,
                            FeatureAttributes, FeatureNNInteractionSc<FeatureLattice>) Features;

    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;


  
    mySystem.setBoxX(128);
    mySystem.setBoxY(128);
    mySystem.setBoxZ(128);
    


    mySystem.setPeriodicX(true);
    mySystem.setPeriodicY(true);
    mySystem.setPeriodicZ(true);

    mySystem.modifyBondset().addBFMclassicBondset();

    mySystem.synchronize(mySystem);
    
    mySystem.setNNInteraction(1,1,eps/1000.0);


//    //Setting up Linear chains by hand:
//    mySystem.modifyMolecules().resize(N);
//
//    //add the polymer chain
//    mySystem.modifyMolecules()[0].setX(0);
//    mySystem.modifyMolecules()[0].setY(0);
//    mySystem.modifyMolecules()[0].setZ(0);
//
//    for(uint32_t i=1;i<N;i++){
//	    mySystem.modifyMolecules()[i].setX(i*2);
//	    mySystem.modifyMolecules()[i].setY(i*2);
//	    mySystem.modifyMolecules()[i].setZ(i);
//
//	    mySystem.modifyMolecules().connect(i-1,i);
//    }




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


    std::ostringstream filename;
    filename << "LinearChain_" << N << "_" << max_mcs << "_" << save_interval << ".bfm";

    std::cout     << "outputfile:       " << filename.str() << std::endl
                  << "max_mcs:       " << max_mcs << std::endl
                  << "save_interval: " << save_interval << std::endl;
    //create the task manager
    TaskManager taskmanager;


    //updater to add Linear Chain:

    taskmanager.addUpdater(new
    UpdaterAddLinearChains<MyIngredients>(mySystem,1,N,1,1),0);

    //updater to add cosolvent:

    //taskmanager.addUpdater(new
    //UpdaterAddCosolvent<MyIngredients>(mySystem,100,5),0);
        //add the simulator
    //the syntax says: the simulator should simulate mySystem
    //every time it is executed, it simulates for 10000 steps
    //The 1 as second argument to addUpdater says that the
    //simulation is to be called in every circle.
    taskmanager.addUpdater(new
    UpdaterSimpleSimulator<MyIngredients,MoveLocalSc>(mySystem,save_interval),1);

    //add the file output, the trajectory file name will be
    // "LinearChain_<ChainLength>.bfm" as defined above.

    taskmanager.addAnalyzer(new
    AnalyzerWriteBfmFile<MyIngredients>(filename.str(),mySystem),1);



    taskmanager.initialize(); 
    if(max_mcs > 0){
    taskmanager.run(max_mcs/save_interval);}
    if(max_mcs == 0){
    taskmanager.run(0);}
    taskmanager.cleanup();

    /* **************************************************************
      * The program should produce among others a file called
      * polymer.bfm. This contains the trajectory. You can visualize
      * the trajectory using the provided LemonadeViewerFLTK from
      * the projects folder.
      * *************************************************************/

}}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
}







