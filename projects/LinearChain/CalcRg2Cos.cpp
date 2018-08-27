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
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <Analyzer_ChainWalking_RG2.h>

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include "catchorg/clara/clara.hpp"

int main(int argc, char* argv[]){


try{
	std::string infile  = "input.bfm";
        int t = 0;


    bool showHelp = false;

    auto parser
     =clara::Opt( infile, "input (=input.bfm)" )
        ["-i"]["--infile"]
        ("BFM-file to load.")
        .required() 
    | clara::Opt(t, "evaluation time")
        ["-e"]["--evt"]
        ("First <integer> MCS will not be considered for calculation of Rg2")
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

    //now set up the system:

    typedef LOKI_TYPELIST_5(FeatureAttributes,FeatureMoleculesIO,FeatureBox,FeatureBondset<>,FeatureExcludedVolumeSc<>) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;

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

    //define the name of the output .dat file:
    std::ostringstream filename;
    filename << infile.substr(0,infile.size()-4) << "_Rg2.dat";

    std::cout     << "inputfile:        " << filename.str() << std::endl
                  << "outputfile:       " << infile << std::endl;

    //create the task manager
    TaskManager taskmanager;



    //add the file output, the trajectory file name will be
    // "LinearChain_<ChainLength>.bfm" as defined above.

    taskmanager.addUpdater(new UpdaterReadBfmFile<MyIngredients>(infile,mySystem,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE),1);
    taskmanager.addAnalyzer(new
    Analyzer_ChainWalking_RG2<MyIngredients>(mySystem,t),1);


    taskmanager.initialize();
    taskmanager.run();
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






