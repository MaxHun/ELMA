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
* The aim of this program is to simulate a linear chain surrounded by a
* (good) cosolvent of certain concentration,, perforing a NN-interaction
* with the polymer of some strength eps.
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
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <UpdaterAddCosolvent.h>

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include "catchorg/clara/clara.hpp"

int main(int argc, char* argv[]){


try{
    /*
    The parameters that can be specified when running the program ought to be:
    chain length N
    Boxsize L
    Number of cosolvent monomers Ncos
    strength of interaction e
    maximum MCS (=1e7)
    save MSC (=5e4)
    */
	uint32_t N  = 256;
	uint32_t max_mcs=100000;
	uint32_t save_interval=500;
	uint32_t eps = 0;
	uint32_t NCos = 1000;
	double   c = 1;// concenttration of the Cosolvent in 1/1000
	uint32_t L = 256;

    bool showHelp = false;

    auto parser

    = clara::Opt( N, "polymer length" )
        ["-n"]["--len"]
        ("Length of the polymer")
        .required()
    | clara::Opt( [&max_mcs](int const m)
        {
         if (m <= 0)
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
    | clara::Opt(L, "box size" )
        ["-l"]["--bsize"]
        ("size of the quadratic box")
        .required()
    | clara::Opt(c, "cosolvent concentration" )
        ["-c"]["--con"]
        ("concentration of the cosolvent in 1/1000")
        .required()
    | clara::Opt(eps, "interaction strength" )
        ["-e"]["--eps"]
        ("interaction energy for polymer-cosolvent interaction in units of 0.001*k*T")
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

    //std::cout << eps << std::endl;


    NCos = (int)(((c/1000)*L*L*L - 8*N)/8); // Calculate NCos, c gets devided by 1000


    std::cout << NCos << std::endl;



    //first set up the random number generator

    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    //now set up the system

    typedef LOKI_TYPELIST_5(FeatureMoleculesIO,FeatureBox,FeatureBondset<>,
                            FeatureAttributes, FeatureNNInteractionSc<FeatureLattice>) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;



    mySystem.setBoxX(L);
    mySystem.setBoxY(L);
    mySystem.setBoxZ(L);


    mySystem.setPeriodicX(true);
    mySystem.setPeriodicY(true);
    mySystem.setPeriodicZ(true);

    // Checking periodicity for filename:
    int perX = mySystem.isPeriodicX();
    int perY = mySystem.isPeriodicY();
    int perZ = mySystem.isPeriodicZ();

    mySystem.modifyBondset().addBFMclassicBondset();

    //set the NNInteraction. The interaction energy is required in units of k*T
    mySystem.setNNInteraction(1,3,-(eps/1000.0));

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
    filename << "CosolvedLC_N" << N << "_c" << c << "_e" << eps << "Per" << perX << perY << perZ << ".bfm";

    std::cout     << "outputfile:       " << filename.str() << std::endl
                  << "max_mcs:       " << max_mcs << std::endl
                  << "save_interval: " << save_interval << std::endl
                  << "cosolvent concentration (1/1000): " << c << std::endl
                  << "intteraction strength (1/1000*kT):" << eps << std::endl
                  << "boxsize: " << L << std::endl
                  << "periodicity(XYZ): " << perX << perY << perZ << std::endl;



    // Print the value of eps:
//    std::cout     << "EPSILON BETRAEGT:"
//                  << eps << "\n"
//                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    //create the task manager
    TaskManager taskmanager;


    //updater to add Linear Chain:

    taskmanager.addUpdater(new
    UpdaterAddLinearChains<MyIngredients>(mySystem,1,N,1,1),0);

    taskmanager.addUpdater(new
    UpdaterAddCosolvent<MyIngredients>(mySystem,NCos,3),0);

    //add the simulator

    taskmanager.addUpdater(new
    UpdaterSimpleSimulator<MyIngredients,MoveLocalSc>(mySystem,save_interval),1);

    //add the file output, the trajectory file name will be
    // "CosolvedLC_<ChainLength>_<NCos>_<eps*1000>.bfm" as defined above.

    taskmanager.addAnalyzer(new
    AnalyzerWriteBfmFile<MyIngredients>(filename.str(),mySystem),1);



    taskmanager.initialize();
    taskmanager.run(max_mcs/save_interval);
    taskmanager.cleanup();



}}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
}







