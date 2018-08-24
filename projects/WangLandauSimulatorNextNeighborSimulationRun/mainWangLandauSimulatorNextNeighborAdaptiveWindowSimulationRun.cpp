

#include <cstring>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include "catchorg/clara/clara.hpp"

#include "FeatureWangLandauNextNeighborAdaptiveWindowSimulationRun.h"
#include "UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun.h"
#include "ReadInHGLnDOS.h"



int main(int argc, char* argv[])
{
	try{
		std::string infile  = "input.bfm";
		std::string outfile = "";
		std::string HGLnDOSfile = "";

		uint64_t max_mcs=100;
		uint32_t save_interval=100;
		uint32_t min_statistic_entries=100;
		double min_histogram = -100.0;
		double max_histogram = +100.0;
		uint32_t bins_histogram = 200;
		double modFactor = 1.01;

		double minWin = -100.0;
		double maxWin = +100.0;

		bool showHelp = false;

		bool noFileWriting = true;

		auto parser
		= clara::Opt( infile, "input (=input.bfm)" )
		["-i"]["--infile"]
			   ("BFM-file to load.")
			   .required()
		| clara::Opt( outfile, "output (=outfile.bfm)" )
		["-o"]["--outfile"]
			   ("BFM-file to save.")
			   .required()
		| clara::Opt(  min_histogram, "min histogram (=-100.0)" )
		["--min"]
			   ("minimal histogram boundary (=-100.0)")
			   .required()
		| clara::Opt(  max_histogram, "max histogram (=+100.0)" )
		["--max"]
			   	("maximal histogram boundary (=+100.0)")
			   	.required()
		| clara::Opt(  bins_histogram, "bins histogram (=+200.0)" )
		["--bins"]
				("bins histogram boundary (=+200.0)")
				.required()
		| clara::Opt( [&max_mcs](uint64_t const m)
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
		| clara::Opt( [&min_statistic_entries](int const b)
			   	{
			   					if (b < 0)
			   					{
			   						return clara::ParserResult::runtimeError("min_statistic_entries must be greater than 0");
			   					}
			   					else
			   					{
			   						min_statistic_entries = b;
			   						return clara::ParserResult::ok(clara::ParseResultType::Matched);
			   					}
			   	}, "min_statistic_entries (=100)")
			   		["--min-statistic"]
			   			   ("(required) <integer> Number of entries for calculating the average at specific energy." )
			   			   .required()
		| clara::Opt(  minWin, "min window (=-100.0)" )
			["--min-win"]
			 ("minimal window boundary (=-100.0)")
			 .required()
		| clara::Opt(  maxWin, "max window (=+100.0)" )
			["--max-win"]
			("maximal window boundary (=+100.0)")
			.required()
		| clara::Opt( HGLnDOSfile, "HGLnDOS" )
			["--HGLnDOS"]
			   ("Histogram file for (logarithmic) DOS to load.")
			   .required()
		| clara::Help( showHelp );

		auto result = parser.parse( clara::Args( argc, argv ) );
		if( !result ) {
			std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
			exit(1);
		}
		else if(showHelp == true)
		{
			std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck and WangLandau in NN-shell" << std::endl
					<< "maximum number of connections per monomer is 6" << std::endl
					<< "Features used: FeatureBondset, FeatureAttributes, FeatureExcludedVolumeSc<FeatureLattice<uint8_t> >, FeatureWangLandauNextNeighbor" << std::endl
					<< "Updaters used: ReadFullBFMFile, SimpleSimulator" << std::endl
					<< "Analyzers used: WriteBfmFile" << std::endl;

			parser.writeToStream(std::cout);
			exit(0);
		}
		else
		{
			std::cout << "infile:        " << infile << std::endl
					<< "outfile:       " << outfile << std::endl
					<< "max_mcs:       " << max_mcs << std::endl
					<< "save_interval: " << save_interval << std::endl
					<< "min_statistic_entries: " << min_statistic_entries << std::endl
					<< "min_histogram: " << min_histogram << std::endl
					<< "max_histogram: " << max_histogram << std::endl
					<< "bins_histogram: " << bins_histogram << std::endl
					<< "min_win: "	<< minWin << std::endl
					<< "max_win: " 	<< maxWin << std::endl
					<< "HGLnDOS:" << HGLnDOSfile << std::endl;
					;
		}
	
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureNNInteractionSc< FeatureLattice >) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureWangLandauNextNeighborAdaptiveWindowSimulationRun< FeatureLattice >) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureWangLandau, FeatureAttributes) Features;

	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	//read-in the histogram of the DOS
	if(HGLnDOSfile.empty())
		std::runtime_error("File HGLnDOS has to be provided. EXITing...\n");


	ReadInHGLnDOS in(min_histogram, max_histogram, bins_histogram, HGLnDOSfile);
	in.readin();

	//copy HGLnDOS into ingredients
	myIngredients.modifyHGLnDOS().reset(min_histogram, max_histogram, bins_histogram);

	for(size_t n=0;n<in.getHGLnDOS().getVectorValues().size();n++){

		if(in.getHGLnDOS().getVectorValues()[n].ReturnN() != 0)
				myIngredients.modifyHGLnDOS().addValue(in.getHGLnDOS().getCenterOfBin(n), in.getHGLnDOS().getVectorValues()[n].ReturnM1());
	}


	// run the simulation and gather the information

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));
	
	taskmanager.addUpdater(new UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<Ing,MoveLocalSc>(myIngredients,
							   minWin, maxWin, min_statistic_entries, max_mcs, save_interval),
					       1);

	if(!outfile.empty())
		taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));


	taskmanager.initialize();
	taskmanager.run();
	taskmanager.cleanup();
	
	outfile=infile+"_final";

	AnalyzerWriteBfmFile<Ing> ABFM(outfile,myIngredients);
	ABFM.initialize();
	ABFM.execute();
	ABFM.cleanup();

	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

