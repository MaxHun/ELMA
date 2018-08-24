#ifndef ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR_AdaptiveWindow_SimulationRun_H
#define ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR_AdaptiveWindow_SimulationRun_H

#include <limits>
#include <iomanip>

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/UpdaterAbstractCreate.h>

#include "Histogram1D.h"
#include "HistogramGeneralStatistik1D.h"

template<class IngredientsType,class MoveType>
class UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun: public UpdaterAbstractCreate<IngredientsType>
{
	
public:

	 typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
	/**
	 * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
	 *
	 * @param ing a reference to the IngredientsType - mainly the system
	 * @param steps MCS per cycle to performed by execute()
	 * @param stepsBeforeBiasCheck interval between checks if histogram has converged im MCS
	 * @param stepsBeforeHistogramUpdate interval between updates on coordinate histogram in MCS
	 * @param convergenceThreshold threshold for relative convergence of histogram
	 * @param fullFlatnessThreshold threshold for absolute convergence of histogram
	 * @param maxAge maximum MCS limit. If the histogram converges earlier, the simulation also stops
	 */
	UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun(IngredientsType& ing
	,double minWindow = -100.0
	,double maxWindow = +100.0
	,uint32_t min_statistic_entries = 100
	,uint64_t maxAge=100
	,uint32_t steps=100
	);
	
	/**
	 * @brief This checks all used Feature and applies all Feature if all conditions are met.
	 *
	 * @details This function runs over \a steps MCS and performs the moves.
	 * It setting the age of the system and prints a simple simple simulation speed
	 * in the number of attempted monomer moves per s (tried and performed monomer moves/s).
	 *
	 * @return True if function is done.
	 */
	bool execute();
	
	
	/**
	 * @brief This function is called \a once in the beginning of the TaskManager.
	 *
	 * @details It´s a virtual function for inheritance.
	 * Use this function for initializing tasks (e.g. init SDL)
	 *
	 **/
	virtual void initialize();
	
	/**
	 * @brief This function is called \a once in the end of the TaskManager.
	 *
	 * @details writes the final state of histogram, bias potential and the
	 * history of the convergence criteria to output files.
	 *
	 **/
	virtual void cleanup(){


		dumpEnergyRG2HGLnDOS(std::string(ingredients.getName() + prefixWindow + "_final_E_HGLnDOS_RG2.dat"), ingredients.getMinWin(), ingredients.getMaxWin());

	};
	
	
private:
	using BaseClass::ingredients;

	//! Specialized move to be used
	MoveType move;

	double calcRG2();

	
	//! write the currently logarithm of Density of States DOS i.e. file named HGLnDOS.dat
	void dumpEnergyRG2HGLnDOS(std::string prefix="", double min= std::numeric_limits<double>::min(), double max= std::numeric_limits<double>::max());

	
	//! flags for setting or unsetting extra output of histogram
	bool writeHistogramProgress;
	
	//! flags for setting or unsetting extra output of bias potential
	bool writePMFProgress;
	
	double minWindow;
	double maxWindow;
	std::string prefixWindow;


	HistogramGeneralStatistik1D HG_Energy_Rg2;

	uint32_t minStatisticEntries;

	//! Max number of mcs to be executed
	uint32_t nsteps;

	//! maximum age of the system
	uint64_t maxSystemAge;

};


////////////////////////////////////////////////////////////////////////////////
// member implementations
///////////////////////////////////////////////////////////////////////////////


// constructor
template<class IngredientsType,class MoveType>
UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<IngredientsType,MoveType>::
UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun(IngredientsType& ing,
				 double _minWindow
				,double _maxWindow
				,uint32_t min_statistic_entries
				,uint64_t maxAge
				,uint32_t steps
				)
:writeHistogramProgress(true)
,writePMFProgress(true)
,minWindow(_minWindow)
,maxWindow(_maxWindow)
,prefixWindow("")
,minStatisticEntries(min_statistic_entries)
,maxSystemAge(maxAge)
,nsteps(steps)
,BaseClass(ing)
{
	ingredients.setWindowState(false, ingredients.getHGLnDOS().getMinCoordinate(), ingredients.getHGLnDOS().getMaxCoordinate());
	HG_Energy_Rg2.reset(ingredients.getHGLnDOS().getMinCoordinate(), ingredients.getHGLnDOS().getMaxCoordinate(), ingredients.getHGLnDOS().getNBins());
}
 
 
//execute simulation
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<IngredientsType,MoveType>::execute()
{
	//some output for simulation performance and state
	//time_t startTimer = std::time(NULL); //in seconds
	//uint64_t systemAge=ingredients.getMolecules().getAge();
	//std::cout<<"mcs "<<systemAge
	//		<< " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) <<std::endl;

	//if max age has reached, stop
	if(ingredients.getMolecules().getAge() >= maxSystemAge) {
		std::cout<<"max age has been reached. maxAge "<<maxSystemAge<<" system age "<<ingredients.getMolecules().getAge()<<std::endl;
		return false;
	}

	for(int n=0;n<nsteps;n++){

		//do one monte carlo sweep
				for(int m=0;m<ingredients.getMolecules().size();m++)
				{
					move.init(ingredients);

					if(move.check(ingredients)==true)
					{
						move.apply(ingredients);
					}

				}
	}
	
	//update age of system and give some more output on simulation progress
	ingredients.modifyMolecules().setAge(ingredients.modifyMolecules().getAge()+nsteps);
	//std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with " << (((1.0*nsteps)*ingredients.getMolecules().size())/(std::difftime(std::time(NULL), startTimer)) ) << " [attempted moves/s]" <<std::endl;
	//std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) << " with " << nsteps << " MCS "<<std::endl;
	//std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with energy " << ingredients.getInternalEnergyCurrentConfiguration(ingredients) <<std::endl;


		//std::cout << " done" << std::endl;
	if(writeHistogramProgress==true)
	{
		if((ingredients.getMolecules().getAge() % 1000000)==0)
		{
			std::cout<<"mcs "<<ingredients.getMolecules().getAge() << std::endl;

			 dumpEnergyRG2HGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_E_HGLnDOS_RG2.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
			 std::cout << "dump E_HGLnDOS_RG2 " << ingredients.getMinWin() << "   " << ingredients.getMaxWin() << std::endl;
		}
	}

		
		//if( (ingredients.isEnergyInWindow() == false ) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) < maxWindow+2*ingredients.getHGLnDOS().getBinwidth()) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) > minWindow-2*ingredients.getHGLnDOS().getBinwidth()) )
			if( (ingredients.isEnergyInWindow() == false ) && (ingredients.getInternalEnergyWithoutCalculation() < maxWindow+2*ingredients.getHGLnDOS().getBinwidth()) && (ingredients.getInternalEnergyWithoutCalculation() > minWindow-2*ingredients.getHGLnDOS().getBinwidth()) )
		{
			std::cout << "RW is in energy window: [" << minWindow << " ; " << maxWindow << "] with " <<  ingredients.getInternalEnergyCurrentConfiguration(ingredients) << std::endl;
			// RW is in energy window
			ingredients.setWindowState(true, minWindow, maxWindow);

			std::stringstream ssprefixWindow;
			ssprefixWindow << "_minWin" << minWindow << "_maxWin" <<  maxWindow;

			prefixWindow = ssprefixWindow.str();

			HG_Energy_Rg2.reset(HG_Energy_Rg2.getMinCoordinate(),HG_Energy_Rg2.getMaxCoordinate(),HG_Energy_Rg2.getNBins());

		}

	double RG2 = calcRG2();
	double energy = ingredients.getInternalEnergyWithoutCalculation();//getInternalEnergyCurrentConfiguration(ingredients);
	HG_Energy_Rg2.addValue(energy, RG2);
    
    size_t counts = 0;
    
    for(size_t n=0;n<HG_Energy_Rg2.getNBins();n++)
		{
			if( (HG_Energy_Rg2.getCenterOfBin(n) >= ingredients.getMinWin()) && (HG_Energy_Rg2.getCenterOfBin(n) <= ingredients.getMaxWin()) )
            {
                counts += HG_Energy_Rg2.getNumCountInBin(n);
               
            }
            
        }
        
    // check for minimal counts and avoid exit of the program
    if(counts < minStatisticEntries)
        return true;
    

	// check if all samples have their desired number of entries
	if( ingredients.isEnergyInWindow() == true )
	{
		bool isStatisticSufficient = true;

		for(size_t n=0;n<HG_Energy_Rg2.getNBins();n++)
		{
			if( (HG_Energy_Rg2.getCenterOfBin(n) >= ingredients.getMinWin()) && (HG_Energy_Rg2.getCenterOfBin(n) <= ingredients.getMaxWin()) )
            {
                // no counts at all in the bin
				//if(HG_Energy_Rg2.getNumCountInBin(n) == 0)
				//		return true;
                
                
				// some bins can not be filled
				if((HG_Energy_Rg2.getNumCountInBin(n) > 0) && (HG_Energy_Rg2.getNumCountInBin(n) < minStatisticEntries))
					{
							isStatisticSufficient = false;

							//std::cout << "need adding entries at bin" << n  << " at E=" << HG_Energy_Rg2.getCenterOfBin(n) << " with " << HG_Energy_Rg2.getNumCountInBin(n) << " / " <<  minStatisticEntries << std::endl;

							return true;
					}
            }
		}

		// minimum number of statistical entries reached - end simulation
		if(isStatisticSufficient == true)
		{
			std::cout << "Statistic reached desired number in energy window: [" << minWindow << " ; " << maxWindow << "] " << std::endl;

			return false;
		}


	}

	return true;
}



// initialize only sets the histogram and bias potential to 0
template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<IngredientsType,MoveType>::initialize()
{

}

template<class IngredientsType, class MoveType>
double UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<IngredientsType,MoveType>::calcRG2()
{
		//std::cout << "SimpleAnalyzer_Rg2.execute() at MCS:" << ingredients.getMolecules().getAge() << std::endl;

		int monomerCounter = 0;

		double Rg2 = 0.0;
		double Rg2_x = 0.0;
		double Rg2_y = 0.0;
		double Rg2_z = 0.0;

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l= k; l < ingredients.getMolecules().size(); l++)
				 if((ingredients.getMolecules()[k].getAttributeTag()==1) && (ingredients.getMolecules()[l].getAttributeTag()==1))
			{
				Rg2_x += (ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX())*(ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX());
				Rg2_y += (ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY())*(ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY());
				Rg2_z += (ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ())*(ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ());


			}

			if((ingredients.getMolecules()[k].getAttributeTag()==1))
				monomerCounter++;
		}
		/*if(monomerCounter != ingredients.getMolecules().size())
		{
			throw std::runtime_error("invalid number of monomers in COM-calculation");
		}*/

		/*if(monomerCounter != 256)
				{
					throw std::runtime_error("invalid number of monomers in COM-calculation");
				}
		*/

		Rg2_x /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_y /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_z /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());

		Rg2 = Rg2_x+Rg2_y+Rg2_z;

		return Rg2;
}


//write current bias potential to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighborAdaptiveWindowSimulationRun<IngredientsType,MoveType>::dumpEnergyRG2HGLnDOS(std::string prefix, double _min, double _max)
{
	std::stringstream filename;
	filename<< prefix;// << "_HGLnDOS" << "_mcs"<<ingredients.getMolecules().getAge()<<".dat";
	std::ofstream file(filename.str().c_str());

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	file << "# histogram HG_E_LnDOS: [" << ingredients.getHGLnDOS().getMinCoordinate() << " ; " << ingredients.getHGLnDOS().getMaxCoordinate() << " ; " << ingredients.getHGLnDOS().getNBins() << " ]" << std::endl;
	file << "# histogram HG_E_RG2: [" << HG_Energy_Rg2.getMinCoordinate() << " ; " << HG_Energy_Rg2.getMaxCoordinate() << " ; " << HG_Energy_Rg2.getNBins() << " ]" << std::endl;
	file << "# " << std::endl;
	file << "# energyE <LnDOS> <Rg2(E)> countsRg2" << std::endl;

	// fill the list
	for(size_t n=0;n<HG_Energy_Rg2.getNBins();n++)
		{
			if( (ingredients.getHGLnDOS().getCenterOfBin(n) >= _min) && (ingredients.getHGLnDOS().getCenterOfBin(n) <= _max) )
				if(HG_Energy_Rg2.getNumCountInBin(n) != 0)
			{
				file << std::setprecision(15) << HG_Energy_Rg2.getCenterOfBin(n) << "\t";
				file << std::setprecision(15) << ingredients.getHGLnDOS().getFirstMomentInBin(n) << "\t";
				file << std::setprecision(15) << HG_Energy_Rg2.getFirstMomentInBin(n) << "\t";
				file << std::setprecision(15) << HG_Energy_Rg2.getNumCountInBin(n) << "\n";

			}
		}

	file.close();

}

#endif
