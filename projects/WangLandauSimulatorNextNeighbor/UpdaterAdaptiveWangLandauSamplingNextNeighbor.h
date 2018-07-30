#ifndef ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR_H
#define ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR_H

#include <limits>
#include <iomanip>

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/UpdaterAbstractCreate.h>

#include "Histogram1D.h"


/**
 * @file
 *
 * @class UpdaterAdaptiveWangLandauSamplingNextNeighbor
 *
 * @brief Adaptive Umbrella Sampling updater along a reaction coordinate
 *
 * @details requires FeatureBilayerDistancePotential or other feature handling a 
 * possibly discrete potential along a reaction coordinate. this feature needs to
 * provide the functions: get/modifyPotential()  returning the potential object, 
 * getCurrentReactionCoordinate() returning the current value of the reaction
 * coordinate. also, the potential object needs to provide the functions 
 * addValue(coordinate, value) and getPMFMaxCoordinate(),getPMFMinCoordinate(),
 * getPMFValues() returning a std::vector<double> of the pmf values
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam <MoveType> name of the specialized move.
 */
template<class IngredientsType,class MoveType>
class UpdaterAdaptiveWangLandauSamplingNextNeighbor: public UpdaterAbstractCreate<IngredientsType>
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
	UpdaterAdaptiveWangLandauSamplingNextNeighbor(IngredientsType& ing
	,uint32_t steps
	,uint32_t stepsBeforeBiasCheck
	,double _initialModificationFactor=1.01
	,uint64_t maxAge=1000000000
	,double minWindow = -100.0
	,double maxWindow = +100.0
	,uint32_t stepsBeforeHistogramUpdate=1000
	,double convergenceThreshold=0.02,double fullFlatnessThreshold=0.1
	,double binCountThreshold=10000.0
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

		std::stringstream itr;
		itr << std::setw(2) << std::setfill('0') << iteration;

		dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_histogram.dat"));
		dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_totalhistogram.dat"));
		dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_HGLnDOS.dat"));
		dumpConvergenceProgress();
		//updateModificationFactor();
	};
	
	
	//! set the threshold for relative convergence of histogram (for update of bias)
	void setRelativeConvergenceThreshold(double threshold){convergenceThreshold=threshold;}
	
	//! return the threshold for relative convergence of histogram (for update of bias)
	double getRelativeConvergenceThreshold() const{return convergenceThreshold;}
	
	//! set the threshold for complete convergence of histogram (for final simulation convergence)
	void setFullConvergenceThreshold(double threshold){fullFlatnessThreshold=threshold;}
	
	//! return the threshold for complete convergence of histogram (for final simulation convergence)
	double getFullConvergenceThreshold() const {return fullFlatnessThreshold;}
	
	//! get the number of pmf updates performed up to now on the bias potential
	uint32_t getNPerformedUpdates()const{return nUpdatesPerformed;}
	
private:
	using BaseClass::ingredients;

	using BaseClass::randomBondvector;

	//! write the currently collected histogram to a file named histogram_mcs.dat
	void dumpHistogram(std::string filePrefix);
	
	//! write the currently collected histogram to a file named histogram_mcs.dat
	void dumpTotalHistogram(std::string filePrefix);

	//! writes the time development of the mean,variance, and square deviation from average to file names convergence.dat
	void dumpConvergenceProgress();
	
	//! update the used bias potential according to collected histogram
	
	void updateModificationFactor();

	//! write the currently logarithm of Density of States DOS i.e. file named HGLnDOS.dat
	void dumpHGLnDOS(std::string prefix="");

	//! check if the histogram has converged according to chosen criteria
	bool histogramConverged();
	
	//! reset the histogram to empty state
	void resetHistogram();
	
	//! check the mean square deviation from a flat histogram to check for full conversion of the simulation
	double histogramFlattness() const;
	
	//! check if mean and variance of histogram have converged. called as part of histogramConverged()
	bool distributionCheck();
	
	//! check if histogram flatness has converged, called as part of histogramConverged()
	bool flatnessCheck();
	
	//! set bin width of coordinate histogram and external potential
	bool setBinWidth();
	
	
	//! System information
	//IngredientsType& ingredients;
	
	//! Specialized move to be used
	MoveType move;
	
	//! Max number of mcs to be executed
	uint32_t nsteps;
	
	//! Number of steps between checks for update of bias
	const uint32_t nStepsBeforeBiasCheck;
	
	//! Number of checks before update of histogram
	const uint32_t nStepsBeforeHistogramUpdate;
	
	//! current counter to decide if update check is needed
	uint32_t counter_nStepsBeforeBiasCheck;
	
	//! current counter to decide if update check is needed
	uint32_t counter_nStepsBeforeHistogramUpdate;
	
	//! counts the number of bias updates performed so far
	uint32_t nUpdatesPerformed;
	
	//! histogram of reaction coordinate
	//Histogram1D histogram;
	
	//values for checking convergence
	double oldMean,oldOldMean;
	double oldVariance,oldOldVariance;
	double oldFlatnesValue,oldOldFlatnessValue;
	
	//thresholds for convergence
	double convergenceThreshold;
	double fullFlatnessThreshold; //for flatness
	double maxBinsThreshold;
	
	//time history of convergence criteria
	std::vector<double> varianceSeries;
	std::vector<double> meanSeries;
	std::vector<double> flatnessSeries;
	
	//! set to true if full convergence threshold has been reached
	bool simulationConverged;
	
	//! flags for setting or unsetting extra output of histogram
	bool writeHistogramProgress;
	
	//! flags for setting or unsetting extra output of bias potential
	bool writePMFProgress;
	
	//! maximum age of the system
	uint64_t maxSystemAge;
	
	double minWindow;
	double maxWindow;
	std::string prefixWindow;
	int iteration;

	double initialModificationFactor;
};


////////////////////////////////////////////////////////////////////////////////
// member implementations
///////////////////////////////////////////////////////////////////////////////


// constructor
template<class IngredientsType,class MoveType>
UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::
UpdaterAdaptiveWangLandauSamplingNextNeighbor(IngredientsType& ing,
				uint32_t steps, 
				uint32_t stepsBeforeBiasCheck,
				double _initialModificationFactor,
				uint64_t maxAge,
				double _minWindow,
				double _maxWindow,
				uint32_t stepsBeforeHistogramUpdate,
				double threshold,
				double fullThreshold,
				double binCountThreshold
				)
://ingredients(ing),
nsteps(steps)
,counter_nStepsBeforeBiasCheck(0)
,counter_nStepsBeforeHistogramUpdate(0)
,nStepsBeforeBiasCheck(stepsBeforeBiasCheck)
,nStepsBeforeHistogramUpdate(stepsBeforeHistogramUpdate)
,writeHistogramProgress(true)
,writePMFProgress(true)
,nUpdatesPerformed(0)
,oldMean(0.0)
,oldOldMean(0.0)
,oldVariance(0.0)
,oldOldVariance(0.0)
,oldFlatnesValue(10.0)
,oldOldFlatnessValue(10.0)
,convergenceThreshold(threshold)
,fullFlatnessThreshold(fullThreshold)
,maxBinsThreshold(binCountThreshold)
,initialModificationFactor(_initialModificationFactor)
,maxSystemAge(maxAge)
,minWindow(_minWindow)
,maxWindow(_maxWindow)
,prefixWindow("")
,iteration(0)
,simulationConverged(false)
,BaseClass(ing)
{
	varianceSeries.resize(0);
	meanSeries.resize(0);
	flatnessSeries.resize(0);
	ingredients.setModificationFactor(initialModificationFactor);
}
 
 
//execute simulation
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::execute()
{
	//some output for simulation performance and state
	time_t startTimer = std::time(NULL); //in seconds
	uint64_t systemAge=ingredients.getMolecules().getAge();
	std::cout<<"mcs "<<systemAge 
		 << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) <<std::endl;
	
	//if max age has reached, stop
	if(systemAge>=maxSystemAge) {
		std::cout<<"max age has been reached. maxAge "<<maxSystemAge<<" system age "<<systemAge<<std::endl;
		return false;
	}
	
	if(simulationConverged) return false;
	
	//simulation loop
	for(int n=0;n<nsteps;n++){
		
		//update bias potential if the histogram has converged well enough
		if(counter_nStepsBeforeBiasCheck==nStepsBeforeBiasCheck){
			counter_nStepsBeforeBiasCheck=0;
			if(writeHistogramProgress==true)
				{

					 dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_histogram.dat"));
					 dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_HGLnDOS.dat"));
				}
			if(histogramConverged()==true){
				
				/*double flat=histogramFlattness();
				if(histogramFlattness()<0.3) convergenceThreshold*=0.5;
				*/
				//write some file output
				dumpConvergenceProgress();
				if(writeHistogramProgress==true){
					std::stringstream filename;
					filename<<ingredients.getName() << prefixWindow  << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_histogram_mcs"
						<<ingredients.getMolecules().getAge()
						<<".dat";
					dumpHistogram(filename.str());

					std::stringstream filenameTotal;
					filenameTotal<<ingredients.getName() << prefixWindow << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_totalhistogram_mcs"<<ingredients.getMolecules().getAge()<<".dat";
					dumpTotalHistogram(filenameTotal.str());
				} 
				
				
				//write some more output
				if(writePMFProgress==true)
				{
					std::stringstream filenametmp;
					filenametmp<<ingredients.getName() << prefixWindow << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_HGLnDOS" << "_mcs"<<ingredients.getMolecules().getAge() << ".dat";

					dumpHGLnDOS(filenametmp.str());

					//dumpHGLnDOS(std::string(ingredients.getName()));
				}
				
				updateModificationFactor();
				iteration++;
				std::cout << "Start Iteration" << iteration << std::endl;

				//if the simulation has converged, end the simulation
				if(simulationConverged==true) return false;
				
				//set histogram to 0
				resetHistogram();
				
			}
			
		}
		
		
		//do one monte carlo sweep
		for(int m=0;m<ingredients.getMolecules().size();m++)
		{
			move.init(ingredients);
			
			if(move.check(ingredients)==true)
			{
				move.apply(ingredients);
			}
			else 
			{ 
				ingredients.rejectMove(ingredients);
			}
		}
		
		if( (ingredients.isEnergyInWindow() == false ) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients)+ingredients.getHGLnDOS().getBinwidth() < maxWindow) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients)-ingredients.getHGLnDOS().getBinwidth() > minWindow) )
		{
			std::cout << "RW is in energy window: [" << minWindow << " ; " << maxWindow << "] with " <<  ingredients.getInternalEnergyCurrentConfiguration(ingredients) << std::endl;
			// RW is in energy window
			ingredients.setWindowState(true, minWindow, maxWindow);

			std::stringstream ssprefixWindow;
			ssprefixWindow << "_minWin" << minWindow << "_maxWin" <<  maxWindow;

			prefixWindow = ssprefixWindow.str();

			// delete all histograms

			ingredients.modifyVisitsEnergyStates().reset(ingredients.getVisitsEnergyStates().getMinCoordinate(),ingredients.getVisitsEnergyStates().getMaxCoordinate(),ingredients.getVisitsEnergyStates().getNBins());
			ingredients.modifyTotalVisitsEnergyStates().reset(ingredients.getTotalVisitsEnergyStates().getMinCoordinate(),ingredients.getTotalVisitsEnergyStates().getMaxCoordinate(),ingredients.getTotalVisitsEnergyStates().getNBins());
			ingredients.modifyHGLnDOS().reset(ingredients.getHGLnDOS().getMinCoordinate(),ingredients.getHGLnDOS().getMaxCoordinate(),ingredients.getHGLnDOS().getNBins());

		}




		//get update on histogram after every sweep. more often would 
		//probably not make much sense, because the chain diffuses slowly
		try{
			if(counter_nStepsBeforeHistogramUpdate==nStepsBeforeHistogramUpdate){
				counter_nStepsBeforeHistogramUpdate=0;
				//histogram.addValue(ingredients.getInternalEnergy(ingredients));
			}
			
		}
		catch(std::runtime_error& e){
			std::stringstream errormessage;
			errormessage<<"UpdaterAdaptiveWangLandauSamplingNextNeighbor: error while altering histogram.\n"
			<<"original message was:\n"<<e.what()<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
		
		//increment counter
		counter_nStepsBeforeBiasCheck++;
		counter_nStepsBeforeHistogramUpdate++;
		
	}
	
	//update age of system and give some more output on simulation progress
	ingredients.modifyMolecules().setAge(ingredients.modifyMolecules().getAge()+nsteps);
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with " << (((1.0*nsteps)*ingredients.getMolecules().size())/(std::difftime(std::time(NULL), startTimer)) ) << " [attempted moves/s]" <<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) << " with " << nsteps << " MCS "<<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with energy " << ingredients.getInternalEnergyCurrentConfiguration(ingredients) <<std::endl;


	if(ingredients.getModificationFactor() < std::exp(std::pow(10,-8)) )
		return false;
	else
		return true;
}


// initialize only sets the histogram and bias potential to 0
template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::initialize()
{
	resetHistogram();
	dumpHistogram(std::string(ingredients.getName() + prefixWindow +  "_initial_histogram.dat"));
	dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow +  "_initial_totalhistogram.dat"));
	dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_initial"));
}


template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::resetHistogram()
{
	size_t nbins=ingredients.getVisitsEnergyStates().getNBins();
	double min = ingredients.getVisitsEnergyStates().getMinCoordinate();
	double max = ingredients.getVisitsEnergyStates().getMaxCoordinate();
	//histogram.reset(-10000,0.4,nbins);
	ingredients.modifyVisitsEnergyStates().reset(min,max,nbins);
}



//checks if mean, variance and deviation from average of the histogram have
//changed more than the percentage threshold. If all has converged, the function
//returns true
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::histogramConverged()
{
	

	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValues();//histogram.getVectorValues();

	double mean= 0.0;
	int entries = 0;
	//double min =currentHistogramState.at(n)

	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		if(currentHistogramState.at(n) != 0.0)
		{
			mean += currentHistogramState.at(n);
			entries++;
		}
	}

	mean = mean/(1.0*entries);
	
	if(entries <= 1)
		return false;

	bool isFlat = true;

	for(size_t n=0;n<currentHistogramState.size();n++)
		{
			if( currentHistogramState.at(n) != 0.0)
			//	if(std::abs((currentHistogramState.at(n)-mean)/mean) > 0.33)
					if(currentHistogramState.at(n)/mean < 0.8)
			{
					isFlat = false;
			}
		}

	return isFlat;

	/*
	bool distributionConverged=distributionCheck();
	bool flatnessConverged=flatnessCheck();
	
	bool converged= distributionConverged && flatnessConverged;
	
	std::cout<<"convergence checks: "<<converged<<std::endl;
		
	return converged;
	*/
	
}

//checks if the mean and variance of the reaction coordinate have changed more
//than percentage threshold. also records these two coordinates in the
//histogram convergence vectors for output in cleanup
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::distributionCheck()
{
	//calculate mean and variance
	double mean=0.0;
	double variance=0.0;
	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValuesNormalized();//histogram.getVectorValuesNormalized();
	
	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		mean+=currentHistogramState.at(n)*double(n);
		variance+=currentHistogramState.at(n)*double(n)*double(n);
	}
	variance-=mean*mean;
	
	//calculate the relative change with respect to the last check
	double meanValueChange=std::fabs((mean-oldMean)/mean);
	double varianceChange=std::fabs((variance-oldVariance)/variance);
	
	//calculate the relative change with respect to the second to last check
	double meanValueChangeSecond=std::fabs((mean-oldOldMean)/mean);
	double varianceChangeSecond=std::fabs((variance-oldOldVariance)/variance);
	
	if(meanValueChange == std::numeric_limits<double>::infinity()
		||meanValueChange == -std::numeric_limits<double>::infinity()
		||meanValueChange == std::numeric_limits<double>::signaling_NaN()
		||varianceChange == std::numeric_limits<double>::infinity()
		||varianceChange == -std::numeric_limits<double>::infinity()
		||varianceChange == std::numeric_limits<double>::signaling_NaN()
		||meanValueChangeSecond == std::numeric_limits<double>::infinity()
		||meanValueChangeSecond == -std::numeric_limits<double>::infinity()
		||meanValueChangeSecond == std::numeric_limits<double>::signaling_NaN()
		||varianceChangeSecond == std::numeric_limits<double>::infinity()
		||varianceChangeSecond == -std::numeric_limits<double>::infinity()
		||varianceChangeSecond == std::numeric_limits<double>::signaling_NaN())
	{
		std::cerr<<"WARNING: UpdaterAdaptiveWangLandauSamplingNextNeighbor finds NaN or inf in convergence check in mcs "
		<<ingredients.getMolecules().getAge()<<"\n";
		return false;
	}
	else if(meanValueChange>convergenceThreshold 
		||varianceChange>convergenceThreshold
		||meanValueChangeSecond>convergenceThreshold
		||varianceChangeSecond>convergenceThreshold){
		std::cout<<"distribution not converged!...\n"
		<<"mean change "<<meanValueChange<<" previous mean change "<<meanValueChangeSecond
		<<"\nvariance change "<<varianceChange<<" previous variance change "<<varianceChangeSecond
		<<"\nthreshold "<<convergenceThreshold<<std::endl;
	        oldOldMean=oldMean;
	        oldOldVariance=oldVariance;
		oldMean=mean;
		oldVariance=variance;
	        varianceSeries.push_back(variance);
	        meanSeries.push_back(mean);
		return false;
	}
	else{
		std::cout<<"distribution converged!...\n"
		<<"mean change "<<meanValueChange<<" previous mean change "<<meanValueChangeSecond
		<<"\nvariance change "<<varianceChange<<" previous variance change "<<varianceChangeSecond
		<<"\nthreshold "<<convergenceThreshold<<std::endl;
		
		oldOldMean=oldMean;
	        oldOldVariance=oldVariance;
		oldMean=mean;
		oldVariance=variance;
	        varianceSeries.push_back(variance);
	        meanSeries.push_back(mean);
		
		return true;
	}
	
}

//checks if the flatness of the histogram, measured by square deviation from 
//mean reached some limiting value (fullConvergenceThreshold), and if it has
//changed more than the percentage threshold (convergenceThreshold) in the
//last two checks. 
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::flatnessCheck()
{
	double currentFlatness=histogramFlattness();
	flatnessSeries.push_back(currentFlatness);
	
	
	double flatnessChange=std::sqrt(std::pow((currentFlatness-oldFlatnesValue)/currentFlatness,2));
	double flatnessChangeSecond=std::sqrt(std::pow((currentFlatness-oldOldFlatnessValue)/currentFlatness,2));
	
	std::cout<<"histogram flatness: "<<currentFlatness<<std::endl;
	std::cout<<"histogram flatness change: "<<flatnessChange<<" old flatness change: "<<flatnessChangeSecond<<std::endl;
	
	oldOldFlatnessValue=oldFlatnesValue;
	oldFlatnesValue=currentFlatness;
	
	if(currentFlatness==0)
	{
		std::cout<<"WARNING...found NaN or Inf in flatness check\n";
	}
		
		if(flatnessChange<=convergenceThreshold && flatnessChangeSecond<=convergenceThreshold && currentFlatness<=fullFlatnessThreshold){
			simulationConverged=true;
		}
// 		else{
// 			std::vector<double> currentHistogram=histogram.getVectorValues();
// 			bool binsHaveReachedMax=true;
// 			for(size_t n=0;n<currentHistogram.size();n++){
// 				binsHaveReachedMax &= (currentHistogram[n]>maxBinsThreshold);
// 			}
// 			
// 			if (binsHaveReachedMax) simulationConverged=true;
// 		}
		
	
	
	
	if(flatnessChange<=convergenceThreshold && flatnessChangeSecond<=convergenceThreshold){
		
		std::cout<<"flatness converged...\n";
		return true;
	}
	else return false;
}

//calculate the measure for the absolute histogram flatness
template<class IngredientsType, class MoveType>
double UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::histogramFlattness()const
{
	
	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValuesNormalized();//.histogram.getVectorValuesNormalized();
	
	//calculate mean
	double mean=1.0/double(currentHistogramState.size());
	double squareDeviation=0.0;
	
	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		squareDeviation+=(currentHistogramState.at(n)-mean)*(currentHistogramState.at(n)-mean)*double(currentHistogramState.size());
	}
	
	squareDeviation=std::sqrt(squareDeviation);
	
	
	return squareDeviation;
	
}

template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::updateModificationFactor()
{
	double newModificationFactor = std::pow(ingredients.getModificationFactor(), 0.5);
	//double newModificationFactor = std::exp(nsteps*1.0/ingredients.modifyMolecules().getAge());
	ingredients.setModificationFactor(newModificationFactor);

	std::cout << "new modification factor :" << std::setprecision(15) << newModificationFactor << std::endl;
}


//write current histogram to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpHistogram(std::string filename)
{

	std::ofstream file(filename.c_str());	

	std::vector<double> currentHistogram=ingredients.getVisitsEnergyStates().getVectorValues();//.histogram.getVectorValues();
	std::vector<double> bins=ingredients.getVisitsEnergyStates().getVectorBins();//histogram.getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getVisitsEnergyStates().getMinCoordinate() << " ; " << ingredients.getVisitsEnergyStates().getMaxCoordinate() << " ; " << ingredients.getVisitsEnergyStates().getNBins() << " ]" << std::endl;
	file << "# " << std::endl;


	for(size_t n=0;n<currentHistogram.size();n++){
		if(currentHistogram[n] != 0)
			file<< std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) << currentHistogram[n]<<"\n";
	}
	file.close();
	
}

//write current histogram to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpTotalHistogram(std::string filename)
{

	std::ofstream file(filename.c_str());

	std::vector<double> currentHistogram=ingredients.getTotalVisitsEnergyStates().getVectorValuesNormalized();//.histogram.getVectorValues();
	std::vector<double> bins=ingredients.getTotalVisitsEnergyStates().getVectorBins();//histogram.getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getTotalVisitsEnergyStates().getMinCoordinate() << " ; " << ingredients.getTotalVisitsEnergyStates().getMaxCoordinate() << " ; " << ingredients.getTotalVisitsEnergyStates().getNBins() << " ]" << std::endl;
	file << "# " << std::endl;

	for(size_t n=0;n<currentHistogram.size();n++){
		if(currentHistogram[n] != 0)
			file << std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) << currentHistogram[n]<<"\n";
	}
	file.close();

}


//write current bias potential to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpHGLnDOS(std::string prefix)
{
	std::stringstream filename;
	filename<< prefix;// << "_HGLnDOS" << "_mcs"<<ingredients.getMolecules().getAge()<<".dat";
	std::ofstream file(filename.str().c_str());

	std::vector<StatisticMoment> currentHGLnDOS=ingredients.getHGLnDOS().getVectorValues();
	std::vector<double> bins=ingredients.getHGLnDOS().getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getHGLnDOS().getMinCoordinate() << " ; " << ingredients.getHGLnDOS().getMaxCoordinate() << " ; " << ingredients.getHGLnDOS().getNBins() << " ]" << std::endl;
	file << "# " << std::endl;

	for(size_t n=0;n<currentHGLnDOS.size();n++){
		if(currentHGLnDOS[n].ReturnN() != 0)	
			file << std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) <<currentHGLnDOS[n].ReturnM1()<<"\n";
	}
	file.close();

}

//write convergence history to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpConvergenceProgress()
{
	std::stringstream filename;
	filename<<"convergence.dat";
	std::ofstream file;
	file.open(filename.str().c_str(),std::ofstream::out);
       
	
	for(size_t n=0;n<meanSeries.size();n++){
		file << std::setprecision(15) << meanSeries.at(n)<<" "<<varianceSeries.at(n)<<" "<<flatnessSeries.at(n)<<"\n";
	}
	
	file.close();
	
}


#endif
