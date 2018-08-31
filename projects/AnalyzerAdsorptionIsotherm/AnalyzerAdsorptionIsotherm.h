/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Ron Dockhorn)
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

#ifndef AnalyzerAdsorptionIsotherm_H
#define AnalyzerAdsorptionIsotherm_H

#include <vector>
#include <string>
#include <utility>      // std::pair
#include <map>
#include <vector>

#include <LeMonADE/utility/Vector3D.h>

#include "StatisticMoment.h"



template<class IngredientsType>
class AnalyzerAdsorptionIsotherm:public AbstractAnalyzer
{
public:
	AnalyzerAdsorptionIsotherm(const IngredientsType& ing,  uint64_t startTime_, std::string dstDir_, int _numCoSolvent);


	virtual ~AnalyzerAdsorptionIsotherm(){

	};

	//typedef typename IngredientsType::molecules_type molecules_type;
	const typename IngredientsType::molecules_type& molecules ;

	const IngredientsType& getIngredients() const {return ingredients;}

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();


	double getInternalEnergyCurrentConfiguration() const;

	double getNumberCoSolventInNNShell() const;

private:

	const IngredientsType& ingredients;

	StatisticMoment Statistic_InternalEnergy;

	StatisticMoment Statistic_NumMonomers;

	uint64_t startTime;

	std::string filename;
	std::string dstdir;

	uint32_t numCoSolvent;
};




/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
AnalyzerAdsorptionIsotherm<IngredientsType>::AnalyzerAdsorptionIsotherm(const IngredientsType& ing,  uint64_t startTime_,  std::string dstDir_, int _numCoSolvent)
:ingredients(ing), molecules(ing.getMolecules()), startTime(startTime_),  dstdir(dstDir_), numCoSolvent(_numCoSolvent)
 {

	Statistic_InternalEnergy.clear();

	Statistic_NumMonomers.clear();

 }

template<class IngredientsType>
void AnalyzerAdsorptionIsotherm<IngredientsType>::initialize()
{

	//execute();

}

template<class IngredientsType>
bool AnalyzerAdsorptionIsotherm<IngredientsType>::execute()
{
	// time of the conformations
	//uint64_t timeInSim =  ingredients.getMolecules().getAge();
	//molecules <-> ingredients.getMolecules()

	if(ingredients.getMolecules().getAge() >= startTime)
	{
		std::cout << "AnalyzerAdsorptionIsotherm.execute() at MCS:" << ingredients.getMolecules().getAge() << std::endl;

		double interalEnergy = getInternalEnergyCurrentConfiguration();

		Statistic_InternalEnergy.AddValue(interalEnergy);

		double numCoSolventInNNShell = getNumberCoSolventInNNShell();

		Statistic_NumMonomers.AddValue(numCoSolventInNNShell);

	}


}

template<class IngredientsType>
double AnalyzerAdsorptionIsotherm<IngredientsType>::getInternalEnergyCurrentConfiguration() const
{
	std::vector<VectorInt3> contactSites;

	//prepare the contact shell, in which the contacts are counted

	contactSites.push_back(VectorInt3(2,0,0));
	contactSites.push_back(VectorInt3(2,1,0));
	contactSites.push_back(VectorInt3(2,0,1));
	contactSites.push_back(VectorInt3(2,1,1));
	contactSites.push_back(VectorInt3(0,2,0));
	contactSites.push_back(VectorInt3(1,2,0));
	contactSites.push_back(VectorInt3(0,2,1));
	contactSites.push_back(VectorInt3(1,2,1));
	contactSites.push_back(VectorInt3(0,0,2));
	contactSites.push_back(VectorInt3(1,0,2));
	contactSites.push_back(VectorInt3(0,1,2));
	contactSites.push_back(VectorInt3(1,1,2));

	contactSites.push_back(VectorInt3(-1,0,0));
	contactSites.push_back(VectorInt3(-1,1,0));
	contactSites.push_back(VectorInt3(-1,0,1));
	contactSites.push_back(VectorInt3(-1,1,1));
	contactSites.push_back(VectorInt3(0,-1,0));
	contactSites.push_back(VectorInt3(1,-1,0));
	contactSites.push_back(VectorInt3(0,-1,1));
	contactSites.push_back(VectorInt3(1,-1,1));
	contactSites.push_back(VectorInt3(0,0,-1));
	contactSites.push_back(VectorInt3(1,0,-1));
	contactSites.push_back(VectorInt3(0,1,-1));
	contactSites.push_back(VectorInt3(1,1,-1));

	VectorInt3 pos;
	double Energy=0.0;

	//loop through all polymer
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();
		pos=ingredients.getMolecules()[n];


		for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
		{
			int32_t latticeEntry=int32_t(ingredients.getLatticeEntry(pos+contactSites[contactNo]));

			if(latticeEntry != 0)
				Energy += ingredients.getNNInteraction(monoType,latticeEntry);
		}
	}

	// factor 0.5 due to the double sum in hamiltonian
	return 0.5*Energy;
}

template<class IngredientsType>
double AnalyzerAdsorptionIsotherm<IngredientsType>::getNumberCoSolventInNNShell() const
{
	std::vector<VectorInt3> contactSites;

	//prepare the contact shell, in which the contacts are counted

	contactSites.push_back(VectorInt3(2,0,0));
	contactSites.push_back(VectorInt3(2,1,0));
	contactSites.push_back(VectorInt3(2,0,1));
	contactSites.push_back(VectorInt3(2,1,1));
	contactSites.push_back(VectorInt3(0,2,0));
	contactSites.push_back(VectorInt3(1,2,0));
	contactSites.push_back(VectorInt3(0,2,1));
	contactSites.push_back(VectorInt3(1,2,1));
	contactSites.push_back(VectorInt3(0,0,2));
	contactSites.push_back(VectorInt3(1,0,2));
	contactSites.push_back(VectorInt3(0,1,2));
	contactSites.push_back(VectorInt3(1,1,2));

	contactSites.push_back(VectorInt3(-1,0,0));
	contactSites.push_back(VectorInt3(-1,1,0));
	contactSites.push_back(VectorInt3(-1,0,1));
	contactSites.push_back(VectorInt3(-1,1,1));
	contactSites.push_back(VectorInt3(0,-1,0));
	contactSites.push_back(VectorInt3(1,-1,0));
	contactSites.push_back(VectorInt3(0,-1,1));
	contactSites.push_back(VectorInt3(1,-1,1));
	contactSites.push_back(VectorInt3(0,0,-1));
	contactSites.push_back(VectorInt3(1,0,-1));
	contactSites.push_back(VectorInt3(0,1,-1));
	contactSites.push_back(VectorInt3(1,1,-1));

	VectorInt3 pos;
	double numCoSolventInShell=0.0;

	//loop through all polymer
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();

		// only for cosolvent
		if( monoType == 3)
		{
			pos=ingredients.getMolecules()[n];


			for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
			{
				int32_t latticeEntry=int32_t(ingredients.getLatticeEntry(pos+contactSites[contactNo]));

				// check if surrounding places are polymer
				if((latticeEntry == 1))
				{
					numCoSolventInShell++;
					break;
				}
			}
		}
	}

	return numCoSolventInShell;
}


struct PathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
void AnalyzerAdsorptionIsotherm<IngredientsType>::cleanup()
{
	std::cout << "File output" << std::endl;
	std::cout<<"Average Values: = " << Statistic_InternalEnergy.ReturnN()   << std::endl;

	std::vector < std::vector<double> > tmpResults;

		tmpResults.resize(14);

		for(int i = 0; i < 14; i++)
			tmpResults[i].resize(1);

		tmpResults[0][0]=ingredients.getNNInteraction(1, 3);
		tmpResults[1][0]=Statistic_InternalEnergy.ReturnM1();
		tmpResults[2][0]=Statistic_InternalEnergy.ReturnM2();
		tmpResults[3][0]=Statistic_InternalEnergy.ReturnM2()-Statistic_InternalEnergy.ReturnM1()*Statistic_InternalEnergy.ReturnM1();

		tmpResults[4][0]=Statistic_InternalEnergy.ReturnM1()/ingredients.getNNInteraction(1, 3);

		tmpResults[5][0]=(Statistic_InternalEnergy.ReturnM1()/ingredients.getNNInteraction(1, 3))/(8.0*(2.0*256+1));

		tmpResults[6][0]=(Statistic_InternalEnergy.ReturnM1()/ingredients.getNNInteraction(1, 3))/(4.0*numCoSolvent);

		tmpResults[7][0]=numCoSolvent;

		tmpResults[8][0]=Statistic_NumMonomers.ReturnM1();
		tmpResults[9][0]=Statistic_NumMonomers.ReturnM2();
		tmpResults[10][0]=Statistic_NumMonomers.ReturnM2()-Statistic_NumMonomers.ReturnM1()*Statistic_NumMonomers.ReturnM1();

		tmpResults[11][0]=Statistic_NumMonomers.ReturnM1()/numCoSolvent;


		tmpResults[12][0]=((8.0*(ingredients.getMolecules().size()))/(1.0*ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()));

		tmpResults[13][0]=8.0*numCoSolvent/(1.0*ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());


		std::stringstream comment;
		comment <<"File produced by analyzer AnalyzerAdsorptionIsotherm\n"
				<<"Analyze AdsorptionIsotherm\n"
				<<"NrDensity c="<<((8.0*(ingredients.getMolecules().size()))/(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ())) << "\n"
				<<"Number CoSolvency=" << numCoSolvent << "\n"
				<<"\n"
				<<"epsilon\t<U>\t<U²>\t<cV>\t<nContacts>\t<nContacts/nContacts,max>\t<nContacts/4nCoSolvent>\t<nCoSolvent>\t<nCoS In NNShell>\t<(nCoS In NNShell)²>\tvar(nCoS In NNShell)\t<(nCoS In NNShell)/nCoSolvent>\t<cAll>\t<cCoSolvent>\n";




	// find the filename without path and extensions
	std::string filenameGeneral= std::string( std::find_if( ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator() ).base(), ingredients.getName().end() );

	std::string::size_type const p(filenameGeneral.find_last_of('.'));
	filenameGeneral = filenameGeneral.substr(0, p);

	std::string filenameRg2_Ree_b2 = filenameGeneral + "_AdsorptionIsotherm.dat";

	std::cout  << " Write output to: " << dstdir  <<"/" << filenameRg2_Ree_b2 << std::endl;

	ResultFormattingTools::writeResultFile(dstdir+"/"+filenameRg2_Ree_b2, this->ingredients, tmpResults, comment.str());


}

#endif /*AnalyzerAdsorptionIsotherm_H*/
