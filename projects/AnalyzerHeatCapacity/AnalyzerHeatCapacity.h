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

#ifndef AnalyzerHeatCapacity_H
#define AnalyzerHeatCapacity_H

#include <vector>
#include <string>
#include <utility>      // std::pair
#include <map>
#include <vector>

#include <LeMonADE/utility/Vector3D.h>

#include "StatisticMoment.h"



template<class IngredientsType>
class AnalyzerHeatCapacity:public AbstractAnalyzer
{
public:
	AnalyzerHeatCapacity(const IngredientsType& ing,  uint64_t startTime_, std::string dstDir_);


	virtual ~AnalyzerHeatCapacity(){

	};

	//typedef typename IngredientsType::molecules_type molecules_type;
	const typename IngredientsType::molecules_type& molecules ;

	const IngredientsType& getIngredients() const {return ingredients;}

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();

private:

	const IngredientsType& ingredients;

	StatisticMoment Statistic_InternalEnergy;

	StatisticMoment Statistic_Rg2;
	StatisticMoment Statistic_Rg2TimesE;

	uint64_t startTime;

	std::string filename;
	std::string dstdir;
};




/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
AnalyzerHeatCapacity<IngredientsType>::AnalyzerHeatCapacity(const IngredientsType& ing,  uint64_t startTime_,  std::string dstDir_)
:ingredients(ing), molecules(ing.getMolecules()), startTime(startTime_),  dstdir(dstDir_)
 {

	Statistic_InternalEnergy.clear();
	Statistic_Rg2.clear();
	Statistic_Rg2TimesE.clear();

 }

template<class IngredientsType>
void AnalyzerHeatCapacity<IngredientsType>::initialize()
{

	//execute();

}

template<class IngredientsType>
bool AnalyzerHeatCapacity<IngredientsType>::execute()
{
	// time of the conformations
	//uint64_t timeInSim =  ingredients.getMolecules().getAge();
	//molecules <-> ingredients.getMolecules()

	if(ingredients.getMolecules().getAge() >= startTime)
	{
		std::cout << "SimpleAnalyzer_Rg2.execute() at MCS:" << ingredients.getMolecules().getAge() << std::endl;



		int monomerCounter = 0;

		double Rg2 = 0.0;
		double Rg2_x = 0.0;
		double Rg2_y = 0.0;
		double Rg2_z = 0.0;

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l= k; l < ingredients.getMolecules().size(); l++)
			//	 if((ingredients.getMolecules()[k].getAttributeTag()==1) && (ingredients.getMolecules()[l].getAttributeTag()==1))
			{
				Rg2_x += (ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX())*(ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX());
				Rg2_y += (ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY())*(ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY());
				Rg2_z += (ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ())*(ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ());


			}

			//if((ingredients.getMolecules()[k].getAttributeTag()==1))
				monomerCounter++;
		}


		Rg2_x /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_y /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_z /= 1.0*(monomerCounter*monomerCounter);//ingredients.getMolecules().size()*ingredients.getMolecules().size());

		Rg2 = Rg2_x+Rg2_y+Rg2_z;


		Statistic_Rg2.AddValue(Rg2);

		double interalEnergy = ingredients.getInternalEnergyCurrentConfiguration(ingredients);

		Statistic_InternalEnergy.AddValue(interalEnergy);

		Statistic_Rg2TimesE.AddValue(Rg2*interalEnergy);

		/*
		//calculate the Bondlength
		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			uint32_t numLinks = ingredients.getMolecules().getNumLinks(k);

			for(uint32_t nl = 0; nl < numLinks; nl++)
			{
				uint32_t idxNeighbor = ingredients.getMolecules().getNeighborIdx(k, nl);

				if(idxNeighbor > k)
			//		if((ingredients.getMolecules()[k].getAttributeTag()==1) && (ingredients.getMolecules()[idxNeighbor].getAttributeTag()==1))
				{
					VectorDouble3 bondlength(ingredients.getMolecules()[idxNeighbor] - ingredients.getMolecules()[k]);
					Statistic_BondLengthSquared.AddValue(bondlength*bondlength);

					Statistic_BondLengthSquared_x.AddValue(bondlength.getX()*bondlength.getX());
					Statistic_BondLengthSquared_y.AddValue(bondlength.getY()*bondlength.getY());
					Statistic_BondLengthSquared_z.AddValue(bondlength.getZ()*bondlength.getZ());
				}
			}


		}
		*/
	}


}

struct PathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
void AnalyzerHeatCapacity<IngredientsType>::cleanup()
{
	std::cout << "File output" << std::endl;
	std::cout<<"Average Values: = " << Statistic_InternalEnergy.ReturnN()   << std::endl;

	std::vector < std::vector<double> > tmpResultsRg2_Ree_b2;

		tmpResultsRg2_Ree_b2.resize(9);

		for(int i = 0; i < 9; i++)
			tmpResultsRg2_Ree_b2[i].resize(1);

		tmpResultsRg2_Ree_b2[0][0]=ingredients.getMolecules().size();


		tmpResultsRg2_Ree_b2[1][0]=Statistic_InternalEnergy.ReturnM2()-Statistic_InternalEnergy.ReturnM1()*Statistic_InternalEnergy.ReturnM1();
		tmpResultsRg2_Ree_b2[2][0]=Statistic_InternalEnergy.ReturnM1();
		tmpResultsRg2_Ree_b2[3][0]=Statistic_InternalEnergy.ReturnM2();

		tmpResultsRg2_Ree_b2[4][0]=Statistic_Rg2.ReturnM1();
		tmpResultsRg2_Ree_b2[5][0]=Statistic_Rg2.ReturnM2();

		tmpResultsRg2_Ree_b2[6][0]=Statistic_Rg2TimesE.ReturnM1()-Statistic_Rg2.ReturnM1()*Statistic_InternalEnergy.ReturnM1();
		tmpResultsRg2_Ree_b2[7][0]=Statistic_Rg2TimesE.ReturnM1();
		tmpResultsRg2_Ree_b2[8][0]=Statistic_Rg2TimesE.ReturnM2();

		std::stringstream comment;
		comment <<"File produced by analyzer AnalyzerHeatCapacity\n"
				<<"Analyze HeatCapacity\n"
				<<"NrDensity c="<<((8.0*(ingredients.getMolecules().size()))/(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ())) << "\n"
				<<"\n"
				<<"N\t<cV>\t<U>\t<U²>\t<Rg²>\t<(Rg²)²>\t<dRg²/dT>*T\t<Rg²*U>\t<(Rg²*U)²>" << "\n";




	// find the filename without path and extensions
	std::string filenameGeneral= std::string( std::find_if( ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator() ).base(), ingredients.getName().end() );

	std::string::size_type const p(filenameGeneral.find_last_of('.'));
	filenameGeneral = filenameGeneral.substr(0, p);

	std::string filenameRg2_Ree_b2 = filenameGeneral + "_HeatCapacity_Rg2Fluctuation.dat";

	std::cout  << " Write output to: " << dstdir  <<"/" << filenameRg2_Ree_b2 << std::endl;

	ResultFormattingTools::writeResultFile(dstdir+"/"+filenameRg2_Ree_b2, this->ingredients, tmpResultsRg2_Ree_b2, comment.str());


}

#endif /*AnalyzerHeatCapacity_H*/
