
/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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

#ifndef LEMONADE_UPDATER_SETUP_COSOLVENT
#define LEMONADE_UPDATER_SETUP_COSOLVENT
/**
 * @file
 *
 * @class UpdaterAddCosolvent
 *
 * @brief Updater to add a cosolvent to the system
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding eigther an empty simulation box for system setup
 * or a prefilled ingredients where the cosolvent shall be added
 * @param NCos_ number of cosolvent molecules that are added to ingredients, each one consisting of 1 monomer
 * @param tagtype_ attribute tag of the monomers
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>

template<class IngredientsType>
class UpdaterAddCosolvent: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  UpdaterAddCosolvent(IngredientsType& ingredients_, uint32_t NCos_, int32_t tagtype_=1);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

//  //! cosolvent consists of unlinked monomers:
//  uint32_t NMonoPerChain = 1;

  //!number of linear chains in the box
  uint32_t NCos;

  //! lattice occupation density
  double density;

  //! bool for execution
  bool wasExecuted;

  //! attribute tag monomers
  int32_t tagtype;
//
//  //!getAttributeTag of odd monomers
//  int32_t type2;

};

/**
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param NCos_ number of chains to be added in the system instead of solvent
* @param NMonoPerChain_ number of monomers in one chain
*/
template < class IngredientsType >
UpdaterAddCosolvent<IngredientsType>::UpdaterAddCosolvent(IngredientsType& ingredients_, uint32_t NCos_, int32_t tagtype_):
BaseClass(ingredients_), NCos(NCos_), density(0), wasExecuted(false),
tagtype(tagtype_)
{}

/**
* @brief initialize function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddCosolvent<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterAddCosolvent" << std::endl;

  // get the target density from the sum of existing monomers and the new added chains
  density=(double)( ingredients.getMolecules().size() + NCos ) * 8  /(double)( ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ() );

  std::cout << "add "<<NCos<<" cosolvent monomers to the box"<<std::endl;

  execute();
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAddCosolvent<IngredientsType>::execute(){
  if(wasExecuted)
    return true;

  std::cout << "execute UpdaterAddCosolvent" << std::endl;

  //loop over chains and chain monomers and build it up
  for(uint32_t i=0;i<(NCos);i++){
    addSingleMonomer(tagtype);
    }


  ingredients.synchronize();
  double lattice_volume(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());
  if(density =! ( (double)(ingredients.getMolecules().size()*8) / lattice_volume ) ){
    std::cout << density << " " <<( (ingredients.getMolecules().size()*8) / lattice_volume)<<std::endl;
    throw std::runtime_error("UpdaterAddCosolvent: number of monomers in molecules does not match the calculated number of monomers!");
  }else{
    std::cout << "real lattice occupation density =" << (8*ingredients.getMolecules().size()) / lattice_volume<<std::endl;
    wasExecuted=true;
    linearizeSystem();
    return true;
  }
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddCosolvent<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_SETUP_COSOLVENT */
