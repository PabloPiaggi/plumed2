/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "PlumedMain.h"
#include "ActionSet.h"
#include "Action.h"
#include "DataGrabbingObject.h"
#include "Value.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/StoreDataVessel.h"

namespace PLMD {

template <class T>
class DataGrabbingObjectTyped : public DataGrabbingObject {
private:
/// A map containing the data we are grabbing
  std::map<std::string,T*> data;
public:
  explicit DataGrabbingObjectTyped(PlumedMain&plumed);
  void setData( const std::string& key, void* outval );
  void finishDataGrab();
};

DataGrabbingObject* DataGrabbingObject::create(unsigned n, PlumedMain& p){
  if(n==sizeof(double)){
    return new DataGrabbingObjectTyped<double>(p);
  } else  if(n==sizeof(float)){ 
    return new DataGrabbingObjectTyped<float>(p);
  } 
  std::string pp; Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

DataGrabbingObject::DataGrabbingObject(PlumedMain&p):
plumed(p)
{
}

DataGrabbingObject::~DataGrabbingObject(){}

bool DataGrabbingObject::activate() const {
  for(unsigned j=0;j<myactions.size();++j) myactions[j]->activate(); 
  if( myactions.size()>0 ) return true;
  return false;
} 

void DataGrabbingObject::grab_shape( ActionSet& a, const std::string& key, int* dims ){
  std::vector<std::string> words=Tools::getWords(key,"\t\n ,");
  vesselbase::ActionWithVessel* myv=a.selectWithLabel<vesselbase::ActionWithVessel*>( words[0] );
  if( myv ){
     plumed_assert( words.size()==1 );
     dims[0]=2; plumed_assert( myv ); dims[1]=myv->getFullNumberOfTasks(); dims[2]=myv->getNumberOfQuantities();
  } else {
     std::vector<Value*> arg; a.interpretArgumentList(words, arg, NULL );
     plumed_assert( arg.size()>0 ); dims[0]=1; dims[1]=arg.size();
  }
}

template <class T>
DataGrabbingObjectTyped<T>::DataGrabbingObjectTyped(PlumedMain&p):
DataGrabbingObject(p)
{
}

template <class T>
void DataGrabbingObjectTyped<T>::setData( const std::string& key, void* outval ){
  if( data.count(key) ) return;
  T* f=static_cast<T*>(outval);
  data.insert(std::pair<std::string,T*>(key,f)); 

  std::vector<std::string> words=Tools::getWords(key,"\t\n ,");
  vesselbase::ActionWithVessel* myv=plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>( words[0] );
  if( myv ){
     bool found=false; myv->buildDataStashes( false, 0.0, NULL );
     for(unsigned i=0;i<myactions.size();++i){
         if( myv->getLabel()==myactions[i]->getLabel() ) found=true;
     }   
     if(!found) myactions.push_back( myv ); 
  } else {
     std::vector<Value*> arg; plumed.getActionSet().interpretArgumentList(words, arg, NULL );
     for(unsigned i=0;i<arg.size();++i){
         bool found=false;
         for(unsigned j=0;j<myactions.size();++j){
             if( (arg[i]->getPntrToAction())->getLabel()==myactions[j]->getLabel() ) found=true;
         }
         if(!found) myactions.push_back( arg[i]->getPntrToAction() );
     }  
  }
}

template <class T>
void DataGrabbingObjectTyped<T>::finishDataGrab(){
  for(typename std::map<std::string,T*>::iterator it=data.begin();it!=data.end();++it){
      T* p = static_cast<T*>( (*it).second );
      std::vector<std::string> words=Tools::getWords((*it).first,"\t\n ,");
      vesselbase::ActionWithVessel* myv=plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>( words[0] );
      if( myv ){
          plumed_assert( words.size()==1 ); vesselbase::StoreDataVessel* vv=NULL;
          for(unsigned k=0;k<myv->getNumberOfVessels();++k){
              vv=dynamic_cast<vesselbase::StoreDataVessel*>( myv->getPntrToVessel(k) );
              if( vv ) break;
          }
          plumed_massert( vv, "error did not find data stash");

          unsigned k=0; std::vector<double> myvals( myv->getNumberOfQuantities() );
          for(unsigned i=0;i<myv->getFullNumberOfTasks();++i){
              vv->retrieveValue( i, true, myvals );
              for(unsigned j=0;j<myv->getNumberOfQuantities();++j){ p[k]=static_cast<T>(myvals[j]); k++; }
          } 
      } else {
          std::vector<Value*> arg; plumed.getActionSet().interpretArgumentList(words, arg, NULL );
          plumed_assert( arg.size()>0 ); for(unsigned i=0;i<arg.size();++i) p[i]=static_cast<T>(arg[i]->get());
      }
  }
}

}

