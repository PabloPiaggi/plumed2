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
#include "tools/Tools.h"
#include "tools/OpenMP.h"
#include <cstring>
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include <set>
#include "config/Config.h"
#include <cstdlib>
#include <stdint.h>
#include "ActionRegister.h"
#include "GREX.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "ActionSet.h"
#include "tools/Log.h"
#include "tools/DLLoader.h"
#include "tools/Communicator.h"
#include "CLToolMain.h"
#include "tools/Stopwatch.h"
#include "tools/Citations.h"
#include "ExchangePatterns.h"
#include "tools/IFile.h"
#include "multicolvar/MultiColvarBase.h"
#include "vesselbase/StoreDataVessel.h"
#include "vesselbase/ActionWithVessel.h"

using namespace std;

enum { SETBOX, SETPOSITIONS, SETMASSES, SETCHARGES, SETPOSITIONSX, SETPOSITIONSY, SETPOSITIONSZ, SETVIRIAL, SETENERGY, SETFORCES, SETFORCESX, SETFORCESY, SETFORCESZ, CALC, PREPAREDEPENDENCIES, SHAREDATA, PREPARECALC, PERFORMCALC, SETSTEP, SETSTEPLONG, SETATOMSNLOCAL, SETATOMSGATINDEX, SETATOMSFGATINDEX, SETATOMSCONTIGUOUS, CREATEFULLLIST, GETFULLLIST, CLEARFULLLIST, READ, ACTION, CLEAR, GETAPIVERSION, GETINTEGERPRECISION, INIT, SETREALPRECISION, SETMDLENGTHUNITS, SETMDENERGYUNITS, SETMDTIMEUNITS, SETNATURALUNITS, SETNOVIRIAL, SETPLUMEDDAT, SETMPICOMM, SETMPIFCOMM, SETMPIMULTISIMCOMM, SETNATOMS, SETTIMESTEP, SETMDENGINE, SETLOG, SETLOGFILE, SETSTOPFLAG, GETEXCHANGESFLAG, SETEXCHANGESSEED, SETNUMBEROFREPLICAS, GETEXCHANGESLIST, RUNFINALJOBS, ISENERGYNEEDED, GETBIAS, SETKBT, SETRESTART };

namespace PLMD{

PlumedMain::PlumedMain():
  comm(*new Communicator),
  multi_sim_comm(*new Communicator),
  dlloader(*new DLLoader),
  cltool(NULL),
  stopwatch(*new Stopwatch),
  grex(NULL),
  initialized(false),
  log(*new Log),
  citations(*new Citations),
  step(0),
  active(false),
  atoms(*new Atoms(*this)),
  actionSet(*new ActionSet(*this)),
  bias(0.0),
  work(0.0),
  exchangePatterns(*new(ExchangePatterns)),
  exchangeStep(false),
  restart(false),
  stopFlag(NULL),
  stopNow(false),
  novirial(false),
  detailedTimers(false)
{
  log.link(comm);
  log.setLinePrefix("PLUMED: ");
  stopwatch.start();
  stopwatch.pause();
  word_map["setBox"]=SETBOX;
  word_map["setPositions"]=SETPOSITIONS;
  word_map["setMasses"]=SETMASSES;
  word_map["setCharges"]=SETCHARGES;
  word_map["setPositionsX"]=SETPOSITIONSX;
  word_map["setPositionsY"]=SETPOSITIONSY;
  word_map["setPositionsZ"]=SETPOSITIONSZ;
  word_map["setVirial"]=SETVIRIAL;
  word_map["setEnergy"]=SETENERGY;
  word_map["setForces"]=SETFORCES;
  word_map["setForcesX"]=SETFORCESX;
  word_map["setForcesY"]=SETFORCESY;
  word_map["setForcesZ"]=SETFORCESZ;
  word_map["calc"]=CALC;
  word_map["prepareDependencies"]=PREPAREDEPENDENCIES;
  word_map["shareData"]=SHAREDATA;
  word_map["prepareCalc"]=PREPARECALC;
  word_map["performCalc"]=PERFORMCALC;
  word_map["setStep"]=SETSTEP;
  word_map["setStepLong"]=SETSTEPLONG;
  word_map["setAtomsNlocal"]=SETATOMSNLOCAL;
  word_map["setAtomsGatindex"]=SETATOMSGATINDEX;
  word_map["setAtomsFGatindex"]=SETATOMSFGATINDEX;
  word_map["setAtomsContiguous"]=SETATOMSCONTIGUOUS;
  word_map["createFullList"]=CREATEFULLLIST;
  word_map["getFullList"]=GETFULLLIST;
  word_map["clearFullList"]=CLEARFULLLIST;
  word_map["read"]=READ;
  word_map["action"]=ACTION;
  word_map["clear"]=CLEAR;
  word_map["getApiVersion"]=GETAPIVERSION;
  word_map["getIntegerPrecision"]=GETINTEGERPRECISION;
  word_map["init"]=INIT;
  word_map["setRealPrecision"]=SETREALPRECISION;
  word_map["setMDLengthUnits"]=SETMDLENGTHUNITS;
  word_map["setMDEnergyUnits"]=SETMDENERGYUNITS;
  word_map["setMDTimeUnits"]=SETMDTIMEUNITS;
  word_map["setNaturalUnits"]=SETNATURALUNITS;
  word_map["setNoVirial"]=SETNOVIRIAL;
  word_map["setPlumedDat"]=SETPLUMEDDAT;
  word_map["setMPIComm"]=SETMPICOMM;
  word_map["setMPIFComm"]=SETMPIFCOMM;
  word_map["setMPImultiSimComm"]=SETMPIMULTISIMCOMM;
  word_map["setNatoms"]=SETNATOMS;
  word_map["setTimestep"]=SETTIMESTEP;
  word_map["setMDEngine"]=SETMDENGINE;
  word_map["setLog"]=SETLOG;
  word_map["setLogFile"]=SETLOGFILE;
  word_map["setStopFlag"]=SETSTOPFLAG;
  word_map["getExchangesFlag"]=GETEXCHANGESFLAG;
  word_map["setExchangesSeed"]=SETEXCHANGESSEED;
  word_map["setNumberOfReplicas"]=SETNUMBEROFREPLICAS;
  word_map["getExchangesList"]=GETEXCHANGESLIST;
  word_map["runFinalJobs"]=RUNFINALJOBS;
  word_map["isEnergyNeeded"]=ISENERGYNEEDED;
  word_map["getBias"]=GETBIAS;
  word_map["setKbT"]=SETKBT;
  word_map["setRestart"]=SETRESTART;
}

PlumedMain::~PlumedMain(){
  stopwatch.start();
  stopwatch.stop();
  if(initialized) log<<stopwatch;
  delete &exchangePatterns;
  delete &actionSet;
  delete &citations;
  delete &atoms;
  delete &log;
  if(grex)  delete grex;
  delete &stopwatch;
  if(cltool) delete cltool;
  delete &dlloader;
  delete &comm;
  delete &multi_sim_comm;
}

/////////////////////////////////////////////////////////////
//  MAIN INTERPRETER

#define CHECK_INIT(ini,word) plumed_massert(ini,"cmd(\"" + word +"\") should be only used after plumed initialization")
#define CHECK_NOTINIT(ini,word) plumed_massert(!(ini),"cmd(\"" + word +"\") should be only used before plumed initialization")
#define CHECK_NULL(val,word) plumed_massert(val,"NULL pointer received in cmd(\"" + word + "\")");

void PlumedMain::cmd(const std::string & word,void*val){

  stopwatch.start();

  std::vector<std::string> words=Tools::getWords(word);
  unsigned nw=words.size();
  if(nw==1) {
    switch(word_map[word]) {
      double d;
      case SETBOX:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setBox(val);
        break;
      case SETPOSITIONS:
        CHECK_INIT(initialized,word);
        atoms.setPositions(val);
        break;
      case SETMASSES:
        CHECK_INIT(initialized,word);
        atoms.setMasses(val);
        break;
      case SETCHARGES:
        CHECK_INIT(initialized,word);
        atoms.setCharges(val);
        break;
      case SETPOSITIONSX:
        CHECK_INIT(initialized,word);
        atoms.setPositions(val,0);
        break;
      case SETPOSITIONSY:
        CHECK_INIT(initialized,word);
        atoms.setPositions(val,1);
        break;
      case SETPOSITIONSZ:
        CHECK_INIT(initialized,word);
        atoms.setPositions(val,2);
        break;
      case SETVIRIAL:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setVirial(val);
        break;
      case SETENERGY:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setEnergy(val);
        break;
      case SETFORCES:
        CHECK_INIT(initialized,word);
        atoms.setForces(val);
        break;
      case SETFORCESX:
        CHECK_INIT(initialized,word);
        atoms.setForces(val,0);
        break;
      case SETFORCESY:
        CHECK_INIT(initialized,word);
        atoms.setForces(val,1);
        break;
      case SETFORCESZ:
        CHECK_INIT(initialized,word);
        atoms.setForces(val,2);
        break;
      case CALC:
        CHECK_INIT(initialized,word);
        calc();
        break;
      case PREPAREDEPENDENCIES:
        CHECK_INIT(initialized,word);
        prepareDependencies();
        break;
      case SHAREDATA:
        CHECK_INIT(initialized,word);
        shareData();
        break;
      case PREPARECALC:
        CHECK_INIT(initialized,word);
        prepareCalc();
        break;
      case PERFORMCALC:
        CHECK_INIT(initialized,word);
        performCalc();
        break;
      case SETSTEP:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        step=(*static_cast<int*>(val));
        atoms.startStep();
        break;
      case SETSTEPLONG:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        step=(*static_cast<long int*>(val));
        atoms.startStep();
        break;
      // words used less frequently:
      case SETATOMSNLOCAL:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setAtomsNlocal(*static_cast<int*>(val));
        break;
      case SETATOMSGATINDEX:
        CHECK_INIT(initialized,word);
        atoms.setAtomsGatindex(static_cast<int*>(val),false);
        break;
      case SETATOMSFGATINDEX:
        CHECK_INIT(initialized,word);
        atoms.setAtomsGatindex(static_cast<int*>(val),true);
        break;
      case SETATOMSCONTIGUOUS:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setAtomsContiguous(*static_cast<int*>(val));
        break;
      case CREATEFULLLIST:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.createFullList(static_cast<int*>(val));
        break;
      case GETFULLLIST:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.getFullList(static_cast<int**>(val));
        break;
      case CLEARFULLLIST:
        CHECK_INIT(initialized,word);
        atoms.clearFullList();
        break;
      case READ:
        CHECK_INIT(initialized,word);
        if(val)readInputFile(static_cast<char*>(val));
        else   readInputFile("plumed.dat");
        break;
      case ACTION:
        CHECK_INIT(initialized,word);
        if(val)readAction(static_cast<char*>(val));
        break;
      case CLEAR:
        CHECK_INIT(initialized,word);
        actionSet.clearDelete();
        break;
      case GETAPIVERSION:
        CHECK_NULL(val,word);
        *(static_cast<int*>(val))=3;
        break;
      case GETINTEGERPRECISION:
        CHECK_NULL(val,word);
        *(static_cast<uint8_t*>(val))=(sizeof(int));
        break;
      // commands which can be used only before initialization:
      case INIT:
        CHECK_NOTINIT(initialized,word);
        init();
        break;
      case SETREALPRECISION:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setRealPrecision(*static_cast<int*>(val));
        break;
      case SETMDLENGTHUNITS:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.MD2double(val,d);
        atoms.setMDLengthUnits(d);
        break;
      case SETMDENERGYUNITS:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.MD2double(val,d);
        atoms.setMDEnergyUnits(d);
        break;
      case SETMDTIMEUNITS:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.MD2double(val,d);
        atoms.setMDTimeUnits(d);
        break;
      case SETNATURALUNITS:
      // set the boltzman constant for MD in natural units (kb=1)
      // only needed in LJ codes if the MD is passing temperatures to plumed (so, not yet...)
      // use as cmd("setNaturalUnits")
        CHECK_NOTINIT(initialized,word);
        atoms.setMDNaturalUnits(true);
        break;
      case SETNOVIRIAL:
        CHECK_NOTINIT(initialized,word);
        novirial=true;
        break;
      case SETPLUMEDDAT:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        plumedDat=static_cast<char*>(val);
        break;
      case SETMPICOMM:
        CHECK_NOTINIT(initialized,word);
        comm.Set_comm(val);
        atoms.setDomainDecomposition(comm);
        break;
      case SETMPIFCOMM:
        CHECK_NOTINIT(initialized,word);
        comm.Set_fcomm(val);
        atoms.setDomainDecomposition(comm);
        break;
      case SETMPIMULTISIMCOMM:
        CHECK_NOTINIT(initialized,word);
        multi_sim_comm.Set_comm(val);
        break;
      case SETNATOMS:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setNatoms(*static_cast<int*>(val));
        break;
      case SETTIMESTEP:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setTimeStep(val);
        break;
      case SETKBT: /* ADDED WITH API==2 */
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        atoms.setKbT(val);
        break;
      case SETRESTART: /* ADDED WITH API==3 */
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        if(*static_cast<int*>(val)!=0) restart=true;
        break;
      case SETMDENGINE:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        MDEngine=static_cast<char*>(val);
        break;
      case SETLOG:
        CHECK_NOTINIT(initialized,word);
        log.link(static_cast<FILE*>(val));
        break;
      case SETLOGFILE:
        CHECK_NOTINIT(initialized,word);
        CHECK_NULL(val,word);
        log.open(static_cast<char*>(val));
        break;
      // other commands that should be used after initialization:
      case SETSTOPFLAG:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        stopFlag=static_cast<int*>(val);
        break;
      case GETEXCHANGESFLAG:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        exchangePatterns.getFlag((*static_cast<int*>(val)));
        break;
      case SETEXCHANGESSEED:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        exchangePatterns.setSeed((*static_cast<int*>(val)));
        break;
      case SETNUMBEROFREPLICAS:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        exchangePatterns.setNofR((*static_cast<int*>(val)));
        break;
      case GETEXCHANGESLIST:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        exchangePatterns.getList((static_cast<int*>(val)));
        break;
      case RUNFINALJOBS:
        CHECK_INIT(initialized,word);
        runJobsAtEndOfCalculation();
        break;
      case ISENERGYNEEDED:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        if(atoms.isEnergyNeeded()) *(static_cast<int*>(val))=1;
        else                       *(static_cast<int*>(val))=0;
        break;
      case GETBIAS:
        CHECK_INIT(initialized,word);
        CHECK_NULL(val,word);
        d=getBias()/(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy());
        atoms.double2MD(d,val);
        break;
      default:
        plumed_merror("cannot interpret cmd(\"" + word + "\"). check plumed developers manual to see the available commands.");
        break;
    }
  } else if(nw==2 && words[0]=="checkAction"){
    int check=0;
    if(actionRegister().check(words[1])) check=1;
    *(static_cast<int*>(val))=check;
  } else if(nw>1 && words[0]=="GREX"){
    if(!grex) grex=new GREX(*this);
    plumed_massert(grex,"error allocating grex");
    std::string kk=words[1];
    for(unsigned i=2;i<words.size();i++) kk+=" "+words[i];
    grex->cmd(kk.c_str(),val);
  } else if(nw>1 && words[0]=="CLTool"){
    CHECK_NOTINIT(initialized,word);
    if(!cltool) cltool=new CLToolMain;
    std::string kk=words[1];
    for(unsigned i=2;i<words.size();i++) kk+=" "+words[i];
    cltool->cmd(kk.c_str(),val);
  } else{
    plumed_merror("cannot interpret cmd(\"" + word + "\"). check plumed developers manual to see the available commands.");
  };
 stopwatch.pause();
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::grab_dimension( const std::string& key, int* ndim, int* dims  ){
  int parlevel=0;
  std::vector<std::string> c;
  c = Tools::getWords(key," ,\n\t",&parlevel);
  if( c[0]=="positions" ){
      *ndim=2; dims[0]=atoms.getNatoms(); dims[1]=3;
  } else if (c[0]=="vessels"){
     *ndim=2;
     string a=c[1];
     vesselbase::ActionWithVessel* mves= this->getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(a);
     if(!mves) plumed_merror("action labelled " +  a + " does not exist or does not have vessels");
     multicolvar::MultiColvarBase* mycolv;
     mycolv = dynamic_cast<multicolvar::MultiColvarBase*>( mves );
     if(!mycolv) plumed_merror("action labeled " + mycolv->getLabel() + " is not a multicolvar");
     if( mycolv->weightWithDerivatives() ){
        dims[1]=mycolv->getNumberOfQuantities();
     } else {
        dims[1]=mycolv->getNumberOfQuantities()-1;
     }
     dims[0]=mycolv->getFullNumberOfTasks();
  } else {
     *ndim=1;
     // Initialize dimension to zero
     dims[0]=0;
     for(unsigned i=0;i<c.size();i++){
        std::size_t dot=c[i].find_first_of('.');
        string a=c[i].substr(0,dot);
        string name=c[i].substr(dot+1);
        // if it contains a dot:
        if(c[i].find(".")!=string::npos){
         	// Take all components
         	if ( name=="*" ){
              ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
              dims[0] += action->getNumberOfComponents();
         	// Take values with a specific name
           } else {
              ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
              if(!action){
         	      std::string str=" (hint! the actions in this ActionSet are: "; 
         	      str+=this->getActionSet().getLabelList()+")";
         	      plumed_merror("cannot find action named " + a +str);
         	   } 
         	   if( !(action->exists(c[i])) ){
         	      std::string str=" (hint! the components in this actions are: "; 
         	      str+=action->getComponentsList()+")";
         	      plumed_merror("action " + a + " has no component named " + name + str);
         	   } ;
              dims[0] += 1;
           }
        // if it doesn't contain a dot:
        } else {
           ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
           if(!action){
         	   std::string str=" (hint! the actions in this ActionSet are: "; 
         	   str+=this->getActionSet().getLabelList()+")";
         	   plumed_merror("cannot find action named " + c[i] + str );
         	}
         	if( !(action->exists(c[i])) ){
         	   std::string str=" (hint! the components in this actions are: "; 
         	   str+=action->getComponentsList()+")";
         	   plumed_merror("action " + c[i] + " has no component named " + c[i] +str);
         	};
           if (action->getNumberOfComponents()!=1) {
         	   plumed_merror("action " + c[i] + " has components! Specify one of them.");
           }
           dims[0] += 1; 
        }
     }
  }
}

void PlumedMain::grab_data( const std::string& key, void* outval ){
  double* p;
  if( atoms.getRealPrecision()==sizeof(double)) p=static_cast<double*>( outval );
  else if( atoms.getRealPrecision()==sizeof(float)) plumed_merror("yeah this needs fixing");    //float* p=static_cast<float*>( outval );
  else plumed_merror("Unknown real precision type");

  int parlevel=0;
  std::vector<std::string> c;
  c = Tools::getWords(key," ,\n\t",&parlevel);

  if( key=="positions" ){
     for(unsigned i=0;i<64*3;++i) p[i]=static_cast<double>(i); 
  } else if (c[0]=="vessels"){
     string a=c[1];
     // Assign action
     vesselbase::ActionWithVessel* mves= this->getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(a);
     if(!mves) plumed_merror("action labelled " +  a + " does not exist or does not have vessels");
     // Create multicolvar
     multicolvar::MultiColvarBase* mycolv;
     mycolv = dynamic_cast<multicolvar::MultiColvarBase*>( mves );
     // Vessel base to store values of multicolvar
     vesselbase::StoreDataVessel* stash=dynamic_cast<vesselbase::StoreDataVessel*>( mves->buildDataStashes( false, 0.0, NULL )  );
     std::vector<double> cvals( mycolv->getNumberOfQuantities() );
     // Calculate
     mves->calculate();
     unsigned k=0;
     for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i){
        stash->retrieveValue( i, true, cvals );
        if( mycolv->weightWithDerivatives() ){
          for(unsigned j=0;j<cvals.size();++j) {
             p[k]= static_cast<double>(cvals[j]);
             k++;
          }
        } else {
          for(unsigned j=1;j<cvals.size();++j) {
             p[k]= static_cast<double>(cvals[j]);
             k++;
          }
        }
     }  
  } else {
     unsigned offsetDim=0;
     for(unsigned i=0;i<c.size();i++){
        std::size_t dot=c[i].find_first_of('.');
        string a=c[i].substr(0,dot);
        string name=c[i].substr(dot+1);
        // If it contains a dot
        if(c[i].find(".")!=string::npos){    // if it contains a dot:
           // Take components from all actions with a specific name
           if ( name=="*" ){
              ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
              for(unsigned i=0;i<action->getNumberOfComponents();++i) p[offsetDim+i]=static_cast<double> (action->getOutputQuantity(i) ) ; 
              offsetDim += action->getNumberOfComponents();
           // Take values with a specific name
           } else {
              ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
              p[offsetDim]=static_cast<double> (action->getOutputQuantity(name));
              offsetDim += 1;
           }
        } else {
           ActionWithValue* action=this->getActionSet().selectWithLabel<ActionWithValue*>(a);
           p[offsetDim]=static_cast<double> (action->getOutputQuantity(0) ) ; 
           offsetDim += 1;
        } 
     }
  }
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::init(){
// check that initialization just happens once
  initialized=true;
  atoms.init();
  if(!log.isOpen()) log.link(stdout);
  log<<"PLUMED is starting\n";
  log<<"Version: "<<config::getVersionLong()<<" (git: "<<config::getVersionGit()<<") compiled on " __DATE__ " at " __TIME__ "\n";
  log<<"Please cite this paper when using PLUMED ";
  log<<cite("Tribello, Bonomi, Branduardi, Camilloni, and Bussi, Comput. Phys. Commun. 185, 604 (2014)");
  log<<"\n";
  log<<"For further information see the PLUMED web page at http://www.plumed-code.org\n";
  log.printf("Molecular dynamics engine: %s\n",MDEngine.c_str());
  log.printf("Precision of reals: %d\n",atoms.getRealPrecision());
  log.printf("Running over %d %s\n",comm.Get_size(),(comm.Get_size()>1?"nodes":"node"));
  log<<"Number of threads: "<<OpenMP::getNumThreads()<<"\n";
  log<<"Cache line size: "<<OpenMP::getCachelineSize()<<"\n";
  log.printf("Number of atoms: %d\n",atoms.getNatoms());
  if(grex) log.printf("GROMACS-like replica exchange is on\n");
  log.printf("File suffix: %s\n",getSuffix().c_str());
  if(plumedDat.length()>0){
    readInputFile(plumedDat);
    plumedDat="";
  }
  atoms.updateUnits();
  log.printf("Timestep: %f\n",atoms.getTimeStep());
  if(atoms.getKbT()>0.0)
    log.printf("KbT: %f\n",atoms.getKbT());
  else {
    log.printf("KbT has not been set by the MD engine\n");
    log.printf("It should be set by hand where needed\n");
  }
  log<<"Relevant bibliography:\n";
  log<<citations;
  log<<"Please read and cite where appropriate!\n";
  log<<"Finished setup\n";
}

void PlumedMain::readInputFile(std::string str){
  plumed_assert(initialized);
  log.printf("FILE: %s\n",str.c_str());
  IFile ifile;
  ifile.link(*this);
  ifile.open(str);
  std::vector<std::string> words;
  while(Tools::getParsedLine(ifile,words) && words[0]!="ENDPLUMED") readInputWords(words);
  log.printf("END FILE: %s\n",str.c_str());
  log.flush();	

  pilots=actionSet.select<ActionPilot*>();
}

void PlumedMain::readAction(std::string str) {
  plumed_assert(initialized);
  int parlevel=0;
  std::vector<std::string> words_tmp;
  words_tmp = Tools::getWords(str,NULL,&parlevel);
  readInputWords(words_tmp);
}

void PlumedMain::readInputWords(const std::vector<std::string> & words){
  plumed_assert(initialized);
  if(words.empty())return;
  else if(words[0]=="ENDPLUMED") return;
  else if(words[0]=="_SET_SUFFIX"){
    plumed_assert(words.size()==2);
    setSuffix(words[1]);
  } else {
    std::vector<std::string> interpreted(words);
    Tools::interpretLabel(interpreted);
    Action* action=actionRegister().create(ActionOptions(*this,interpreted));
    if(!action){
      log<<"ERROR\n";
      log<<"I cannot understand line:";
      for(unsigned i=0;i<interpreted.size();++i) log<<" "<<interpreted[i];
      log<<"\n";
      exit(1);
    };
    action->checkRead();
    actionSet.push_back(action);
  };

  pilots=actionSet.select<ActionPilot*>();
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::exit(int c){
  comm.Abort(c);
}

Log& PlumedMain::getLog(){
  return log;
}

void PlumedMain::calc(){
  prepareCalc();
  performCalc();
}

void PlumedMain::prepareCalc(){
  prepareDependencies();
  shareData();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// here we have the main steps in "calc()"
// they can be called individually, but the standard thing is to
// traverse them in this order:
void PlumedMain::prepareDependencies(){

  stopwatch.start("1 Prepare dependencies");

// activate all the actions which are on step
// activation is recursive and enables also the dependencies
// before doing that, the prepare() method is called to see if there is some
// new/changed dependency (up to now, only useful for dependences on virtual atoms,
// which can be dynamically changed).

// First switch off all actions
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
     (*p)->deactivate();
     (*p)->clearOptions();
  }

// for optimization, an "active" flag remains false if no action at all is active
  active=false;
  for(unsigned i=0;i<pilots.size();++i){
    if(pilots[i]->onStep()){
      pilots[i]->activate();
      active=true;
     }
  };

// also, if one of them is the total energy, tell to atoms that energy should be collected
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive()){
      if((*p)->checkNeedsGradients()) (*p)->setOption("GRADIENTS");
    }
  }

  stopwatch.stop("1 Prepare dependencies");
}

void PlumedMain::shareData(){
// atom positions are shared (but only if there is something to do)
  if(!active)return;
  stopwatch.start("2 Sharing data");
  if(atoms.getNatoms()>0) atoms.share();
  stopwatch.stop("2 Sharing data");
}

void PlumedMain::performCalc(){
  waitData();
  justCalculate();
  justApply();
}

void PlumedMain::waitData(){
  if(!active)return;
  stopwatch.start("3 Waiting for data");
  if(atoms.getNatoms()>0) atoms.wait();
  stopwatch.stop("3 Waiting for data");
}

void PlumedMain::justCalculate(){
  if(!active)return;
  stopwatch.start("4 Calculating (forward loop)");
  bias=0.0;
  work=0.0;

  int iaction=0;
// calculate the active actions in order (assuming *backward* dependence)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    std::string actionNumberLabel;
    if(detailedTimers){
      Tools::convert(iaction,actionNumberLabel);
      actionNumberLabel="4A "+actionNumberLabel+" "+(*p)->getLabel();
      stopwatch.start(actionNumberLabel);
    }
    ActionWithValue*av=dynamic_cast<ActionWithValue*>(*p);
    ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(*p);
    {
      if(av) av->clearInputForces();
      if(av) av->clearDerivatives();
    }
    {
      if(aa) aa->clearOutputForces();
      if(aa) if(aa->isActive()) aa->retrieveAtoms();
    }
    if((*p)->isActive()){
      if((*p)->checkNumericalDerivatives()) (*p)->calculateNumericalDerivatives();
      else (*p)->calculate();
      // This retrieves components called bias 
      if(av) bias+=av->getOutputQuantity("bias");
      if(av) work+=av->getOutputQuantity("work");
      if(av)av->setGradientsIfNeeded();	
      ActionWithVirtualAtom*avv=dynamic_cast<ActionWithVirtualAtom*>(*p);
      if(avv)avv->setGradientsIfNeeded();	
    }

    if(detailedTimers) stopwatch.stop(actionNumberLabel);
    iaction++;
  }
  stopwatch.stop("4 Calculating (forward loop)");
}

void PlumedMain::justApply(){
  
  if(!active)return;
  int iaction=0;
  stopwatch.start("5 Applying (backward loop)");
// apply them in reverse order
  for(ActionSet::reverse_iterator p=actionSet.rbegin();p!=actionSet.rend();++p){
    if((*p)->isActive()){

      std::string actionNumberLabel;
      if(detailedTimers){
        Tools::convert(iaction,actionNumberLabel);
        actionNumberLabel="5A "+actionNumberLabel+" "+(*p)->getLabel();
        stopwatch.start(actionNumberLabel);
      }

      (*p)->apply();
      ActionAtomistic*a=dynamic_cast<ActionAtomistic*>(*p);
// still ActionAtomistic has a special treatment, since they may need to add forces on atoms
      if(a) a->applyForces();

      if(detailedTimers) stopwatch.stop(actionNumberLabel);
    }
    iaction++;
  }

// this is updating the MD copy of the forces
  if(detailedTimers) stopwatch.start("5B Update forces");
  if(atoms.getNatoms()>0) atoms.updateForces();
  if(detailedTimers) stopwatch.stop("5B Update forces");

  if(detailedTimers) stopwatch.start("5C Update");
// update step (for statistics, etc)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive() && (*p)->checkUpdate()) (*p)->update();
  }
  if(detailedTimers) stopwatch.stop("5C Update");
// Check that no action has told the calculation to stop
  if(stopNow){
     if(stopFlag) (*stopFlag)=1;
     else plumed_merror("your md code cannot handle plumed stop events - add a call to plumed.comm(stopFlag,stopCondition)");
  }  
  stopwatch.stop("5 Applying (backward loop)");

// flush by default every 10000 steps
// hopefully will not affect performance
  if(step%10000==0){
    fflush();
    log.flush();
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p) (*p)->fflush();
  }
}

void PlumedMain::load(const std::string& ss){
  if(DLLoader::installed()){
     string s=ss;
     size_t n=s.find_last_of(".");
     string extension="";
     string base=s;
     if(n!=std::string::npos && n<s.length()-1) extension=s.substr(n+1);
     if(n!=std::string::npos && n<s.length())   base=s.substr(0,n);
     if(extension=="cpp"){
       string cmd="plumed mklib "+s;
       log<<"Executing: "<<cmd;
       if(comm.Get_size()>0) log<<" (only on master node)";
       log<<"\n";
       if(comm.Get_rank()==0) system(cmd.c_str());
       comm.Barrier();
       base="./"+base;
     }
     s=base+"."+config::getSoExt();
     void *p=dlloader.load(s);
     if(!p){
       log<<"ERROR\n";
       log<<"I cannot load library "<<ss<<"\n";
       log<<dlloader.error();
       log<<"\n";
       this->exit(1);
     }
     log<<"Loading shared library "<<s.c_str()<<"\n";
     log<<"Here is the new list of available actions\n";
     log<<actionRegister();
  } else plumed_merror("loading not enabled, please recompile with -D__PLUMED_HAS_DLOPEN");
}

double PlumedMain::getBias() const{
  return bias;
}

double PlumedMain::getWork() const{
  return work;
}

FILE* PlumedMain::fopen(const char *path, const char *mode){
  std::string mmode(mode);
  std::string ppath(path);
  std::string suffix(getSuffix());
  std::string ppathsuf=ppath+suffix;
  FILE*fp=std::fopen(const_cast<char*>(ppathsuf.c_str()),const_cast<char*>(mmode.c_str()));
  if(!fp) fp=std::fopen(const_cast<char*>(ppath.c_str()),const_cast<char*>(mmode.c_str()));
  plumed_massert(fp,"file " + ppath + " cannot be found");
  return fp;
}

int PlumedMain::fclose(FILE*fp){
  return std::fclose(fp);
}

std::string PlumedMain::cite(const std::string&item){
  return citations.cite(item);
}

void PlumedMain::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    (*p)->flush();
  }
}

void PlumedMain::insertFile(FileBase&f){
  files.insert(&f);
}

void PlumedMain::eraseFile(FileBase&f){
  files.erase(&f);
}

void PlumedMain::stop(){ 
  stopNow=true;
}

void PlumedMain::runJobsAtEndOfCalculation(){
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
      (*p)->runFinalJobs();
  }
} 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
