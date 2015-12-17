/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "ActionSet.h"
#include "ActionWithValue.h"
#ifdef __PLUMED_HAS_CREGEX 
#include <cstring>
#include "regex.h"
#endif

using namespace std;
namespace PLMD{

ActionSet::ActionSet(PlumedMain&p):
plumed(p){
  (void) plumed; // to suppress warning about "unused plumed"
}

ActionSet::~ActionSet()
{
  for(int i=size()-1;i>=0;i--) delete (*this)[i];
}

void ActionSet::clearDelete(){
  for(int i=size()-1;i>=0;i--) delete (*this)[i];
  clear();
}


std::string ActionSet::getLabelList() const{
  std::string outlist;
  for(const_iterator p=begin();p!=end();++p){
    outlist+=dynamic_cast<Action*>(*p)->getLabel()+" ";
  };
  return  outlist;
}

std::vector<std::string> ActionSet::getLabelVector() const{
  std::vector<std::string> outlist;
  for(const_iterator p=begin();p!=end();++p){
    outlist.push_back(dynamic_cast<Action*>(*p)->getLabel());
  };
  return  outlist;
}

void ActionSet::interpretArgumentList(const std::vector<std::string>& c, std::vector<Value*>&arg, Action* myaction ) const {

  for(unsigned i=0;i<c.size();i++){
      // is a regex? then just interpret it. The signal is () 
      std::size_t found1 = c[i].find("(");
      std::size_t found2 ;
      if(found1!=std::string::npos){
        found2=c[i].find(")",found1+1,1); // find it again
        if(found2!=std::string::npos){
                // start regex parsing
#ifdef __PLUMED_HAS_CREGEX 
                // take the string enclosed in quotes and put in round brackets 
                std::string myregex=c[i].substr(found1,found2-found1+1);
                if( myaction ) myaction->log.printf("  Evaluating regexp for this action: %s \n",myregex.c_str());
                int errcode;
                regex_t *preg = (regex_t*)malloc(sizeof(regex_t)); // pointer to the regular expression
                regmatch_t *pmatch;
                if ((errcode=regcomp(preg, myregex.c_str() ,REG_EXTENDED|REG_NEWLINE))) { // compile the regular expression
                        char* errbuf;
                        size_t errbuf_size;
                        // one can check the errors asking to regerror
                        errbuf_size = regerror(errcode, preg, NULL, 0);
                        if (!(errbuf=(char*)malloc(errbuf_size))) {
                                plumed_merror("cannot allocate the buffer for error detection in regexp!");
                        };
                        regerror(errcode, preg, errbuf, errbuf_size);
                        if( myaction ) myaction->error(errbuf);
                }
                plumed_massert(preg->re_nsub==1,"I can parse with only one subexpression");
                pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*preg->re_nsub);
                // select all the actions that have a value  
                std::vector<ActionWithValue*> all=select<ActionWithValue*>();
                if( all.empty() && myaction ) myaction->error("your input file is not telling plumed to calculate anything");
                for(unsigned j=0;j<all.size();j++){
                                std::string thisargument=all[j]->getLabel();
                                std::vector<std::string> ss=all[j]->getComponentsVector();
                                for(unsigned  k=0;k<ss.size();++k){
                                        thisargument=ss[k];
                                        unsigned ll=strlen(ss[k].c_str())+1;
                                        char*str;
                                        str=new char [ll];
                                        strcpy(str,ss[k].c_str());
                                        char *ppstr=str;
                                        if(!regexec(preg, ppstr , preg->re_nsub, pmatch, 0)) {
                                                if( myaction ) myaction->log.printf("  Something matched with \"%s\" : ",ss[k].c_str());
                                                do {
                                                                if (pmatch[0].rm_so != -1) {    /* The regex is matching part of a string */
                                                                        char *submatch;
                                                                        size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
                                                                        submatch = (char*)malloc(matchlen+1);
                                                                        strncpy(submatch, ppstr+pmatch[0].rm_so, matchlen+1);
                                                                        submatch[matchlen]='\0';
                                                                        if( myaction ) myaction->log.printf("  subpattern %s\n", submatch);
                                                                        // this is the match: try to see if it is a valid action 
                                                                        std::string putativeVal(submatch);
                                                                        if( all[j]->exists(putativeVal) ){
                                                                                arg.push_back(all[j]->copyOutput(putativeVal));
                                                                                if( myaction ) myaction->log.printf("  Action %s added! \n",putativeVal.c_str());
                                                                       }
                                                                        free(submatch);
                                                                };
                                                                ppstr += pmatch[0].rm_eo;       /* Restart from last match */
                                                } while(!regexec(preg,ppstr,preg->re_nsub,pmatch,0));
                                        }
                                        delete [] str;
                                }
                };
                regfree(preg);
                free(preg);
                free(pmatch);
#else
                plumed_merror("Regexp support not compiled!");
#endif
        }else{
                plumed_merror("did you want to use regexp to input arguments? enclose it between two round braces (...) with no spaces!");
        }
      }else{
        std::size_t dot=c[i].find_first_of('.');
        string a=c[i].substr(0,dot);
        string name=c[i].substr(dot+1);
        if(c[i].find(".")!=string::npos){    // if it contains a dot:
          if(a=="*" && name=="*"){
             // Take all values from all actions
             std::vector<ActionWithValue*> all=select<ActionWithValue*>();
             if( all.empty() && myaction ) myaction->error("your input file is not telling plumed to calculate anything");
             for(unsigned j=0;j<all.size();j++){
               for(int k=0;k<all[j]->getNumberOfComponents();++k) arg.push_back(all[j]->copyOutput(k));
             }
          } else if ( name=="*"){
             // Take all the values from an action with a specific name
             ActionWithValue* action=selectWithLabel<ActionWithValue*>(a);
             if(!action){
                   std::string str=" (hint! the actions in this ActionSet are: ";
                   str+=getLabelList()+")";
                   if( myaction ) myaction->error("cannot find action named " + a + str);
             }
             for(int k=0;k<action->getNumberOfComponents();++k) arg.push_back(action->copyOutput(k));
          } else if ( a=="*" ){
             // Take components from all actions with a specific name
             std::vector<ActionWithValue*> all=select<ActionWithValue*>();
             if( all.empty() && myaction ) myaction->error("your input file is not telling plumed to calculate anything");
             unsigned nval=0;
             for(unsigned j=0;j<all.size();j++){
                std::string flab; flab=all[j]->getLabel() + "." + name;
                if( all[j]->exists(flab) ){ arg.push_back(all[j]->copyOutput(flab)); nval++; }
             }
             if(nval==0 && myaction ) myaction->error("found no actions with a component called " + name );
          } else {
             // Take values with a specific name
             ActionWithValue* action=selectWithLabel<ActionWithValue*>(a);
             if(!action){
                   std::string str=" (hint! the actions in this ActionSet are: ";
                   str+=getLabelList()+")";
                   if( myaction ) myaction->error("cannot find action named " + a +str);
             }
             if( !(action->exists(c[i])) ){
                   std::string str=" (hint! the components in this actions are: ";
                   str+=action->getComponentsList()+")";
                   if( myaction ) myaction->error("action " + a + " has no component named " + name + str);
             } ;
             arg.push_back(action->copyOutput(c[i]));
          }
        } else {    // if it doesn't contain a dot
          if(c[i]=="*"){
             // Take all values from all actions
             std::vector<ActionWithValue*> all=select<ActionWithValue*>();
             if( all.empty() && myaction ) myaction->error("your input file is not telling plumed to calculate anything");
             for(unsigned j=0;j<all.size();j++){
               for(int k=0;k<all[j]->getNumberOfComponents();++k) arg.push_back(all[j]->copyOutput(k));
             }
          } else {
             ActionWithValue* action=selectWithLabel<ActionWithValue*>(c[i]);
             if(!action){
                   std::string str=" (hint! the actions in this ActionSet are: ";
                   str+=getLabelList()+")";
                   if( myaction ) myaction->error("cannot find action named " + c[i] + str );
             }
             if( !(action->exists(c[i])) ){
                   std::string str=" (hint! the components in this actions are: ";
                   str+=action->getComponentsList()+")";
                   if( myaction ) myaction->error("action " + c[i] + " has no component named " + c[i] +str);
             };
             arg.push_back(action->copyOutput(c[i]));
          }
        }
      }
  }
}




}
