#ifndef MJM_VERSATILE_POINT_H__
#define MJM_VERSATILE_POINT_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"
#include "mjm_worm_blob.h"
#include "mjm_collections.h"
//#include "mjm_tokenized_collections.h"
#include "mjm_canned_methods.h"


#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>

/*
@software{,
  author = {Michael J Marchywka},
  title = {},
abstract=(),
institution={},
license={Knowledge sir should be free to all },
publisher={Mike Marchywka},
email={marchywka@hotmail.com},
authorid={orcid.org/0000-0001-9237-455X},
  filename = {},
  url = {},
  version = {0.0.0},
  date-started = {}
}
*/

// Wed 12 Jul 2023 05:51:26 AM EDT
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_versatile_point   
// QUICKCOMPILE  g++  -Wall -Wno-misleading-indentation  -std=gnu++11 -DTEST_MJM_VERSATILE_POINT -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_versatile_point.h  -o mjm_versatile_point.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_versatile_point("mjm_versatile_point" , "  ");

template <class Tr,int Nc, int Nv>
class mjm_versatile_point 
{
 typedef mjm_versatile_point Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

// TYPEDEF 
enum { BAD=~0};
//typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;


// API

public:
D & x;
D & y;
// y will be invalide in Nc less than 2 but should compile :)
// partial specialization? 
mjm_versatile_point(): x(m_x[0]),y(m_x[1]) {}
mjm_versatile_point(const Myt & that): x(m_x[0]),y(m_x[1]) {
*this=that; }
Myt & operator=(const Myt & that ) 
{
memcpy(&m_x[0],&that.m_x[0],sizeof(m_x));
return *this; 
}
template <class Tt> Myt & operator=(const Tt & that ) 
{
x=that.x;
y=that.y;
//memcpy(&m_x[0],&that.m_x[0],sizeof(m_x));
return *this; 
}
~mjm_versatile_point() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
int myatoi(const StrTy & s )const   { return mjm_canned_methods::myatoi(s.c_str()); }
int myatoi(const char * c) const  { return mjm_canned_methods::myatoi(c); }
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


// MEMBERS
D m_x[Nc];
int m_v[Nv];

}; // mjm_versatile_point

//////////////////////////////////////////////

template <class Tr,int Nc, int Nv>
class mjm_versatile_point_map : public std::map<typename Tr::StrTy, mjm_versatile_point< Tr,Nc,Nv > >  
{
 typedef mjm_versatile_point_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_versatile_point< Tr,Nc,Nv> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_versatile_point_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }

private:

void Init()
{
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;


}
return ss.str();
// return Dump(flags); 

}




private:

}; // mjm_versatile_point_map




////////////////////////////////////////////
#ifdef  TEST_MJM_VERSATILE_POINT
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_versatile_point <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_VERSATILE_POINT "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_versatile_point<Tr>  Myt;
//Myt x(argc,args);
Myt x;

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_VERSATILE_POINT_H__ 
