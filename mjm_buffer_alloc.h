#ifndef MJM_BUFFER_ALLOC_H__
#define MJM_BUFFER_ALLOC_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
#include "mjm_strings.h"

#include "mjm_canned_methods.h"
#include "mjm_collections.h"
#include "mjm_string_kvp.h"


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


// Fri Mar 17 20:51:56 EDT 2023
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_buffer_alloc   
// g++  -Wall -std=gnu++11 -DTEST_MJM_BUFFER_ALLOC -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_buffer_alloc.h  -o mjm_buffer_alloc.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_buffer_alloc("mjm_buffer_alloc"
, "  ");

template <class Tr>
class mjm_buffer_alloc 
{
 typedef mjm_buffer_alloc Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

typedef mjm_thread_util<Tr> ThreadTy;
typedef pthread_t ThreadId;

typedef mjm_string_kvp<Tr> StrKvp;
//typedef mjm_so_loader<Tr> Loader;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line ;
typedef mjm_string_base_params<Tr> BaseParams;
enum {BAD=~0};

typedef unsigned char Data;

class _os_ptr
{

public:
_os_ptr(): m_ptr(0), m_sz(BAD), m_pc(BAD),m_avail(0) {}
// LOL do not put in a vector etc 
~_os_ptr() {delete[] m_ptr; }
IdxTy avail() const { return m_avail;}
bool alloc(const IdxTy sz)
{
if (m_ptr) { MM_ERR(" realloc "<<MMPR(dump()))
delete [] m_ptr;
 }
m_ptr= new Data[sz];
m_sz=sz;
m_avail=m_sz;
m_pc=0;
return true;
} // alloc  
Data * next(const IdxTy sz) 
{if (m_avail<sz) return 0; 
Data * p= m_ptr+m_pc; 
m_pc+=sz;
m_avail-=sz;
return p; 
} // next
private:
_os_ptr operator=(const _os_ptr & x  ) const { return *this; } 
_os_ptr(const _os_ptr & x  ) {}
Data * m_ptr;
IdxTy m_sz;
IdxTy m_pc;
IdxTy m_avail;

}; // _os_ptr;
typedef _os_ptr Ptr;


class _user_block
{

public:

Data * m_ptr;
IdxTy m_size;
Ptr * m_os_block;


}; // _user_block

typedef _user_block UserPtr;

// serial number 
typedef std::map<IdxTy, Ptr> OSBlocks;
// idx is size 
typedef std::map<IdxTy, std::vector<UserPtr > > UserBlocks;
typedef std::map<Data* , UserPtr > UsedBlocks;


public:
mjm_buffer_alloc() {m_mutex_vector.resize(1);}
~mjm_buffer_alloc() {}
Data * request(const IdxTy sz) 
{
Data * p=0; 
EnterSerial(0);
auto ii=m_free.lower_bound(sz);
if (ii==m_free.end()) NewFreeBlock(sz);
//if there is no free block, see if there is one in os block.. 
ii=m_free.lower_bound(sz);
p=0;

// this must be of non-zero size.. 
auto & ve=(*ii).second;
p=ve.back().p;
 m_inuse[p]=ve.back();;
ve.pop_back();
if (ve.size()==0) m_free.erease(ii); 
ExitSerial(0);

return p;
} // request
// call with mutex 
void NewFreeBlock(const IdxTy sz)
{
Data * p=0; 
if (!m_os.size()) AddOS();
// hpefully this works lol
p = m_os[m_os.size()-1].next(sz);
if (p==0) { AddOS(); p = m_os[m_os.size()-1].next(sz); } 
UserPtr up;
up.p=p;
up.sz=sz;
up.m_os_block= &m_os[m_os.size()-1]; 
m_free[sz].push_back();

} // NewFreeBlock 


IdxTy release(Data * p)
{
EnterSerial(0);
const IdxTy sz=m_inuse[p];
if (sz==0) MM_ERR(" try to free zero size or missing pt "<<MMPR(p))
m_free[sz].push_back(p);
ExitSerial(0);
return 0;
} // release
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
void AddOS()
{
// does thiw work??? 
m_os[m_os.size()].alloc(1<<20);

}
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
OSBlocks m_os;
UserBlocks  m_free;
UsedBlocks  m_inuse ;

}; // mjm_buffer_alloc

//////////////////////////////////////////////

template <class Tr>
class mjm_buffer_alloc_map : public std::map<typename Tr::StrTy, mjm_buffer_alloc< Tr > >  
{
 typedef mjm_buffer_alloc_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_buffer_alloc< Tr> >   Super;
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
mjm_buffer_alloc_map() {}
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

}; // mjm_buffer_alloc_map




////////////////////////////////////////////
#ifdef  TEST_MJM_BUFFER_ALLOC
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
typedef tester_< mjm_buffer_alloc <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_BUFFER_ALLOC "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_buffer_alloc<Tr>  Myt;
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

#endif // MJM_BUFFER_ALLOC_H__ 
