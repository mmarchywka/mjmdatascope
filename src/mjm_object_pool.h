#ifndef MJM_OBJECT_POOL_H__
#define MJM_OBJECT_POOL_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

//#include "mjm_block_matrix.h"
#include "mjm_exbool_vector.h"
#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"
#include "mjm_worm_blob.h"
#include "mjm_collections.h"
//#include "mjm_tokenized_collections.h"
#include "mjm_canned_methods.h"

#include "mjm_pawnoff.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"
#include "mjm_generic_iterators.h"



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
<one line to give the program's name and a brief idea of what it does.>


Conceived and written by Mike Marchywka from 2023 to present.
See dates in individual code pieces as they were 
generated from my wizards. 
Copyright (C) <year> <name of author>


This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
@software{,
  author = {Michael J Marchywka},
  city = {Jasper GA 30143 USA},
  title = {},
abstract={},
institution={},
license={Knowledge sir should be free to all },
publisher={Mike Marchywka},
email={marchywka@hotmail.com},
authorid={orcid.org/0000-0001-9237-455X},
  filename={mjm_object_pool.h},
  url = {},
  version = {0.0.0},
  date-started={2023-12-24},
}
*/

// Sun 24 Dec 2023 08:17:46 PM EST
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_object_pool   
// QUICKCOMPILE  g++  -MMD -MF mjm_object_pool.deps  -Wall -Wno-misleading-indentation  -std=gnu++11 -DTEST_MJM_OBJECT_POOL -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_object_pool.h  -o mjm_object_pool.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_object_pool("mjm_object_pool" , "  ");

template <class Tr,class To=unsigned int >
class mjm_object_pool 
{
 typedef mjm_object_pool Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef To ObjectTy;
// TYPEDEF 
enum { BAD=~0};
enum { SZ=sizeof(ObjectTy) };
enum { RWLOCK=0,NLOCKS=3};
//typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_exbool_vector<Tr,IdxTy> AllocBits;

// API

public:
mjm_object_pool() {Init(); }
mjm_object_pool(const IdxTy n) {Init(); Alloc(n);}
mjm_object_pool(const StrTy & sin,const IdxTy flags) {Init(sin,flags); }
mjm_object_pool(const Ragged & r,const IdxTy start, const IdxTy first,const IdxTy flags ) {Init(r,start,first,flags);}
ObjectTy * operator[](const IdxTy n){ return m_ptr+n; }
const IdxTy size() const { return m_sz; }
void load(const StrTy & sin,const IdxTy flags) {Init(sin,flags); }
void load(const Ragged & r,const IdxTy start, const IdxTy first,const IdxTy flags ) {Init(r,start,first,flags);}
void set_size(const IdxTy n) { Alloc(n); } 
IdxTy alloc() { 
EnterSerial(RWLOCK);
IdxTy n=m_bits.first(false); 
if (n==BAD)
{
// Never hold locks during unprediclable operations
// like IO ... but in this case you have a fatal flaw alread 
// dump m_bits to a stringstgream.. 
MM_ERR(" rejected alloc "<<MMPR(m_bits.dump()))
ExitSerial(RWLOCK);
 return BAD;
}
m_bits.set(n,true);
ExitSerial(RWLOCK);
m_ptr[n]=ObjectTy();
// call dtor an dctor? lol 
return n; 
} // alloc 
IdxTy alloc(ObjectTy & x) { 

EnterSerial(RWLOCK);
IdxTy n=m_bits.first(false); 
//MM_ERR(MMPR2(__FUNCTION__,n))
if (n==BAD)
{
MM_ERR(" rejected alloc "<<MMPR(m_bits.dump()))
ExitSerial(RWLOCK);
 return BAD;
}
m_bits.set(n,true);
ExitSerial(RWLOCK);
//MM_ERR(MMPR(dump()))
// call dtor an dctor? lol 
m_ptr[n]=ObjectTy();
m_ptr[n]=x; 
return n; 
} // alloc 
void free(const IdxTy n) 
{
//MM_ERR(MMPR2(__FUNCTION__,n))
if (n==BAD) return; 
// call dtor? lol 
//m_ptr[n].~ObjectTy();
EnterSerial(RWLOCK);
m_bits.set(n,false);
ExitSerial(RWLOCK);
//MM_ERR(MMPR(dump()))
} // free


void save(const StrTy & fn,const StrTy &s) {Save(fn,s); }

// so now the problem of not calling object dtor lol 
// need to dtor if object has pointers etc 
~mjm_object_pool() {
//MM_DIE(" figure out dtor lol")
Free();
delete [] m_ptr; }

StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
int myatoi(const StrTy & s )const   { return mjm_canned_methods::myatoi(s.c_str()); }
int myatoi(const char * c) const  { return mjm_canned_methods::myatoi(c); }
static bool Bit(const IdxTy f, const IdxTy b) { return  ((f>>b)&1)!=0; }
// should loop over map now 
static void Set(IdxTy& f, const IdxTy b,const bool x) //const  
    { if (x) f|=(1<<b); else f&=((~1)<<b); }
StrTy Dump(const IdxTy flags=0) {Ss ss;  
EnterSerial(RWLOCK);
ss<<MMPR2(m_sz,m_ptr);
ss<<MMPR(m_bits.dump());
ExitSerial(RWLOCK);
return ss.str(); 
}
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


void Save(const StrTy & fn,const StrTy &s) {
// std::ofstream ofs(fn);

 } // Save
void Init(const Ragged & r, const IdxTy start=0, const IdxTy first=0, const IdxTy flags=0  )
{
Init();
const IdxTy sz=r.size();
for(IdxTy i=start; i<sz; ++i)
{
//const Line & l=r[i];
//const IdxTy len=l.size();


}  // i 

} // Init 
void Init(const StrTy  & sin,const IdxTy flags =0  )
{
Init();
BaseParams kvp(sin);
} // Init 

void Init()
{
m_ptr=0;
m_sz=0;
m_mutex_vector=MutexVector(NLOCKS);
} // Init
void Alloc(const IdxTy sz)
{
Free();
delete[] m_ptr;
m_ptr= new ObjectTy[sz];
// could call default ctor but typilcally user just
// assigns anyway 
//::memset(m_ptr,0,sz*SZ);
m_sz=sz;
m_bits=AllocBits(m_sz);
MM_ERR(MMPR2(m_sz,SZ))
}// Alloc
void Copy(const Myt & that )
{
Free();
delete[] m_ptr;
m_ptr= new ObjectTy[that.m_sz];
m_sz=that.sz;
//::memcpy(m_ptr,that.m_ptr,sz*SZ);
MM_ILOOP(i,m_sz){  m_ptr[i]=that.m_ptr[i];}
m_bits=that.m_bits;
}// Copy 
void Free()
{
MM_ILOOP(i,m_sz) 
{
if (m_bits[i]) m_ptr[i].~ObjectTy();
}
m_bits.clear();
} // Free 


// MEMBERS
ObjectTy * m_ptr;
IdxTy m_sz;
AllocBits m_bits;
}; // mjm_object_pool

//////////////////////////////////////////////

template <class Tr>
class mjm_object_pool_map : public std::map<typename Tr::StrTy, mjm_object_pool< Tr > >  
{
 typedef mjm_object_pool_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_object_pool< Tr> >   Super;
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
mjm_object_pool_map() {}
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

}; // mjm_object_pool_map




////////////////////////////////////////////
#ifdef  TEST_MJM_OBJECT_POOL
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
typedef tester_< mjm_object_pool <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_OBJECT_POOL "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_MSG(ss.str())

}
typedef mjm_ragged_table Ragged;
int main(int argc,char **args)
{
about();
typedef mjm_object_pool<Tr>  Myt;
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

if (cmd=="loadragged") {
 	const IdxTy start=atoi(cip.wif(2).c_str()); 
	const IdxTy first=atoi(cip.wif(3).c_str()); 
	const IdxTy flags=atoi(cip.wif(4).c_str()); 
Ragged r; r.load(cip.p1); x.load(r,start,first,flags); }
if (cmd=="load") {x.load(cip.p1,atoi(cip.p2.c_str())); }
if (cmd=="save") {x.save(cip.p1,cip.p2); }
if (cmd=="quit") break;
// NB this does not work in gneral when errors are disabled
//if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="dump") { auto wtf=x.dump();  MM_ERR(wtf) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_OBJECT_POOL_H__ 
