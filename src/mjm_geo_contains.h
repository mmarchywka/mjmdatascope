#ifndef MJM_GEO_CONTAINS_H__
#define MJM_GEO_CONTAINS_H__

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

// Tue 11 Jul 2023 04:21:02 AM EDT
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_geo_contains   
// QUICKCOMPILE  g++  -Wall -Wno-misleading-indentation  -std=gnu++11 -DTEST_MJM_GEO_CONTAINS -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_geo_contains.h  -o mjm_geo_contains.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_geo_contains("mjm_geo_contains" , "  ");

template <class Tr>
class mjm_geo_contains 
{
 typedef mjm_geo_contains Myt;
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
mjm_geo_contains() {}
~mjm_geo_contains() {}
template< class Tv> IdxTy coords(D & a,  D & b, const D & x, const D & y, const Tv & _v1,const Tv & _v2,const Tv & _v3, const IdxTy flags) const
{
IdxTy rc=0;
const bool rot1=Bit(flags,0);
// in the only case when this worked, a=-dleta, b=.5 the
// contains_only would have worked well. 
const bool rot2=Bit(flags,1);
// IdxTy _p1=p1;
// IdxTy _p2=p2;
// IdxTy _p3=p3;
Tv v1,v2,v3;
if (rot1) { v1=_v2; v2=_v3; v3=_v1; }
else if (rot2) { v1=_v3; v2=_v1; v3=_v2; }
else  { v1=_v1; v2=_v2; v3=_v3; }

D det=((v2.y-v3.y)*(v1.x-v3.x)+(v3.x-v2.x)*(v1.y-v3.y));
//m_det=1.0/((v2.y-v3.y)*(v1.x-v3.x)+(v3.x-v2.x)*(v1.y-v3.y));
if (det==0) { MM_ERR(" det is zero "<<MMPR(det))}
det=1.0/det;
 //{ MM_ERR(" det is zero "<<MMPR(m_det))}
// TODO these need to be correclated properly
// to avoid a+b= 1 + delta insteado fo - delta
D xa=x-v3.x;
D ya=y-v3.y;
// bad roundoff on a+b
//a=(v2.y-v3.y)*(x-v3.x)+(v3.x-v2.x)*(y-v3.y); // p1 
//b=(v3.y-v1.y)*(x-v3.x)+(v1.x-v3.x)*(y-v3.y); // p2 
// this helped but not the real solution... 
D a0y=v2.y-v3.y; //  
D a1y=v3.y-v1.y; //  
D a0x=v3.x-v2.x; //  
D a1x=(v1.x-v3.x); //  
//D anb=(v2.y-v3.y+v3.y-v1.y)*xa +(v3.x-v2.x+v1.x-v3.x)*ya; //  
//D amb=(v2.y-v3.y-v3.y+v1.y)*xa +(v3.x-v2.x-v1.x+v3.x)*ya; //  
D anb=(a0y+a1y)*xa+(a0x+a1x)*ya;
D amb=(a0y-a1y)*xa+(a0x-a1x)*ya;
a=.5*(anb+amb);
b=.5*(anb-amb);
a=a*det;
b=b*det;
return rc;
} // coords


// going from 1 to 2, b<0 p on RIGHT and b>0 p on LEFT
template< class Tv> D  contains_only( const D & x, const D & y, const Tv & v1,  const Tv & v2, const IdxTy flags) const
{
const D xp= v1.x*v2.y-v2.x*v1.y;
const D dx=v2.x-v1.x;
const D dy=v2.y-v1.y;
return xp+y*dx-x*dy;
} // contains_only
template< class Tp> IdxTy orient(Tp & p, const IdxTy flags) const
{
IdxTy rc=0;
const bool cw=Bit(flags,0);
D b1=contains_only(p[0].x,p[0].y,p[1],p[2],0);
if (cw) b1=0-b1;
// b1>0 ccw
if (b1>=0) return (b1>0)?0:1;
const auto p0=p[1];
p[1]=p[2]; p[2]=p0; 
return rc; 
} // orient 
template< class Tv> bool  contains_only( const D & x, const D & y, const Tv & v1, const Tv & v2, const Tv & v3, const IdxTy flags) const
{
const D b1=contains_only(x,y,v1,v2,flags);
const D b2=contains_only(x,y,v2,v3,flags);
const D b3=contains_only(x,y,v3,v1,flags);
MM_ERR(MMPR3(b1,b2,b3))
// point is on LEFT of all segments 
return (b1>=0)&&(b2>=0)&&(b3>=0);
} // contains_only 

template< class Tv> bool  contains( const D & x, const D & y, const Tv & v1, const Tv & v2, const Tv & v3, const IdxTy flags) const
{
const bool print_ab=!true; // Bit(flags,0);
//const bool print_ab= Bit(flags,0);
//const bool print_abo=!true; // Bit(flags,0);
const bool print_abo=!true; // Bit(flags,0);
//const bool print_ab_mid=Bit(flags,1);
D a,b;
coords(a,b,x,y,v1,v2,v3,0);
if (print_ab) {  MM_ERR(MMPR4(a,b,x,y)) }
if (print_ab) { contains_only(x,y,v1,v2,v3,flags); }

//if (print_ab_mid) if (( a>0)&&(a<1)&&(b>0)&&(b<1)){  MM_ERR(MMPR4(a,b,x,y)<<MMPR((a+b-1.0))) }
//const bool overall= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
//const bool pieces= ((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
//const bool sum= ((a+b)<=1); 
const bool overall1= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0);
// still needed, not sure anything helped lol 
//if (pieces&&!sum&&((a+b)<1.1))  coords(a,b,x,y,verts,1);
if (!overall1&&((a+b)<1.1))
 { coords(a,b,x,y,v1,v2,v3,1);
const bool overall2= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0);
if (!overall2&&((a+b)<1.1)) coords(a,b,x,y,v1,v2,v3,2);
}
//else return pieces;
const bool overall= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0);
if (overall&&(a>.1)&&(a<.9)) {
if (print_abo) {  MM_ERR(MMPR4(a,b,x,y)) }
if (print_abo) { contains_only(x,y,v1,v2,v3,flags); }

} 
return overall; 
//return ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
} // contains






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



}; // mjm_geo_contains

//////////////////////////////////////////////

template <class Tr>
class mjm_geo_contains_map : public std::map<typename Tr::StrTy, mjm_geo_contains< Tr > >  
{
 typedef mjm_geo_contains_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_geo_contains< Tr> >   Super;
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
mjm_geo_contains_map() {}
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

}; // mjm_geo_contains_map




////////////////////////////////////////////
#ifdef  TEST_MJM_GEO_CONTAINS
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
typedef tester_< mjm_geo_contains <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_GEO_CONTAINS "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_geo_contains<Tr>  Myt;
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

#endif // MJM_GEO_CONTAINS_H__ 