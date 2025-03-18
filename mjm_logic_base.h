#ifndef MJM_LOGIC_BASE_H__
#define MJM_LOGIC_BASE_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
// this requires c++11 or gnu++11 which conflicts with others.. 
//#include "mjm_rational.h"
#include "mjm_instruments.h"
//#include <algorithm>
#include <string>
//#include <streams>
#include <sstream>
/* skeletion copy and paste although better, add to a code geneartor with list of fields 
class fl_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
public:
fl_params( const StrTy & nm) : Super(nm) {}
fl_params() : Super() {}
//IdxTy output_offset() const { return m_map.get_uint("output_offset",output_freq()-1); } // freq-1;
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
//D Vt() const { return m_map.get_double("Vt",.0259); }
//bool enable_trapping() const { return m_map.get_bool("enable_trapping",!true); }

StrTy to_string(const StrTy & sep=" ") const
{
// should be geneated code, do not edit  wtf ?
Ss ss;
//ss<<"show_iter_status="<<show_iter_status()<<sep;
//ss<<"always_dump_state="<<always_dump_state()<<sep;
return ss.str();
}
}; // fl_params
*/



class mjm_logic_base 
{
public:
class Tr
{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef std::string StrTy;
typedef std::stringstream Ss;

};

typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef ReadWriteMap RWMap;
mjm_logic_base( const StrTy & nm) : m_map(nm) {}
mjm_logic_base() {} 

void load(const StrTy & nm) { m_map.load(nm); } 
void set(const StrTy & nm) { m_map.set_line(nm); } 
void set(const StrTy & nm, const D & val)
{
Ss ss;
ss.precision(64);
ss<<nm<<" "<<val;
set(ss.str());
}
// TODO FIXME 2017-09-09 why was this not here before?
// in theory this is not exact 
void set(const StrTy & nm, const StrTy & val)
{
Ss ss;
ss<<nm<<" "<<val;
set(ss.str());
}


void set(const StrTy & nm, const int  & val)
{
Ss ss;
ss<<nm<<" "<<val;
set(ss.str());
}
void set(const StrTy & nm, const unsigned int  & val)
{
Ss ss;
ss<<nm<<" "<<val;
set(ss.str());
}


// wtf? who dpeneds on this shot? 
//StrTy  get(const StrTy & nm) { return m_map.get_string(nm," "); } 
StrTy  get(const StrTy & nm) { return m_map.get_string(nm,""); } 

//template <class Tx> void config(Tx & dest, const StrTy & cmd, int  & i, int argc, char ** args)
 void config( const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
const StrTy nm=StrTy(args[i]);
MM_ERR(" loading file "<<nm<<" based on "<<s)
load(nm);
++i; // consume param and cmd
}

}
//void config_set( int  & i, int argc, char ** args)
template <class Tcmd> void config_set(const Tcmd & cmd,  int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
//if (s=="-set-param")
if (s==StrTy(cmd))
{
++i;
if (argc<=i) return;
const StrTy nm=StrTy(args[i]);
set(nm);
MM_ERR(" setting "<<nm<<" based on "<<s)
++i;
}

}
/*
const StrTy & operator[](const StrTy & key ) const
{
auto ii=m_map.find(key);
if (ii!=m_map.end()) return (*ii).second;
return m_default;
}
*/
 
template <class Tx> void get_all(Tx & keys, Tx & values) const
{ m_map.get_all(keys,values); }
StrTy dump()
{
Ss ss;
std::vector<StrTy> k,v;
get_all(k,v);
MM_SZ_LOOP(i,k,sz)
{
ss<<k[i]<<"="<<v[i]<<CRLF;

}
return ss.str();
}

RWMap & map() { return m_map;}
const RWMap & map() const  { return m_map;}

protected:
RWMap m_map;
StrTy m_default;
}; // mjm_logic_base 

#endif

