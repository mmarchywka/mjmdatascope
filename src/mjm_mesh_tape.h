#ifndef MJM_MESH_TAPE_H__
#define MJM_MESH_TAPE_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_ff_band.h"
#include "mjm_string_kvp.h"
//#include "mjm_dscope_rawfifo.h"
//#include "mjm_ff_mesh_adapt.h"
//#include "mjm_etchfront.h"
//#include "mjm_pawnoff.h"
//#include "mjm_blob.h"
//#include "mjm_affine_mesh_hood.h"
#include "mjm_collections.h"
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
/*
2023-07-06  create from mjm_ff_phase_bound.h to pull out tape classes


*/




mjm_global_credits::credit __credit__mjm_mesh_tape("mjm_mesh_tape" , "  ");

template <class Tr> class mjm_tape_set
{
 typedef mjm_tape_set Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

// TYPEDEF 
typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line ;

//typedef mjm_ff_band<Tr> SurfaceBand;
//typedef mjm_fem_band<Tr> SurfaceBand;

//typedef mjm_pawnoff<Tr>  Hand;
//typedef mjm_blob<Tr>  Blob;

int myatoi(const char * c)  { return ::strtol(c,0,0); }
int myatoi(const StrTy & s) { return myatoi(s.c_str()); } 
static bool Bit(const IdxTy f, const IdxTy b)  { return  ((f>>b)&1)!=0; }


// API

public:
// description of "point tape" or line of points conformal
// to surface or anticipated surface 
typedef std::vector<D> Facs;
mjm_tape_set() : space(0),gspace(0), m_stagger(~0),m_flags(0) {} 
mjm_tape_set(const D & s, const D & g) : space(s),gspace(g),m_stagger(~0),m_flags(0)  {} 
mjm_tape_set(const D & s, const D & g,const D & f) 
	: space(s),gspace(g),m_stagger(~0), m_flags(0) {factors.push_back(f); } 
mjm_tape_set(const Line & l, const IdxTy start)  {
m_stagger= ~0; 
space=0;
gspace=0;
m_flags=0;
bool loading_fac=true;
const IdxTy sz=l.size();
if (sz>start) space=atof(l[start].c_str()); 
if (sz>(start+1)) gspace=atof(l[start+1].c_str()); 
for(IdxTy i=start+2; i<sz; ++i)
{
if (l[i]=="zero") { m_flags|=1; m_new_dofs.push_back(0); continue; }
if (l[i]=="stagger") 
	{ if ((i+1)<sz) m_stagger=myatoi(l[i+1]); ++i; continue; }
if (l[i]=="dofs") { m_flags|=1; loading_fac=false;  continue; }
if (l[i]=="dof") 
	{m_flags|=1; 
	if ((i+1)<sz) m_new_dofs.push_back(atof(l[i+1].c_str())); 
	++i; continue; }
if (loading_fac) factors.push_back(atof(l[i].c_str())); 
else m_new_dofs.push_back(atof(l[i].c_str())); 
} //i 
MM_ERR(MMPR2(__FUNCTION__,dump()))
}  // ctor

IdxTy stagger() const  { return m_stagger; }
IdxTy stagger(const IdxTy n ) const 
{ if (m_stagger!= IdxTy(~0)) return m_stagger;
return n;
}
IdxTy size() const { return factors.size(); }
void add(const D & f) { factors.push_back(f); } 
const Facs &  efactors() const { return factors; } 
const D &  efactor(const IdxTy i) const { return factors[i]; } 
const bool set_dofs() const { return Bit(m_flags,0); } 
D dof_val(const IdxTy i) const
{
const IdxTy sz=m_new_dofs.size();
if (i<sz) return m_new_dofs[i];
if (sz==1) return m_new_dofs[0];
return -1;
} // dof_val
StrTy dump() const
{
Ss ss;
ss<<MMPR4(space,gspace,m_stagger,m_flags);
MM_SZ_LOOP(i,factors,szi) { ss<<MMPR(factors[i]); }  
MM_SZ_LOOP(i,m_new_dofs,szfd) { ss<<MMPR(m_new_dofs[i]); }  
return ss.str();
} // dump 
D space,gspace;
IdxTy m_stagger;
IdxTy m_flags;
private:

Facs  factors,m_new_dofs;

}; // mjm_tape_set
//typedef _tape_set TapeSet;
template <class Tr> class mjm_ff_band_desc
{

 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line ;

//typedef mjm_ff_band<Tr> SurfaceBand;
//typedef mjm_fem_band<Tr> SurfaceBand;

//typedef mjm_pawnoff<Tr>  Hand;
//typedef mjm_blob<Tr>  Blob;

int myatoi(const char * c)  { return ::strtol(c,0,0); }
int myatoi(const StrTy & s) { return myatoi(s.c_str()); } 
static bool Bit(const IdxTy f, const IdxTy b)  { return  ((f>>b)&1)!=0; }




public:

mjm_ff_band_desc():m_f1(0),m_f2(0)  {}
mjm_ff_band_desc(const Line & l, const IdxTy start)  {
m_f1=0;
m_f2=0;
const IdxTy sz=l.size();
if (sz>start) m_f1=atof(l[start].c_str()); 
if (sz>(start+1)) m_f2=atof(l[start+1].c_str()); 
for(IdxTy i=start+2; i<sz; ++i)
{
//if (l[i]=="zero") { m_flags|=1; m_new_dofs.push_back(0); continue; }
//if (l[i]=="stagger") 
//	{ if ((i+1)<sz) m_stagger=myatoi(l[i+1]); ++i; continue; }
//if (l[i]=="dofs") { m_flags|=1; loading_fac=false;  continue; }
//if (l[i]=="dof") 
//	{m_flags|=1; 
//	if ((i+1)<sz) m_new_dofs.push_back(atof(l[i+1].c_str())); 
//	++i; continue; }
//if (loading_fac) factors.push_back(atof(l[i].c_str())); 
//else m_new_dofs.push_back(atof(l[i].c_str())); 
} //i 
MM_ERR(MMPR2(__FUNCTION__,dump()))
}  // ctor
StrTy dump() const { Ss ss; ss<<MMPR2(m_f1,m_f2); return ss.str(); } 
D m_f1,m_f2;

}; // mjm_ff_band_desc;
template <class Tr> class mjm_tape_vector
{
 typedef mjm_tape_vector Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

// TYPEDEF 
typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line ;

// TODO fix this and generalize 
// index src nbumber to set and index 
typedef std::pair<IdxTy, IdxTy> TapeIdx;
typedef std::vector<TapeIdx> Tapes;
typedef mjm_ff_band_desc<Tr> SurfaceBand;
typedef std::vector<SurfaceBand> SurfaceBands;

int myatoi(const char * c)  { return ::strtol(c,0,0); }
int myatoi(const StrTy & s) { return myatoi(s.c_str()); } 
static bool Bit(const IdxTy f, const IdxTy b)  { return  ((f>>b)&1)!=0; }


public:
typedef mjm_tape_set<Tr> V;
typedef std::vector<V> Vec;
typedef typename Vec::const_iterator Itor;
typedef SurfaceBand surface_band_type;
mjm_tape_vector(){  } 
mjm_tape_vector(const D space, const D gspace){ basic(space,gspace); } 
Itor begin() const { return m_v.begin();}
Itor end() const { return m_v.end();}
//const surface_band_type & operator()(const IdxTy n)const { return m_bands[n];}
const IdxTy srcmin() const { return 10; } 
const surface_band_type & band(const IdxTy n)const { return m_bands[n];}
const IdxTy bands()const { return m_bands.size();}
void clear() { m_bands.clear(); m_v.clear(); m_tapes.clear(); } 
IdxTy size() const { return m_v.size(); } 
void push_back( const V& x) { 
const IdxTy i=m_v.size();
MM_SZ_LOOP(j,x,xsz) { m_tapes.push_back(TapeIdx(i,j)); } 
m_v.push_back(x); 
}
IdxTy vfactor(const IdxTy i ) const
{
if ( i>=m_tapes.size()) return 0; 
const auto & idx=m_tapes[i];
return m_v[idx.first].efactor(idx.second); 
}
IdxTy dof_val(const IdxTy i ) const
{
if ( i>=m_tapes.size()) return 0; 
const auto & idx=m_tapes[i];
return m_v[idx.first].dof_val(idx.second); 
}
void dof_val(D & c, const IdxTy i ) const
{
if ( i>=m_tapes.size()) return ; 
const auto & idx=m_tapes[i];
if (! m_v[idx.first].set_dofs()) return ; 
//MM_ERR(" wtf "<<MMPR3(i,idx.first,idx.second))
//MM_SZ_LOOP(j,m_v,asd) { MM_ERR(" sizes "<<MMPR3(j,m_v[j].size(),m_v.size()))}
const D & _c= m_v[idx.first].dof_val(idx.second); 
if (_c<0) return ;
c=_c;
}


IdxTy load(const Ragged & r, const IdxTy flags=0)
{
MM_SZ_LOOP(i,r,smz)
{
const Line & l=r[i];
const IdxTy len=l.size();
if (len<2) continue;
if (l[0]=="band" ){  SurfaceBand sb(l,1); m_bands.push_back(sb);     continue;  }
if (l[0]!="tapeset" ) continue;
if ((len>3)&&(l[1]=="basic" )) 
	{ basic(myatoi(l[2]),myatoi(l[3])); continue;} 
push_back(V(l,1)); 

} // i
MM_ERR(" loaded ragged tape set "<<MMPR2(r.size(),size()))
MM_ERR(MMPR(dump()))
return 0;
} // load 
StrTy dump() const
{
Ss ss;
MM_SZ_LOOP(i,m_tapes,tsz) { ss<<" m_tapes "<<MMPR3(i,m_tapes[i].first, m_tapes[i].second); } 
MM_SZ_LOOP(i,m_v,vsz) { 
ss<<" m_v "<<MMPR2(i,m_v[i].dump());
 } 

return ss.str();
} // dump

void basic(const D space, const D gspace)
{
V x(space,gspace);
x.add(1.0);
x.add(2.0);
x.add(-1.0);
x.add(-2.0);
x.add(-4.0);
x.add(-6.0);
x.add(-8.0);
x.add(-10.0);
push_back(x);

} // basic
Vec m_v;
Tapes m_tapes;
SurfaceBands m_bands;
}; // mjm_tape_vector

#endif // MJM_MESH_TAPE_H__ 
