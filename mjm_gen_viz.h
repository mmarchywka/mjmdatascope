#ifndef MJM_GEN_VIZ_H__
#define MJM_GEN_VIZ_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
#include "mjm_generic_iterators.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>

namespace mjm_gen_viz_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel;
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // 


// define input and consistent inverse output formats for supporting
// data structures. This is csv with quote and escape supporting
// ivector elements with decorated integers or footnotes.  
template <class Tr=mjm_gen_viz_traits::Tr > class form_def
{
typedef form_def<Tr>  Myt;
//typedef collections_traits::Tr Tr;
protected:
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;


public:
form_def(): m_sep(",") {}
// saving operations 
StrTy operator()(const StrTy & s)const  { return StrTy("\"")+s+StrTy("\""); } 
StrTy operator()(const IdxTy & i)const  {Ss s; s<<i; return s.str(); } 


StrTy operator()(const std::vector<IdxTy> & v)const  
{Ss s;s<<"\""; MM_LOOP(ii,v) {if (ii!=v.begin()) s<<sep();   s<<(*ii); } s<<"\""; return s.str(); } 


StrTy operator()(const std::vector<IdxTy> & v, const std::vector<StrTy> & footnotes)const  
{Ss s;s<<"\""; 
const IdxTy szf=footnotes.size();
MM_SZ_LOOP(i,v,sz)
//MM_LOOP(ii,v) 
{if (i!=0) s<<sep();   s<<v[i]; if (i<szf) { s<<footnotes[i];}  } 
s<<"\""; return s.str(); } 



StrTy operator()(const int & i)const  {Ss s; s<<i; return s.str(); } 
 // loading ops
Myt &  operator()(StrTy & d, const StrTy & s)  { d=s; return *this; } 
Myt &  operator()(unsigned int  & d, const StrTy & s)  { d=myatoi(s); return *this; } 


Myt &  operator()(std::vector< unsigned int>   & d, const StrTy & s)  
{
StrTy r=imbue_taboo(s);
  Ss ss(r); 
while ( true) {StrTy x=""; ss>>x; if (x=="") break;   d.push_back( myatoi(x)) ; } return *this; } 

Myt &  operator()(std::vector< unsigned int>   & d, std::vector<StrTy> &footnotes, const StrTy & s)  
{
StrTy r=imbue_taboo(s);
  Ss ss(r); 
while ( true) {StrTy x=""; ss>>x; if (x=="") break;   
const IdxTy sz=x.length();
const char * p=x.c_str();
//char y[sz+1];
StrTy note="";
for (IdxTy i=0; i<sz; ++i) 
{
	const char c=p[i];  
	if (!::isdigit(c)&&(c!='-') ) 
	{ 
		note=StrTy(p+i); break; 
	}  
} 
d.push_back( myatoi(x)); 
footnotes.push_back(note);
}
 return *this; } 



Myt &  operator()( int  & d, const StrTy & s)  { d=myatoi(s); return *this; } 
Myt &  operator()( double  & d, const StrTy & s)  { Ss ss(s); ss>>d; return *this; } 

const StrTy&  sep() const { return m_sep; } 

StrTy imbue_taboo(const StrTy & s) const
{const IdxTy sz=s.length(); char r[sz+1]; r[sz]=0;  const char sep=m_sep.c_str()[0]; 
for(IdxTy i=0; i<sz; ++i) { char c=s[i]; if (c==sep) r[i]=' '; else r[i]=c; }
return StrTy(r);
}


int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }


void setup(CommandInterpretter & li) {
li.set_split(1,m_sep.c_str()[0]);

}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}
StrTy m_sep;
Ss m_ss;

}; // form_def

class mjm_output_warp
{
typedef mjm_output_warp Myt;
protected:
typedef mjm_gen_viz_traits::Tr  Tr;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

typedef std::map<IdxTy, IdxTy> level_map;
public:
mjm_output_warp():m_xpitch(1),m_ypitch(1), m_x_zed(0),m_y_zed(0),m_levels(0),m_max_x(1000),m_max_y(1000)  {}
mjm_output_warp(const D & xp, const D & yp, const D &xz, const D &yz)
:m_xpitch(xp),m_ypitch(yp), m_x_zed(xz),m_y_zed(yz),m_levels(0),m_max_x(1000),
m_max_y(1000)  {}
virtual ~mjm_output_warp() {}
void operator()(D & x, D & y, const IdxTy & xin, const IdxTy & yin)
{
(*this)(x,y,D(xin),D(yin));
}
void operator()(IdxTy & x, IdxTy & y, const IdxTy & xin, const IdxTy & yin)
{
D xx,yy;
(*this)(xx,yy,D(xin),D(yin));
x=IdxTy(xx);
y=IdxTy(yy);
}
// 2025 compiler coplains about hiding just a warnin fix later 
#if 0
// normally coord xform should be fast but here vitual is ok 
virtual void operator()(D & x, D & y, const D & xin, const D & yin)
{
x=xin*m_xpitch+m_x_zed;
y=yin*m_ypitch+m_y_zed;
}
virtual void operator()(D & xin, D & yin)
{
D x=xin*m_xpitch+m_x_zed;
D y=yin*m_ypitch+m_y_zed;
xin=x; yin=y;
}
#endif


void bounds(IdxTy  & x, IdxTy & y, const level_map & m)
{
D xx,yy;
bounds(xx,yy,m);
x=IdxTy(xx);
y=IdxTy(yy);
}//template <class Tm > 
virtual void bounds(D & x, D & y, const level_map & m)
{
m_levels=&m;
IdxTy maxsz=0;
IdxTy maxlevel=0;
//MM_LOOP(ii,lc.m_level_sz)
MM_LOOP(ii,m)
{
if ((*ii).first>maxlevel) maxlevel=(*ii).first;
if ((*ii).second>maxsz) maxsz=(*ii).second;
}
//y=maxsz; x=maxlevel;
(*this)(x,y,maxlevel,maxsz); 
m_max_x=maxlevel;
m_max_y=maxsz;
}

protected:
D  m_xpitch, m_ypitch, m_x_zed,m_y_zed;
const level_map* m_levels;
D m_max_x,m_max_y;
}; // mjm_output_warp

class mjm_circular_warp : public mjm_output_warp
{
typedef mjm_output_warp Super;
typedef mjm_circular_warp Myt;
public:
mjm_circular_warp(const D & xp, const D & yp, const D &xz, const D &yz)
:Super(xp,yp,xz,yz),m_xc(0),m_yc(0) {}
// warp into a circule with x as radius and y as some measure of the angle.
// need ymax(x) t get angle right 
// 2025 compiler complains fix when error 
#if 0
void operator()(IdxTy & x, IdxTy & y, const IdxTy & xin, const IdxTy & yin)
{

D xx,yy;
(*this)(xx,yy,D(xin),D(yin));
x=IdxTy(xx);
y=IdxTy(yy);
// Super()(x,y,xin,yin);

 } 
#endif
// 2025-08 compiler complains this hides 
// ../../mjm/hlib/mjm_gen_viz.h:180:14: warning: ‘virtual void mjm_output_warp::operator()(D&, D&)’ was hidden [-Woverloaded-virtual=]
#if 0 
virtual void operator()(D & x, D & y, const D & xin, const D & yin)
{
const D r=xin; // -m_xc; 
auto ii=(*m_levels).find(IdxTy(xin));
const D ymax=(ii!=(*m_levels).end())?(*ii).second:1; // (*m_levels)[IdxTy(xin)]; // .size();
const D theta= 2.0*M_PI*yin/(1.0+ymax);
x= m_xpitch*(r*::cos(theta)+m_xc); //   xin*m_xpitch+m_x_zed;
y=m_ypitch*(r*::sin(theta)+m_yc); // yin*m_ypitch+m_y_zed;
MM_ERR(MMPR4(x,y,m_xpitch,m_ypitch))
}
#endif

void bounds(IdxTy  & x, IdxTy & y, const level_map & m)
{ Super::bounds(x,y,m); } 
virtual void bounds(D & x, D & y, const level_map & m)
{
Super::bounds(x,y,m);
m_xc=m_max_x/2;
m_yc=m_max_y/2;
x=m_xpitch*m_max_x;
y=m_ypitch*m_max_y;
}
private:
D m_xc,m_yc;

}; // mjm_cirular_warp


class mjm_gen_viz 
{
typedef mjm_gen_viz Myt;
protected:
typedef mjm_gen_viz_traits::Tr  Tr;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;


public:
mjm_gen_viz(): m_dmel(0) {}
virtual ~mjm_gen_viz() {}

protected:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }

void DMel(const StrTy & e) { }
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}


Dmel * m_dmel;
Ss m_ss;
//To m_obj;

};  // mjm_collection


class tree_layout
{
typedef mjm_gen_viz_traits::Tr Tr;  
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
class pos { 
public: pos(const IdxTy _id, const IdxTy _x, const IdxTy _y)
{ id=_id; x=_x; y=_y; }
 pos() {}
IdxTy id,x,y; 
};
public:
typedef pos position;
void operator()(const std::vector<IdxTy> & v)
{
const IdxTy sz=v.size();
const IdxTy cnt=m_level_sz[sz];
const IdxTy id=v[sz-1];
++m_level_sz[sz];
m_layout[id]=pos(id,sz,cnt);
// this needs to be done when the parent is positioned
// or work solely on id
//if (sz>1) m_connect.push_back(
if (sz>1) m_parent[id]=v[sz-2];
}
std::map<IdxTy, pos> m_layout;
std::map<IdxTy , IdxTy> m_level_sz;
std::map<IdxTy , IdxTy> m_parent;

}; // tree_layout
#endif

