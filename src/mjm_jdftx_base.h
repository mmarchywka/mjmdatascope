#ifndef MJM_TCL_JDFTX_BASE_H__
#define MJM_TCL_JDFTX_BASE_H__


#include "mjm_dft_collection_data.h"
#include <stdlib.h>
//#include <tcl.h>
//#include "mjm_tcl_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
// this should be moved and unified.
//#include <cblas.h>
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>

//#include "mjm_sample_spectrum.h"
#include "mjm_io.h"
// the crap that needs this shold be moved.
//#include "prob_current.h"
#include "mjm_interchanges.h"
/*
*/
class mjm_jdftx_typedefs
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;


static OsTy & info() { return std::cout; } 


};
class mjm_jdftx_file : public mjm_jdftx_typedefs
{
typedef mjm_jdftx_file Myt;

protected:
typedef mjm_collection_data MetaData;
typedef std::vector<StrTy> Labels;
typedef std::vector<D> Data;
typedef std::pair<Labels,Data> Line;
typedef std::vector<Line> Lines;

Lines m_data;
MetaData m_meta;
public:
IdxTy size() const { return m_data.size(); } 
const Data & get(const IdxTy atom) const { return m_data[atom].second; }
Data & get(const IdxTy atom) { return m_data[atom].second; }
// no idea why not addede before.
const Data & operator()(const IdxTy atom) const { return m_data[atom].second; }
const D & get(const IdxTy atom, const IdxTy coord) const
{ return m_data[atom].second[coord]; }
const StrTy & get_label(const IdxTy atom, const IdxTy coord) const
{ return m_data[atom].first[coord]; }

void perturb(const IdxTy atom, const IdxTy coord, const D delta)
{
m_data[atom].second[coord]+=delta;
}
void compile()
{
m_meta.find_neighbors(*this,2);
m_meta.find_local_coords(*this,0);
}
// this should point to a specific ion object otherwise hard to deref
class _iterator
{
typedef _iterator Myt;
IdxTy m_atom,m_coord,m_size,m_idx;
public:
_iterator(const mjm_jdftx_file & that)
:m_atom(0),m_coord(0), m_size(that.size()*3),m_idx(0){}
_iterator(const IdxTy natoms )
:m_atom(0),m_coord(0), m_size(natoms*3),m_idx(0){}
void clear() { m_atom=0; m_coord=0; m_idx=0; } 
Myt & operator++() {inc(); return *this; }
operator bool() const { return m_idx<m_size; } 
IdxTy operator *() const { return m_idx; } 
operator IdxTy() const { return m_idx; } 
IdxTy atom() const { return m_atom; }  
IdxTy coord() const { return m_coord; }
void atom_coord(IdxTy &a , IdxTy & c) const { a=m_atom; c=m_coord;}  
IdxTy idx() const {return m_idx;}
private:
void inc()
{
++m_coord; if ( m_coord==3) { m_coord=0; ++m_atom;}
++m_idx;

}
};

typedef _iterator iterator;

// hopefully inlining works LOL
const D distance(const iterator & a1 , const iterator &a2) const
{ return distance(a1.atom(),a2.atom()); }
const D distance(const IdxTy a1 , const IdxTy a2) const
{ return ::sqrt(distance(a1,a2)); }
const D distance2(const IdxTy a1 , const IdxTy a2) const
{
const Data & s1=m_data[a1].second;
const Data & s2=m_data[a2].second;
const D dx=s2[0]-s1[0];
const D dy=s2[1]-s1[1];
const D dz=s2[2]-s1[2];

return (dx*dx+dy*dy+dz*dz);

}
const D distance(const IdxTy a1 , const Myt & that, const IdxTy a2) const
{

//return ::sqrt(dx*dx+dy*dy+dz*dz);
return ::sqrt(distance2(a1,that,a2));

}
const D distance2(const IdxTy a1 , const Myt & that, const IdxTy a2) const
{
const Data & s1=m_data[a1].second;
const Data & s2=that.m_data[a2].second;
const D dx=s2[0]-s1[0];
const D dy=s2[1]-s1[1];
const D dz=s2[2]-s1[2];

//return ::sqrt(dx*dx+dy*dy+dz*dz);
return (dx*dx+dy*dy+dz*dz);

}




// this needs a better impl if done a lot not
// even sure that 2nd nn works right now and 
// users will fail for diatomics LOL
IdxTy nearest(const IdxTy n,const IdxTy nn) const
{ return m_meta.nearest(n,nn); }

void zero(Data & d, const IdxTy sz) const
{
d.clear();
while ( d.size()!=sz) d.push_back(0);
}
// move amount delta along direction coord for atom 
// the coords are defined by nearest neighbor dir, second nn
// projected into plane perp, and cross product of these
// this code needs to go somewhere else though
void perturb_natural(const IdxTy atom, const IdxTy coord, const D delta)
{
Data dn=natural(atom,coord);
Data & a=m_data[atom].second;
//info()<<MM_MARK<<" have vectors "<<atom<<" to "<<n1<<" and "<<n2<<CRLF;
mjm::add3(a,a,dn,delta,(D)1);

}

Data natural(const IdxTy atom, const IdxTy coord) const 
{ return m_meta.natural(atom,coord); }

IdxTy get_coord(const IdxTy atom, const Myt &  that) const
{ return m_meta.get_coord(atom,*this,that); }
// this needs to be resolved along the stuid unit vectors...

typedef std::vector<D> resolve_result;
template<class Tf> resolve_result resolve(
const IdxTy thisatom,  const IdxTy thiscoord, const IdxTy atom, const IdxTy coord, const Myt & iox, const Tf& fx, const Tf& fzed) const
{ // resolve force along the displacement direction
const Data & s1=m_data[thisatom].second;
const Data & s2=iox.m_data[thisatom].second;
const Data  v=natural(atom,coord);
const Data  ud1=natural(thisatom,thiscoord);
// this should be normalized and dotted but what the hexk...
const Data & f1=fx.get(atom);
const Data & f2=fzed.get(atom);
Data d1; zero(d1,3);
Data d2; zero(d2,3);
Data dsum; zero(dsum,3);
mjm::add3(d1,s1,s2,(D)1.0,(D)-1.0);
mjm::add3(d2,f1,f2,(D)1,(D)-1);
mjm::add3(dsum,f1,f2,(D)1,(D)1);
const D fqual=mjm::dot(d2,d2)/mjm::dot(dsum,dsum);
//info()<<MM_MARK<<" FFFFUUUU "<<atom<<" "<<((d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2]))<<" "<<((d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]))<<CRLF;
info()<<MM_MARK<<" FOQUALITUY atom="<<atom<<"."<<coord<<" indir= "<<mjm::dot(ud1,d1)<<" forcex= "<< fqual<<CRLF;

resolve_result rr;
const D dist=mjm::dot(ud1,d1);
// this relies on displacement ebing along a given direction 
//D wtf=::sqrt(mjm::dot(d2,d2));
//if (mjm::dot(d2,v)<0) wtf=-wtf;
//return wtf/::sqrt(mjm::dot(d1,d1));
//return mjm::dot(d2,v)/::sqrt(mjm::dot(d1,d1));

//return mjm::dot(d2,v)/(mjm::dot(ud1,d1));
// first the force constant,
rr.push_back(mjm::dot(d2,v)/(dist));
// distance
rr.push_back(dist);
// the background force amount
rr.push_back(mjm::dot(f2,f2));
// observed force dotted into background, 
rr.push_back(mjm::dot(d2,f2));
// raw force- 
rr.push_back(mjm::dot(f1,v));
// bg force in diorection of basis
rr.push_back(mjm::dot(f2,v));



return rr;

//return ::sqrt((d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2])/(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]));
}

void zero_means()
{
info()<<MM_MARK<<" setting means to zero a bit of a kkluge"<<CRLF;
const IdxTy sz=size();
if ( sz<1) return ;
iterator ii(*this);
Data d; zero(d,3);
while (ii)
{
mjm::add3(d,d,get(ii.atom()),(D)1,(D)1);
++ii;
}
d[0]=d[0]/sz;
d[1]=d[1]/sz;
d[2]=d[2]/sz;

info()<<MM_MARK<<"FORCE KLUGE OFFSETS "<<d[0]<<" "<<d[1]<<" "<<d[2]<<CRLF;
ii.clear();
while (ii)
{
Data & r=get(ii.atom());
mjm::add3(r,r,d,(D)-1.0,(D)1);
++ii;
}


}

};

 

class mjm_jdftx_forces : public mjm_jdftx_file
{
typedef mjm_jdftx_forces Myt;
public:
void load(IsTy * is)
{
IdxTy i=0;
D x,y,z,flag;
while ( !is->eof()&&is->good())
{
//(*is)>>x>>y;
bool add_this_one=true;
if ( !(!is->eof()&&is->good())) break;
if ( add_this_one)
{
StrTy ion,species;
(*is)>>ion;
if ( ion!="force") 
{
info()<<MM_MARK<<" ignoring "<<ion<<" which should be force ..."<<CRLF;
continue;
}
(*is)>>species>>x>>y>>z>>flag;
Line l;
l.first.push_back(ion);
l.first.push_back(species);
l.second.push_back(x);
l.second.push_back(y);
l.second.push_back(z);
l.second.push_back(flag);
m_data.push_back(l);
//std::cout<<MM_MARK<<" adding energfy "<<i<<" energ "<<e<<" "<<i<<CRLF;
//add(x,y);
++i;
} // add
}

}
template <class Os> void store(Os & os,const StrTy & hdr)
{
const IdxTy sz = m_data.size();
info()<<MM_MARK<<" storing size is "<<sz<<CRLF;
IdxTy i=0;
const StrTy sep=" ";
if ( hdr.length()!=0) os<<hdr<<CRLF;
while ( i<sz)
{
const Line & l= m_data[i];
os<<l.first[0]<<sep<<l.first[1];
os<<sep<<l.second[0]<<sep<<l.second[1]<<sep<<l.second[2]<<sep<<l.second[3];
os<<CRLF;
//os<<m_data[i][0]<<" "<<m_data[i][1]<<CRLF;
++i;
}
}




};

class mjm_jdftx_ionpos : public mjm_jdftx_file
{
typedef mjm_jdftx_ionpos Myt;
public:
const StrTy species(const IdxTy atom ) const { return m_data[atom].first[1]; } 
void load(IsTy * is)
{
IdxTy i=0;
D x,y,z,flag;
while ( !is->eof()&&is->good())
{
//(*is)>>x>>y;
bool add_this_one=true;
if ( !(!is->eof()&&is->good())) break;
if ( add_this_one)
{
StrTy ion,species;
(*is)>>ion;
if ( ion!="ion") 
{
info()<<MM_MARK<<" ignoring "<<ion<<" which should be ion..."<<CRLF;
continue;
}
(*is)>>species>>x>>y>>z>>flag;
Line l;
l.first.push_back(ion);
l.first.push_back(species);
l.second.push_back(x);
l.second.push_back(y);
l.second.push_back(z);
l.second.push_back(flag);
m_data.push_back(l);
//std::cout<<MM_MARK<<" adding energfy "<<i<<" energ "<<e<<" "<<i<<CRLF;
//add(x,y);
++i;
} // add
}

}
template <class Os> void store(Os & os,const StrTy & hdr)
{
const IdxTy sz = m_data.size();
info()<<MM_MARK<<" storing size is "<<sz<<CRLF;
IdxTy i=0;
const StrTy sep=" ";
const std::streamsize oldsz=os.precision(18);
if ( hdr.length()!=0) os<<hdr<<CRLF;
while ( i<sz)
{
const Line & l= m_data[i];
os<<l.first[0]<<sep<<l.first[1];
os<<sep<<l.second[0]<<sep<<l.second[1]<<sep<<l.second[2]<<sep<<l.second[3];
os<<CRLF;
//os<<m_data[i][0]<<" "<<m_data[i][1]<<CRLF;
++i;
}

os.precision(oldsz);
}




};

#if 0 
/*
template <class SigClass> class mjm_tcl_base


*/
#endif

#endif

