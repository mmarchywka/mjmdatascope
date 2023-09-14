#ifndef MJM_DFT_COLLECTION_DATA_H__
#define MJM_DFT_COLLECTION_DATA_H__

/*
En masse, calculate a bunch of things about a collection
such as a moelcule. This included nearest neighors and local coord
systems. 

*/

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
#include <algorithm>

//#include "mjm_sample_spectrum.h"
#include "mjm_io.h"
// the crap that needs this shold be moved.
//#include "prob_current.h"
#include "mjm_interchanges.h"
/*
*/
class mjm_collection_data_typedefs
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

/*
// read mostly matrix
// the vectors are used for convenience here but terrible
// for everything else
template <class Td, int N, int M > class  rm_matrix
{
typedef unsigned int IdxTy;
typedef Td D;
enum { m_n=N, m_m=M};
public:

};

*/



// this is only done once, brute force is ok even over 
// long molecules
template <class To> class distance_compare : 
public mjm_collection_data_typedefs
{
typedef To Collection;
typedef D VTy;
VTy m_origin[3];
const Collection & m_tgt;
template <class Tx> D r( const Tx & x) const
{
const D dx=x[0]-m_origin[0];
const D dy=x[1]-m_origin[1];
const D dz=x[2]-m_origin[2];
return dx*dx+dy*dy+dz*dz;
}
template <class Tx, class Ty> D dx( const Tx & x, const Ty & y ) const
{
return  0; 
}
public:
// these will be atoms or actually we should be sorting INDEXES
// and be bound to some object.  
distance_compare(const D & x, const D& y, const D & z, const Collection & c)
:m_tgt(c)
{m_origin[0]=x; m_origin[1]=y; m_origin[2]=z;} 
distance_compare(const IdxTy  i, const Collection & c) : m_tgt(c)
{m_origin[0]=c(i)[0]; m_origin[1]=c(i)[1]; m_origin[2]=c(i)[2]; } 



int operator()(const IdxTy & a, const IdxTy & b) const
{
const D del=r(m_tgt.get(a))-r(m_tgt.get(b));
if ( del>0) return 1;
if (del<0) return -1;
return 0; 

}

};  //disance_compare

class mjm_collection_data : public mjm_collection_data_typedefs
{
typedef mjm_collection_data Myt;
protected:

// too much effort for now... 
//typedef rm_matrix<D,3,3> CoordSystem;
//typedef std::vector<D> Array;

//typedef std::vector<Array> MatTy;

typedef mjm_2D_data MatTy;
typedef MatTy::array_type Array;

typedef std::vector<IdxTy> VectorProxy;
typedef MatTy CoordSystem;
typedef std::vector<CoordSystem> LocalCoords;

typedef std::vector<IdxTy> Nearest;
typedef std::vector<Nearest> NearestVector;


NearestVector m_neighbors;
LocalCoords m_local_coords;


void Init(VectorProxy & vp, const IdxTy n)
{
vp.clear();
while (vp.size()<n) vp.push_back(vp.size());

}
public:
// find the first n nearest
template <class Ty> void find_neighbors(const Ty & list, const IdxTy n)
{
typedef distance_compare<Ty> Compare;

const IdxTy sz=list.size();
if ( sz<=n) info()<<MM_MARK<<" FATAL ERROR list is too small for needed neighbors"<<CRLF;
// stop with for, just use while itor
for(IdxTy i=0; i<sz ; ++i)
{ 
// i'm not sure this helps unless list items are large and vp
// contains the coordinates as now it has two lists to 
// deref on each compare, however idx+coords would be good
VectorProxy vp;
Init(vp,sz);
Compare c(i,list);
std::sort(vp.begin(), vp.end(),c);
Nearest nn;
IdxTy j=0;
// should eliminate self from nearest although it should be easier to 
// just bump by one. 
while ( nn.size()<n)
{
// should only be at zero 
while ( vp[j]==i) ++j;
info()<<MM_MARK<<"NAYBURS  for atom "<<i<<" adding nieghbor "<<j<<" as atom "<<vp[j]<<CRLF; 
 nn.push_back(vp[j]);
++j;
}
m_neighbors.push_back(nn);

} // i

}
const Array & natural(const IdxTy atom, const IdxTy coord) const
{
return m_local_coords[atom][coord];
}

// this is stupid need to remove for rm_matrix
void zero(Array & d, const IdxTy sz) const
{
d.clear();
while ( d.size()!=sz) d.push_back(0);
}

const IdxTy nearest(const IdxTy atom, const IdxTy nn) const
{ return m_neighbors[atom][nn]; } 
// move amount delta along direction coord for atom 
// the coords are defined by nearest neighbor dir, second nn
// projected into plane perp, and cross product of these
template <class Ty> void find_local_coords(const Ty & list, const IdxTy param=0)
{
const IdxTy sz=list.size();
// stop with for, just use while itor
// needs to determine if sz<3....
for(IdxTy i=0; i<sz ; ++i)
{ 
// create a local coord system for each 
const IdxTy n1=m_neighbors[i][0];
const IdxTy n2=m_neighbors[i][1];
//info()<<MM_MARK<<" mnaural coods for "<<atom<<" to "<<n1<<" and "<<n2<<CRLF;
// bneed to check for diatomic molecule  LOL
// want to use read mostly arrays here. .. 
const Array & s1=list.get(i);
// could use directly?
const Array & s2=list.get(n1);
// need to span 3D, even if planar molecule 
// so can't use 3 nn's...
const Array & s3=list.get(n2);
MatTy m(3,3);
//m.set_size(3,3);
// this is dumb, need a matrix class.
Array &d1=m(0); // zero(d1,3);
Array &d2=m(1); //zero(d2,3);
Array &d3=m(2);; // zero(d3,3);
mjm::add3(d1,s1,s2,1,-1);
mjm::add3(d2,s1,s3,1,-1);
mjm::normalize(d1);
mjm::cross(d3,d2,d1);
// now fix d2 to make all rothogoanl
mjm::cross(d2,d3,d1);
mjm::normalize(d2);
mjm::normalize(d3);
//m.push_back(d1); m.push_back(d2); m.push_back(d3);
//fix this 
//m(0)=d1; m(1)=d2; m(2)=d3;
m_local_coords.push_back(m);

}


}
// given an offset, find the higest vector dot
template < class Ty> IdxTy get_coord(const IdxTy atom, const Ty & i1, const Ty  & that) const
{
const Array & s1=i1.get(atom);
const Array & s2=that.get(atom);
Array d1; zero(d1,3);
mjm::add3(d1,s1,s2,(D)1.0,(D)-1.0);

IdxTy m=0; D mag=0,val=0;
for (IdxTy i=0; i<3; ++i)
{
const Array  v=natural(atom,i);
 D dd=mjm::dot(d1,v);

dd*=dd;
if ( dd>mag) { m=i; mag=dd;val=mjm::dot(d1,v); } 
}
info()<<MM_MARK<<" RESOLVE COORD m="<<m<<" val="<<val<<CRLF;

return m;
}
// this needs to be resolved along the stuid unit vectors...


};


#endif

