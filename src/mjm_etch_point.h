#ifndef MJM_ETCHPOINT_H__
#define MJM_ETCHPOINT_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

//#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
//#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"

#include "mjm_canned_methods.h"
//#include "mjm_dscope_dgram.h"
//#include "mjm_ff_band.h"
//#include "mjm_etch_curve.h"

//#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"

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



template <class Tr>
class mjm_etch_point
{
 typedef mjm_etch_point Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef std::vector<StrTy> Words;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////

typedef mjm_string_base_params<Tr> BaseParams;
// traits? 
enum {BAD=~0};
//typedef mjm_dscope_dgram<Tr> DgramID;
//typedef typename DgramID::data_type Data;
//typedef mjm_thread_util<Tr> ThreadTy;
//typedef pthread_t ThreadId;
//typedef typename ThreadTy::template ThParam<Myt> ThParam;

typedef Myt _Point;
typedef Myt P;
public:
// s is a parameterization... 
mjm_etch_point(): x(0),y(0),z(0),s(0),v(0) {}
mjm_etch_point(const D & _x, const D & _y): x(_x),y(_y),z(0),s(0),v(0) {}
D r2() const 
{ return x*x+y*y+z*z; } // r()
// fixed,  this is confusing as factors in exponents get messed up .
D r2( const _Point & that) const
{
_Point dr=that-(*this);
return dr.r2();
}
D r() const  { return sqrt(r2()); } 
D r( const _Point & that) const { return sqrt(r2(that)); } 
// actually cos or 1 if aligned, -1 if opposing.. 
// not strictly dot cos, need to respect orientation
// wrt normal.. 
D arc( const _Point & that,const _Point & n ) const 
{ 
const D  dir=(((*this)+that)*n);
// in a perfectly planar case, dir is zero but no ambiguity
if (dir==0) return 0; 

D dot=(*this)*that;
D a1=r();
D a2=that.r();
// don't bother on return code for now...
if ((a1==0)||(a2==0)) {MM_ERR(" zero in arc "<<MMPR4(dump(),that.dump(),a1,a2))}
// this is now between zero and one... 
 D xcos= .5*(1+dot/(a1*a2)); 
// dir<0 is an exposed peak 
return (dir<0)?(-xcos):xcos; 

}  // arc


// this and "that" is next ( ds>0) with n=z x s
_Point n( const _Point & that) const
{
_Point pn; // =*this;
pn.y= (that.x-x);
pn.x= (y-that.y);
const D r =(*this).r(that);
if (r<1e-10)
{
MM_ERR(" same point maybe "<<MMPR4(x,y,that.x,that.y))
}
pn.y=pn.y/r;
pn.x=pn.x/r;
if (std::isnan(pn.x)||std::isnan(pn.y))
{
MM_ERR(" git nan  "<<MMPR(r)<<MMPR4(x,y,that.x,that.y))


} 
return pn;
}

_Point   operator-(const _Point & that) const
{
_Point pn; // =*this;
pn.x=x-that.x;
pn.y=y-that.y;
pn.z=z-that.z;
pn.v=v-that.v;
return pn;
} // -


_Point   operator+(const _Point & that) const
{
_Point pn; // =*this;
pn.x=x+that.x;
pn.y=y+that.y;
pn.z=z+that.z;
pn.v=v+that.v;
return pn;
} // -


_Point &  operator+=(const _Point & that) //const
{
_Point& pn=*this;
pn.x=x+that.x;
pn.y=y+that.y;
pn.z=z+that.z;
pn.v=v+that.v;
//MM_ERR(MMPR3(pn.x,x,that.x))
//return pn;
return *this;
} // -=


_Point &  operator-=(const _Point & that) //const
{
_Point& pn=*this;
pn.x=x-that.x;
pn.y=y-that.y;
pn.z=z-that.z;
pn.v=v-that.v;

//return pn;
return *this;
} // -=

_Point   operator*(const D & that) const
{
_Point pn;
pn.x=x*that;
pn.y=y*that;
pn.z=z*that;
pn.v=v*that;

return pn;
} // * 

D operator*(const _Point & that) const
{
D dot=x*that.x;
dot+=y*that.y;
dot+=z*that.z;
return dot;
} // * 
// find the values for the cross of the two lines defeined by
// (p0,p1) and (p2,p3) with a and b the distance from P0 and P2
// in units of p1-p0 and p2-p3. If they are both between 0 and 1 
// segs cross
static const bool do_cross( const P & ij) 
 { 
const bool parallel= (std::isnan(ij.x)||std::isnan(ij.y));
if (parallel) return false;
return (ij.x>=0) &&(ij.x<=1) &&(ij.y>=0)&& (ij.y<=1); 
} // do_cross



static P cross(const P & p0, const P & p1, const P & p2, const P & p3) //const
{
/* 
const D b0=p3.x-p1.x;
const D b1=p3.y-p1.y;
const D a00=p0.x-p1.x;
const D a01=p2.x-p3.x;
const D a10=p0.y-p1.y;
const D a11=p2.y-p3.y;
*/

const D b0=p2.x-p0.x;
const D b1=p2.y-p0.y;
const D a00=p1.x-p0.x;
const D a01=p2.x-p3.x;
const D a10=p1.y-p0.y;
const D a11=p2.y-p3.y;


const D det=a00*a11-a10*a01;
// this should be due to parallel ie either no intersect
// or equal for now go ahead and return nan..
if (det==0) 
{
if (false) MM_ERR(" zero det "<<MMPR4(p0.dump(),p1.dump(),p2.dump(),p3.dump()))
}
//MM_ERR(MMPR(det))
P res((a11*b0-a01*b1)/det, (a00*b1-a10*b0)/det);
return res;
}// cross
StrTy dump() const
{
Ss ss;
ss<<MMPR4(x,y,z,s)<<MMPR(v);
return ss.str();
} // dump


D x,y,z,s,v;

}; //mjm_etch_point 


#endif // MJM_ETCHPOINT_H__ 
