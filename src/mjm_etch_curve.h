#ifndef MJM_ETCHFCURVE_H__
#define MJM_ETCHFCURVE_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

//#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"

#include "mjm_canned_methods.h"
//#include "mjm_dscope_dgram.h"
//#include "mjm_ff_band.h"

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



template <class Tr,class _Point>
class mjm_etch_curve
{
 typedef mjm_etch_curve Myt;
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


public:
typedef _Point point_type;
mjm_etch_curve() : m_flags(0) {}
const bool closed() const { return Bit(m_flags,0); }
void closed() { m_flags|=1; }
void closed(const bool x ) {if (x)  m_flags|=1; else opened();  }
void opened() { m_flags&=~1; }
IdxTy size() const { return m_points.size(); } 
typedef std::vector<point_type> PointVec;
PointVec m_points;
PointVec m_normals; // should point in oriented dir
std::vector<D> m_radii; // Danger Will Robinson 1/r to avoid zero 
std::vector<D> m_arc; // amount of exposure  in cose theta ( 1 =buried )  
// needs to accomadate splits and joins ZZ
IdxTy m_flags;

// calculate radii from p and n ... 

D invrad(const IdxTy i,const IdxTy j) const
{

const D dotp=(m_normals[i]*m_normals[j]);
//const D theta=acos(m_normals[i]*m_normals[j]);
D theta=acos(dotp);
if (std::isnan(theta))
{
//MM_ERR(" bad there"<< MMPR3((m_normals[i]*m_normals[j]),m_normals[i].dump(),m_normals[j].dump())<<MMPR(m_points.size())<<MMPR4(i,j,m_normals.size(),m_radii.size()));
if ((dotp>=1) &&(dotp<1.1))  theta=0;
else
if ((dotp<=-1) &&(dotp> - 1.1))  theta=M_PI;

else 
{
MM_ERR(" bad there"<< MMPR3((m_normals[i]*m_normals[j]),m_normals[i].dump(),m_normals[j].dump())<<MMPR(m_points.size())<<MMPR4(i,j,m_normals.size(),m_radii.size()));

throw "doa";
} // else
} // isnan

return theta/m_points[j].r(m_points[i]); 
} // invrad
void recenter( const D & xc, const D & yc) {
D xcz=0;
D ycz=0;
MM_LOOP(ii,m_points)
{
xcz+=(*ii).x;
ycz+=(*ii).y;
} // ii 
xcz/=m_points.size();
ycz/=m_points.size();
// MM_ERR(" centering "<<MMPR4(xc,yc,xcz,ycz))
MM_LOOP(ii,m_points)
{
(*ii).x-=xcz+xc;
(*ii).y-=ycz+yc;

} // ii 
} // receneter
void makex()
{
const IdxTy n=200;
D f=6.28/n;

for(IdxTy i=0; i<n; ++i)
{
D x=1.0*i;
point_type p(x/n,2*sin(f*x));
m_points.push_back(p);
} // i 

} // make 
void make()
{
const IdxTy n=300;
D step=1.5/n;
for(int i=0; i<int(n); ++i)
{
D x=step*i;
bool ss=(x>.3)&&(x<.7);

D y=ss?(sqrt(.2*.2-(x-.5)*(x-.5))):0;
point_type p(x,y);
m_points.push_back(p);
} // i 

} // make
Myt & push_back(const D & x, const D & y)
{ m_points.push_back(point_type(x,y)); return *this; } 
void makesine()
{
const IdxTy n=200;
D f=6.28;
// need approx length of sine wave... 
D step=1.0/n;
D x=0;
D y=0;
const IdxTy leaders=10;
for(IdxTy i=0; i<leaders;  ++i)
{
D x=step*i;
D y=0;
point_type p(x,y);
m_points.push_back(p);

}
// need this crap for equal spacing.. 
for(IdxTy i=leaders; i<n; ++i)
{
x=step*i;
D y=sin(f*(x-step*leaders));
point_type p(x,y);
m_points.push_back(p);
} // i 
// TODO fails with zero points. 
const D xlast=m_points.back().x;
const D ylast=m_points.back().y;

for(IdxTy i=0; i<leaders;  ++i)
{
point_type p(xlast+step*(i+1),ylast);
m_points.push_back(p);

} // i 


} 


void assmake()
{
const IdxTy n=200;
D f=6.28;
// need approx length of sine wave... 
D stotal=0;
for(IdxTy i=0; i<n; ++i) { D dx=1.0/n; D dy=f*dx*cos(f*i); stotal+=dx*sqrt(1.0+dy*dy); } 
const D step=stotal/n;
D x=0;
D y=0;
const IdxTy leaders=10;
for(IdxTy i=0; i<leaders;  ++i)
{
D x=step*i;
D y=0;
point_type p(x,y);
m_points.push_back(p);

}
// need this crap for equal spacing.. 
for(IdxTy i=leaders; i<n; ++i)
{
D s=step*i;
D der=f/n*cos(f*(i-leaders));
D der2=-f*f/n*sin(f*(i-leaders));
D d1=sqrt(1+ der*der);
D dx=step/d1;
x=x+dx;
D y=sin(f*(x-step*leaders));
point_type p(x,y);
m_points.push_back(p);
} // i 
// TODO fails with zero points. 
const D xlast=m_points.back().x;
const D ylast=m_points.back().y;

for(IdxTy i=0; i<leaders;  ++i)
{
point_type p(xlast+step*(i+1),ylast);
//m_points.push_back(p);

} // i 


} 

void makexxx()
{
const IdxTy n=200;
D f=6.28;
// need approx length of sine wave... 
D stotal=0;
for(IdxTy i=0; i<n; ++i) { D dx=1.0/n; D dy=f*dx*cos(f*i); stotal+=dx*sqrt(1.0+dy*dy); } 
const D step=stotal/n;
D x=0;
D y=0;
for(IdxTy i=0; i<n; ++i)
{
D s=step*i;
D der=f/n*cos(f*i);
D der2=-f*f/n*sin(f*i);
D d1=sqrt(1+ der*der);
//D d2=sqrt(1+ 1.0/der/der);
D dx=step/d1;
x=x+dx;
D y=sin(f*x);
point_type p(x,y);
m_points.push_back(p);

D xn=-der/d1;
D yn=1.0/d1;
point_type pn(xn,yn);
m_normals.push_back(pn);
// the radius of curvature is a scalar
// in the normal direction
// https://en.wikipedia.org/wiki/Radius_of_curvature#
// do not take abs value as sign matters wrt normal
// do 1/r as it can go to infinity at flat locations. 
m_radii.push_back(der2/(d1*d1*d1));

} // i 
} // make 





}; // mjm_etch_curve 



#endif // MJM_ETCH_CURVE_H__ 
