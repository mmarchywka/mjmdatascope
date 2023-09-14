#ifndef MJM_FF_OBJECTS_H__
#define MJM_FF_OBJECTS_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
//#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_pawnoff.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"
#include "mjm_worm_blob.h"
#include "mjm_affine_mesh_hood.h"
#include "mjm_generic_iterators.h"

#include "mjm_canned_methods.h"
//#include "mjm_dscope_dgram.h"

#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"

#include <map> 
#include <set> 
#include <vector> 
#include <deque> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>


// manually coped from mjm_ff_mesh.h 2023-07-08
// Sat 27 May 2023 08:37:15 AM EDT

mjm_global_credits::credit __credit__mjm_ff_objects("mjm_ff_objects"
, "  ");
#if 0
template <class Tr>
class mjm_ff_mesh 
{
 typedef mjm_ff_mesh Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;
#endif

// FIXME doh put this somwhere lol 

namespace {
//static 
int myatoi(const StrTy & s )  { return mjm_canned_methods::myatoi(s.c_str()); }
//static 
int myatoi(const char * c)  { return mjm_canned_methods::myatoi(c); }
typedef unsigned int IdxTy;
bool Bit(const IdxTy f, const IdxTy b)  { return  ((f>>b)&1)!=0; }

}; // namespace

template <class Tr>
class mjm_ff_node
{
// typedef mjm_ff_mesh Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;

public:
mjm_ff_node() : x(0),y(0),b(0) {}
mjm_ff_node(const Line& l ) : x(atof(l[0].c_str())), y(atof(l[1].c_str())), b((l.size()<3)?1:myatoi(l[2])) {
//MM_ERR(MMPR2(x,l[0]))
}
D dist2(const mjm_ff_node & that) const 
	{return (that.x-x)*(that.x-x)+(that.y-y)*(that.y-y); }  
D dist(const mjm_ff_node & that) const {return sqrt(dist2(that)); }  

Os & write(Os & os ) { os << x<<" "<<y<<" "<<b; return os; } 
mjm_ff_node operator-(const mjm_ff_node & that) const
{
mjm_ff_node x=*this;
x.x-=that.x;
x.y-=that.y;

return x; 
}
D  x,y;
IdxTy b;

}; // mjm_ff_node


#if 0
public:
typedef _node vertex_type;
typedef vertex_type V;
typedef std::vector<vertex_type> Verticies;
private:
#endif


template <class Tr, class Myt>
class mjm_ff_triangle
{

 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;



public:
typedef mjm_ff_node<Tr> vertex_type;
typedef vertex_type V;
typedef std::vector<vertex_type> Verticies;

mjm_ff_triangle() : p1(0),p2(0),p3(0),reg(0), m_det(0) {}
mjm_ff_triangle(const Line& l ) 
: p1(myatoi(l[0].c_str()))
, p2(myatoi(l[1].c_str()))
, p3(myatoi(l[2].c_str()))
, reg((l.size()<4)?0:myatoi(l[3])), m_det(0) {}
Os & write(Os & os ) { os << p1<<" "<<p2<<" "<<p3<<" "<<reg; return os; } 
// .5*|a x b|
static D area( const V & p1, const V & p2, const V & p3, const IdxTy flags)
{
D cx= (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x); 
return .5*cx;
} 
// TODO would be better to move offset somewhere.. 
// this shouldn't complte anyway as there is no m_vertcies member.. 
D area(const Myt & m ) const 
	{ return area(m.m_verticies[p1-1],m.m_verticies[p2-1],m.m_verticies[p3-1],0); } 

bool  had_vertex(const IdxTy n) { return (p1==n)||(p2==n)||(p3==n); }
IdxTy shared(const mjm_ff_triangle & y)
{
std::vector<IdxTy> hits(9);
IdxTy i=0;
hits[i]=(p1==y.p1); ++i;
hits[i]=(p1==y.p2); ++i;
hits[i]=(p1==y.p3); ++i;

hits[i]=(p2==y.p1); ++i;
hits[i]=(p2==y.p2); ++i;
hits[i]=(p2==y.p3); ++i;

hits[i]=(p3==y.p1); ++i;
hits[i]=(p3==y.p2); ++i;
hits[i]=(p3==y.p3); ++i;
IdxTy sum =0;
MM_LOOP(ii,hits) { if ((*ii)) ++sum; } 
return sum; 
} // shared


// find the 2 barycentric coods for x,y in this triangle
// give the list of verticies.
//
// TODO container needs the mesh or something to avoid
// confusion and passing redundant crap 
// p3 is 1-a-b, 
template< class Tv> IdxTy coords(D & a,  D & b, const D & x, const D & y, const Tv & verts, const IdxTy flags) const 
{
IdxTy rc=0;
const bool rot1=Bit(flags,0);
// in the only case when this worked, a=-dleta, b=.5 the
// contains_only would have worked well. 
const bool rot2=Bit(flags,1);
 IdxTy _p1=p1;
 IdxTy _p2=p2;
 IdxTy _p3=p3;
if (rot1) {  _p1=p2; _p2=p3; _p3=p1; } 
if (rot2) {  _p1=p3; _p2=p1; _p3=p2; } 
const auto & v1=verts[_p1-1];
const auto & v2=verts[_p2-1];
const auto & v3=verts[_p3-1];
//MyMatrix m(3,3);
 // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
//if (m_det==0)
{
D det=((v2.y-v3.y)*(v1.x-v3.x)+(v3.x-v2.x)*(v1.y-v3.y));
//m_det=1.0/((v2.y-v3.y)*(v1.x-v3.x)+(v3.x-v2.x)*(v1.y-v3.y));
//if (det==0)
// { MM_ERR(" det is zero "<<MMPR(det))}
m_det=1.0/det;
} // m_det
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
a=a*m_det;
b=b*m_det;
return rc;
} // coords

// TODO could be a lot faster only checking until
// bad etc 
// missing zero by 1e-13 sometimes although the standard
// may be picking the wrong one didn't boher to check
// since it don't really matter now. 
// the triangles are oriented so a loop should all point to point
// in same dir with perp interesting line. 
// The trick is the boarder triangles and making sure errors
// are same in competing adjacent ones doh 
template< class Tv> D  contains_only( const D & x, const D & y, const Tv & v1,  const Tv & v2, const IdxTy flags) const 
{
const D xp= v1.x*v2.y-v2.x*v1.y;
const D dx=v2.x-v1.x;
const D dy=v2.y-v1.y;
return xp+y*dx-x*dy;
} // contains_only
template< class Tv> bool  contains_only( const D & x, const D & y, const Tv & verts, const IdxTy flags) const 
{
const auto & v1=verts[p1-1];
const auto & v2=verts[p2-1];
const auto & v3=verts[p3-1];
const D b1=contains_only(x,y,v1,v2,flags);
const D b2=contains_only(x,y,v2,v3,flags);
const D b3=contains_only(x,y,v3,v1,flags);
MM_ERR(MMPR3(b1,b2,b3))
return (b1>=0)&&(b2>=0)&&(b3>=0);
} // contains_only 

template< class Tv> bool  contains( const D & x, const D & y, const Tv & verts, const IdxTy flags) const 
{
const bool print_ab=!true; // Bit(flags,0);
//const bool print_ab= Bit(flags,0);
//const bool print_abo=!true; // Bit(flags,0);
const bool print_abo=!true; // Bit(flags,0);
//const bool print_ab_mid=Bit(flags,1);
D a,b;
coords(a,b,x,y,verts,0);
if (print_ab) {  MM_ERR(MMPR4(a,b,x,y)) }
if (print_ab) { contains_only(x,y,verts,flags); } 

//if (print_ab_mid) if (( a>0)&&(a<1)&&(b>0)&&(b<1)){  MM_ERR(MMPR4(a,b,x,y)<<MMPR((a+b-1.0))) }
//const bool overall= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
//const bool pieces= ((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
//const bool sum= ((a+b)<=1); 
const bool overall1= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
// still needed, not sure anything helped lol 
//if (pieces&&!sum&&((a+b)<1.1))  coords(a,b,x,y,verts,1);
if (!overall1&&((a+b)<1.1)) 
 { coords(a,b,x,y,verts,1);
const bool overall2= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
if (!overall2&&((a+b)<1.1)) coords(a,b,x,y,verts,2);
}
//else return pieces;
const bool overall= ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
if (overall&&(a>.1)&&(a<.9)) {
if (print_abo) {  MM_ERR(MMPR4(a,b,x,y)) }
if (print_abo) { contains_only(x,y,verts,flags); } 

} 
return overall; 
//return ((a+b)<=1) &&((a+b)>=0)&&(a>=0)&&(a<=1)&&(b<=1)&&(b>=0); 
} // contains
template< class Tv> StrTy dump(const Tv & verts, const IdxTy flags) const 
{
Ss ss;
const auto & v1=verts[p1-1];
const auto & v2=verts[p2-1];
const auto & v3=verts[p3-1];
ss<<MMPR4(v1.x,v1.y,v2.x,v2.y)<<MMPR4(v3.x,v3.y,reg,m_det);
return ss.str();
} // dump 

IdxTy p1,p2,p3,reg;
mutable D m_det;
}; // mjm_ff_triangle

template <class Tr> class mjm_ff_edge
{


 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;


public:
mjm_ff_edge() : p1(0),p2(0),b(0) {}
mjm_ff_edge(const Line& l ) 
: p1(myatoi(l[0].c_str()))
, p2(myatoi(l[1].c_str()))
, b((l.size()<3)?0:myatoi(l[2])) {}
Os & write(Os & os ) { os << p1<<" "<<p2<<" "<<b; return os; } 

IdxTy p1,p2,b;
}; // mjm_ff_edge



template <class Tr> class mjm_ff_path
{

 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;





//typedef std::vector<IdxTy> Pts;
typedef std::deque<IdxTy> Pts;
typedef std::map<IdxTy, IdxTy> Index;
public:
mjm_ff_path(): m_id(BAD) {}
void id(const IdxTy n) { m_id=n; } 
// need to add either beginning or end.. 
void add(const IdxTy p, const IdxTy n ) { m_p[p].push_back(n); } 
void padd(const IdxTy p, const IdxTy n ) { m_p[p].push_front(n); } 
IdxTy  add_edge(const IdxTy i1, const IdxTy i2)
{
if (!true)
{
if (m_p.size()==0) m_p.push_back(Pts()); 
m_p[0].push_back(i1);
m_p[0].push_back(i2);
return 0;
}

IdxTy szp=m_p.size();
for(IdxTy p=0; p<szp; ++p)
{
auto& pp=m_p[p];
IdxTy sz=pp.size();
if (sz==1) { MM_ERR(" size is one "<<MMPR3(i1,i2, pp[0])) }
if (sz==0) { MM_ERR(" size is ZERO "<<MMPR2(i1,i2)) }
const IdxTy first=pp.front(); // pp[0];
const IdxTy last=pp.back(); // pp[sz-1];
if (last==i1) {  add(p,i2);  return 0; } 
else if (last==i2) {  add(p,i1);  return 0; } 
else if (first==i2){  padd(p,i1);  return 0; } 
else if (first==i1) {  padd(p,i2); return 0; }  
} // p 
m_p.push_back(Pts());
add(szp,i1); add(szp,i2);

return 0; 
} // add_edge
void concat(Pts& d, const Pts& s, const IdxTy flags)
{
const bool prep=Bit(flags,0);
const bool rev=Bit(flags,1);
if (prep&&!rev) { MM_LOOP(ii,s) {d.push_front(*ii); } }
if (prep&&rev) { MM_SZ_LOOP(i,s,szs) {d.push_front(s[szs-1-i]); } }
if (!prep&&rev) { MM_SZ_LOOP(i,s,szs) {d.push_back(s[szs-1-i]); } }
else { MM_LOOP(ii,s) {d.push_back(*ii); } }

} // concat 


IdxTy concat_all(const IdxTy flags)
{
IdxTy sz=0;
do
{
 sz=m_p.size();
concat(flags);
} while ( sz!=m_p.size()); 
return 0; 
}// concat_all

IdxTy concat(const IdxTy flags)
{
//typedef std::deque<Pts>  PathSeq;
//PathSeq p2;
std::vector<Pts> p2;
IdxTy szp=m_p.size();
for(IdxTy p=0; p<szp; ++p)
{
auto pp=m_p[p];
 // each entry should have size at least two..
if (pp.size()<2) { MM_ERR(" short chain "<<MMPR2(p,pp.size())) }
const IdxTy front=pp[0];
const IdxTy back=pp[pp.size()-1]; 
IdxTy szp2=p2.size();
IdxTy j=0;
for(; j<szp2; ++j)
{
auto pj=p2[j];
if (pj[0]==back){  concat(pj,pp,1); break ; } 
else if (pj[0]==front){  concat(pj,pp,3); break ; } 
else if (pj[pj.size()-1]==front){  concat(pj,pp,0); break ; } 
else if (pj[pj.size()-1]==back){  concat(pj,pp,2); break ; } 
} // j 
if (j==szp2) p2.push_back(pp);
} // p 
MM_ERR(" concat done "<<MMPR2(m_p.size(),p2.size())) 
m_p=p2;
MM_SZ_LOOP(i,m_p,mpsz)
{
MM_ERR(MMPR2(i,m_p[i].size()))
}

return 0;
} // concat 
/*_edge() : p1(0),p2(0),b(0) {}
_edge(const Line& l ) 
: p1(myatoi(l[0].c_str()))
, p2(myatoi(l[1].c_str()))
, b((l.size()<3)?0:myatoi(l[2])) {}
Os & write(Os & os ) { os << p1<<" "<<p2<<" "<<b; return os; } 

IdxTy p1,p2,b;
*/
IdxTy m_id;
std::vector<Pts> m_p;
Index m_idx;

}; // mjm_ff_path





template <class Tr,class Myt> class mjm_ff_dof
{

 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;


typedef mjm_ff_node<Tr> vertex_type;
typedef vertex_type V;
typedef std::vector<vertex_type> Verticies;

public:
typedef D value_type;
typedef std::vector<value_type> Vals;

mjm_ff_dof() : m_mesh(0) {}
const value_type  &  operator[](const IdxTy n) const { return m_vals[n]; }
value_type  &  operator[](const IdxTy n) { return m_vals[n]; }
IdxTy load(const StrTy & fn, const IdxTy flags)
{
std::ifstream is(fn);
return load(is,flags);
}
IdxTy load(IsTy & is, const IdxTy flags)
{
Ragged r;
r.load(is);
return load(r,flags); 
} // load

IdxTy load(const Ragged & r, const IdxTy flags)
{
const bool ff_format=Bit(flags,4);
const bool has_xy=!ff_format&&!Bit(flags,0);
const IdxTy col=has_xy?2:0;
const IdxTy start=(ff_format)?1:0;

const IdxTy sz=r.size();
MM_ERR(" loading dof "<<MMPR4(flags,has_xy,col,sz)<<MMPR(ff_format))
if (sz==0) return 0;
if (ff_format)
{
const IdxTy nodes=atoi(r[0][0].c_str()); 
if (m_mesh) if (nodes!=m_mesh->m_verticies.size())
{
MM_ERR(" sizes wrong in dof "<<MMPR2(nodes,m_mesh->m_verticies.size()))
} // size check 
} // ff_format
for(IdxTy i=start; i<sz; ++i)
{
const Line & l=r[i];
if (l.size()<=col) { MM_ERR(" line too short "<<MMPR2(i,l.size())) } 
else 
{
m_vals.push_back(atof(l[col].c_str())); 
if (has_xy&&m_mesh)
{
const V & p=m_mesh->m_verticies[i];
const D x=atof(l[0].c_str());
const D y=atof(l[1].c_str());
//MM_ERR(MMPR3(x,y,atof(l[col].c_str())));
const D dx=x-p.x;
const D dy=y-p.y;
// TODO wtf is this so  bad 
if(( dx*dx>(1e-5))||( dy*dy>(1e-5)) )
	MM_ERR(" bad coords "<<MMPR4(x,y,p.x,p.y)<<MMPR2(dx,dy))
} // m_mesh 


} // add
} // i 
if (m_mesh)
{
MM_ERR(MMPR3(m_mesh->m_verticies.size(),m_vals.size(),m_name))
}
else MM_ERR(" no mesh specificied ");
return 0;
} // load. 
IdxTy save(const StrTy & fn, const IdxTy flags)
{
std::ofstream os(fn);
return save(os,flags);
} // save 
IdxTy save(OsTy & os, const IdxTy flags)
{
// flags need to be concictent with the mesh flags
// tha also have a precision bit at zero
const bool have_mesh=(m_mesh!=0);
const bool ff_format=Bit(flags,4);
const bool write_coords=!ff_format&&!Bit(flags,1);
const bool write_size=ff_format||Bit(flags,2);
MM_ERR(MMPR4(have_mesh,write_coords,write_size,ff_format)<<MMPR(flags))
IdxTy nprec_coord=12;
IdxTy nprec_v=12;

if (write_size) os<<m_vals.size()<<CRLF;
// Ss ss; ss<<std::setprecision(m_nprec); ss<<x;
if (have_mesh)
{
MM_SZ_LOOP(i,m_vals,sz)
{
const V & p=m_mesh->m_verticies[i];
os<<std::setprecision(nprec_coord); 
if (write_coords) { os<<p.x<<" "<<p.y<<" "; } 
os<<std::setprecision(nprec_v); 
os<<m_vals[i];
os<<CRLF;
} // i 
} // have_mesh 
else
{
os<<std::setprecision(nprec_v); 
MM_SZ_LOOP(i,m_vals,sz)
{
os<<m_vals[i];
os<<CRLF;
} // i 

} // no mesh 

return 0; 
} //save 

IdxTy size() const { return m_vals.size(); } 
void clear()  {  m_vals.clear(); } 
void resize(const IdxTy n ) {  m_vals.resize(n); } 
void add(const D & v) { m_vals.push_back(v); } 
StrTy m_name;
Vals m_vals;
Myt * m_mesh;
}; // mjm_ff_dof


template <class Tr> 
class mjm_ff_seg
{
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;

public:
mjm_ff_seg() 
	: x(0),y(0),a(-1),b(0),dsx(0),dsy(0),p1(0),p2(0),t(0),seg_flag(0) {}
//_seg(const Line& l ) : x(atof(l[0].c_str())), y(atof(l[1].c_str())), b((l.size()<3)?1:myatoi(l[2])) {
//MM_ERR(MMPR2(x,l[0]))
//}
template<class Tp>
//mjm_ff_seg(const mjm_ff_node & p) 
mjm_ff_seg(const Tp & p) 
	: x(p.x),y(p.y),a(-1),b(p.b),dsx(0),dsy(0),p1(0),p2(0),t(0),seg_flag(0) {}
void set(const IdxTy _p1, const IdxTy _p2, const IdxTy _t)
{ p1=_p1; p2=_p2; t=_t; } 
D dist2(const mjm_ff_seg & that) const 
	{return (that.x-x)*(that.x-x)+(that.y-y)*(that.y-y); }  
D dist(const mjm_ff_seg & that) const {return sqrt(dist2(that)); }  
bool terminal() const { return Bit(seg_flag,0); } 
mjm_ff_seg& terminal(const bool x ) 
{ if (x) seg_flag|=1; else seg_flag=~1; return *this; } 
bool closed() const { return Bit(seg_flag,1); } 
mjm_ff_seg& closed(const bool x ) 
{ if (x) seg_flag|=2; else seg_flag=~2; return *this; } 
StrTy dump(const IdxTy flags=0) const 
{ Ss ss; 
const bool xyb=Bit(flags,0);
if (xyb) { write(ss);  } 
else { ss<<" "<<MMPR4(p1,p2,t,a)<<MMPR(seg_flag); } 
return ss.str(); 
} 
Os & write(Os & os )const { os << x<<" "<<y<<" "<<b; return os; } 
void ds(const mjm_ff_seg & that)
{
dsx=x-that.x;
dsy=y-that.y;
}
// location and fraction from p1, gradient perp to ds ;
D  x,y,a,grad;
IdxTy b; // unknown for now. 
D dsx,dsy; // direction of arc if known 
// if just one point, second one is zero 
// t is ZERO based so always valid... doh? 

// now v(p1)<v(p2) 
IdxTy p1,p2,t; // point 1, point 2, and triangle  
IdxTy seg_flag;

}; // mjm_ff_seg

// the segs are elements/
template <class Tr,class Myt> 
class mjm_ff_trail : public std::vector<mjm_ff_seg<Tr> >
{
typedef mjm_ff_trail<Tr,Myt> TrailTy;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;


typedef mjm_ff_seg<Tr> Seg;
typedef std::vector<Seg> Super;
typedef std::vector<IdxTy> IdxVal;
typedef std::set<IdxTy> SidxVal;
typedef std::map<IdxTy, IdxVal > Index;
typedef std::map<IdxTy, IdxTy > IdxMap;
typedef std::map<IdxTy, SidxVal > Sindex;
typedef std::map<IdxTy, IdxTy > HistoMap;
public:
typedef mjm_ff_node<Tr> _node;
typedef mjm_ff_seg<Tr> _seg;
typedef typename Super::value_type value_type;
template <class Tc> class xthenyf
{public:
mjm_ff_trail * m;
Tc  * p;
#define XYF(a) (*m)[(*p)[a].front()] 
bool operator()(const IdxTy a, const IdxTy b) 
	{ const bool eq= ( XYF(a).x==XYF(b).x)  ;
	 if (eq) return XYF(a).y<XYF(b).y ;
	return  (XYF(a).x<XYF(b).x );
}
}; // xthenyf
#undef XYF
template <class Tc> class xtheny
{public:
//mjm_ff_trail * p;
Tc  * p;
bool operator()(const IdxTy a, const IdxTy b) 
	{ const bool eq= ((*p)[a].x==(*p)[b].x );
	 if (eq) return ((*p)[a].y<(*p)[b].y );
	return  ((*p)[a].x<(*p)[b].x );
}
};

IdxTy size() const { 
if (Super::size()!=m_trail_flags.size())
	MM_ERR(" sizes do not match "<<MMPR2(Super::size(),m_trail_flags.size()))
return Super::size(); } 
IdxTy trails() const { return m_trail_heads.size(); } 
void push_back(const _node & ni, const IdxTy fl=0) 
	{ push_back(_seg(ni)); m_trail_flags.push_back(fl); } 
void push_back(const _seg & n, const IdxTy fl=0 ) 
	{ Super::push_back(n); m_trail_flags.push_back(fl);  } 
void clear() { Super::clear(); m_trail_flags.clear(); m_trail_heads.clear(); } 
void head(const IdxTy tr, const IdxTy loc ) { m_trail_heads[tr]=loc;}

void new_trail() { m_trail_heads[m_trail_heads.size()]=size(); } 
IdxTy sorted(mjm_ff_trail & y, const IdxTy flags)
{
if (size()==0) return 0;
const D tol=1e-9;
// TODO doubled nodes but ordered later. 
const bool remove_ident=!false; // Bit(flags,0);
std::vector<IdxTy> order(size());
for(IdxTy i=0; i<size(); ++i) order[i]=i;
xtheny<mjm_ff_trail>  xty; xty.p=this;
std::sort(order.begin(),order.end(),xty);
y.resize(size());
for(IdxTy i=0; i<size(); ++i) y[i]=(*this)[order[i]];

//y=*this;

if (remove_ident)
{
 mjm_ff_trail z;
z.push_back(y[0]);
for(IdxTy i=1; i<size(); ++i) 
{
const auto & p1=y[i];
const auto & p2=z.back(); 
const D dx=p1.x-p2.x;
const D dy=p1.y-p2.y;

D sx=p1.x+p2.x;
D sy=p1.y+p2.y;
if (sx==0) sx=1;
if (sy==0) sy=1;
if ( fabs(dx/sx)<tol) if ( fabs(dy/sy)<tol) continue;
//if ( y[i].dist(z.back())<tol) continue;
z.push_back(y[i]);
}
y=z;
} // remove_ident
index();
// intereseting but not useful and slo
//adjacents();
// this needs to operate on y not "this" 
y.m_mesh=m_mesh;
y.match_seg_vertex(y);
//match_seg_vertex_fck(y);
return 0; 
} // sorted
void index(const IdxTy flags=0)
{
// add for constraint...
clear_indexes();
MM_SZ_LOOP(i,(*this),sz)
{
const Seg & s=(*this)[i];
m_st[s.t].push_back(i);
// index zero vertex entries too, note that tri is zero based
m_sv[s.p1].push_back(i);
m_sv[s.p2].push_back(i);
} // i 
MM_LOOP(ii,m_st) { ++m_sth[(*ii).second.size()]; }
MM_LOOP(ii,m_sv) { ++m_svh[(*ii).second.size()]; }

MM_LOOP(ii,m_sth) { MM_ERR(" sth "<<MMPR2((*ii).first,(*ii).second)) }
MM_LOOP(ii,m_svh) { MM_ERR(" svh "<<MMPR2((*ii).first,(*ii).second)) }
MM_ERR(" exacts "<<MMPR(m_svh[0]))

} // index 
// try to line up with seg numbers...
typedef std::deque<IdxTy>  Partial;
typedef std::vector<Partial> Partials;

// index is vertex number, 1 based, value is triangles sharing that vertex
void index_tri(Sindex& idx) const 
{
MM_SZ_LOOP(i,m_mesh->m_tri,sz)
{
const auto & t=m_mesh->m_tri[i];
idx[t.p1].insert(i+1);
idx[t.p2].insert(i+1);
idx[t.p3].insert(i+1);

} // i 

} // index_tri
// y was reordered so indexing "this" only workds if
// called on instance y... 
void index_edges(Index & edges) const 
{ 
MM_SZ_LOOP(i,(*this),sz)
{ 
const Seg& si=(*this)[i];
// could get all twice.. 
edges[si.p1].push_back(i);
edges[si.p2].push_back(i);
// now const 
//si.b=0; // this will be used for a "positioned" marker 
} // i 

} // index_edges
void coalesce(mjm_ff_trail & y, Partials & frags, Sindex & tri) const 
{
Partials fnew,fnewer;
condense_frags(fnew,frags,tri);
// sort frags so combining them is not a mess
// although there should be very few of them
// if the path selection worked but now there are 
// hundreds...
std::vector<IdxTy> order(fnew.size());
for(IdxTy i=0; i<fnew.size(); ++i) order[i]=i;
//xthenyf<_trail>  xty; xty.m=this;  xty.p=&fnew;
xthenyf<Partials>  xty; xty.m=&y;  xty.p=&fnew;
std::sort(order.begin(),order.end(),xty);
//fnewer.resize(size());
for(IdxTy i=0; i<order.size(); ++i) fnewer.push_back(fnew[order[i]]);
fnew=fnewer;

//fnew=frags;
const IdxTy fragssz=frags.size();
const IdxTy fnewsz=fnew.size();
//const IdxTy sz=(*this).size();
MM_ERR(MMPR3(fragssz,fnewsz,size()))
mjm_ff_trail tn; /// (sz);
MM_SZ_LOOP(i,fnew,xxx){ 
tn.new_trail(); 
//const bool closed=(*ii).back().closed();
//MM_SZ_LOOP(j,fnew[i],yyy) {tn.push_back((*this)[fnew[i][j]],closed?1:0); }  } 

// keep gradient pointing to right as curve traversed...
// orient
bool fwd=true;
const auto & fni=fnew[i];
if (fni.size()>1)
{
const auto & s1=(*this)[fni[0]];
const auto & s2=(*this)[fni[1]];
// p1 dof smaller than p2, 
const auto & p11=m_mesh->m_verticies[s1.p1-1];
const auto & p12=m_mesh->m_verticies[s1.p2-1];
const auto & p21=m_mesh->m_verticies[s2.p1-1];
const auto & p22=m_mesh->m_verticies[s2.p2-1];
const auto nup1=p12-p11;
const auto nup2=p22-p21;
const D dx=s2.x-s1.x;
const D dy=s2.y-s1.y;
const D ns1=nup1.x*dy-nup1.y*dx;
const D ns2=nup2.x*dy-nup2.y*dx;
MM_ERR(MMPR2(ns1,ns2))

if (ns1<0) fwd=false;


} // fni.size
if (fwd)
{ MM_SZ_LOOP(j,fnew[i],yyy) {tn.push_back((*this)[fnew[i][j]]); }  } 
else
{ MM_SZ_LOOP(j,fnew[i],yyy) {tn.push_back((*this)[fnew[i][yyy-1-j]]); }  } }


MM_LOOP(ii,fnew)
{
const bool closed=circular((*ii),tri);
MM_ERR(MMPR4((*ii).size(),(*ii).front(),(*ii).back(),closed))
const IdxTy fr=(*ii).front();
const IdxTy ba=(*ii).back();
const IdxTy term=(ba>fr)?ba:fr;
// TODO FIXME wth?
//if (closed) (*ii).front().closed(true; 
if (closed) tn[term].closed(true); 
tn[term].terminal(true); 
} // check circles 

//(*this)=tn;
y=tn; 
} // coalesce

SidxVal intersect(const SidxVal &s1,const SidxVal &s2,const SidxVal &s3) const 
{
SidxVal x;
MM_LOOP(ii,s1)
{
const IdxTy c2=s2.count(*ii);
if (c2==0) continue;
const IdxTy c3=s3.count(*ii);
if (c3==0) continue;
x.insert(*ii);
} 
return x; 
} // intersect 

// TODO indexing with cpoint will add zed length entries
// if not already there no const
IdxTy align( Sindex & tri, const Seg & s1, const Seg & s2,const IdxTy flags=0) const 
{
const bool disjoint=Bit(flags,0);
const IdxTy p11=s1.p1;
const IdxTy p12=s1.p2;
const IdxTy p21=s2.p1;
const IdxTy p22=s2.p2;
IdxTy cpoint,p1,p2;
if (p11==p21) {  cpoint=p11; p1=p12; p2=p22; } 
else if (p11==p22){  cpoint=p11; p1=p12; p2=p21; }
else if (p12==p21){  cpoint=p12; p1=p11; p2=p22; } 
else if (p12==p22){  cpoint=p12; p1=p11; p2=p21; } 
else {
if (!disjoint) { MM_ERR(" no common point "<<MMPR4(p11,p12,p21,p22)) } 
return 0; } 
// find intersection of 3 sets, tri[cpoint], tri[p1], tri[p2]
SidxVal sis=intersect(tri[cpoint],tri[p1],tri[p2]);
// number of triangles shared by s1 and s2
return sis.size(); 
} // align 

////////////////////////////////////////////////////////
// this is called  on instance  y


void match_seg_vertex(mjm_ff_trail & y, const IdxTy flags=0) const 
{
if (size()==0) { MM_ERR(" nothing to do ")  return ; } 
if (size()==0) { MM_ERR(" nothing to do ")  return ; } 
Sindex tri; // set of all triangles  contaiint key pt
//Index edges; // key=vertex, val= vector of segs having that vertex
Partials frags; // vector of partially assembled contigs 
std::set<IdxTy> vadded; // set of indexes already added ( zero base) 
index_tri(tri); 
//index_edges(edges);
//IdxTy iszmax=0;
//MM_LOOP(ii,edges) { IdxTy n=(*ii).second.size(); if (n>iszmax) iszmax=n; } 

// add any leftovers as singletons that should glue things
for(IdxTy i=0; i<size(); ++i) // y.size?
{
if (vadded.count(i)) { continue; } 
Partial p;
p.push_back(i);
vadded.insert(i);
add_frag(frags,p,tri);

} // i 

coalesce(y,frags,tri); // finally join any separate units
} // match_seg_vertex




// return 1 for front 2 for back 
IdxTy see_add( const Seg & s ,const Seg & pf, const Seg & pb, Sindex & tri) const 
{
//IdxTy align( Sindex & tri, const Seg & s1, const Seg & s2) const 
const IdxTy hitsf=align(tri,s,pf,1);
const IdxTy hitsb=align(tri,s,pb,1);
const IdxTy thhit=0;
if (hitsf>thhit)
{
if (( pf.p1!=0) &&((pf.p1==s.p1)||(pf.p1==s.p2))) { return 1; }
if (( pf.p2!=0) &&((pf.p2==s.p1)||(pf.p2==s.p2))) { return 1; }
}
if (hitsb>thhit)
{
if (( pb.p1!=0) &&((pb.p1==s.p1)||(pb.p1==s.p2))) { return 2; }
if (( pb.p2!=0) &&((pb.p2==s.p1)||(pb.p2==s.p2))) { return 2; }
}
return 0;
} // see_Add




bool circular( Partial &g, Sindex & tri) const
{
if (g.size()==0) return false  ;
const Seg & sf=(*this)[g.front()-0];
const Seg & sb=(*this)[g.back()-0];
IdxTy ac=see_add(sb,sf,sf,tri);
return (ac!=0);
} // circular 

// either push g onto f or append to a mating entry in f
void add_frag( Partials & f,  Partial &g, Sindex & tri) const
{
//const bool only_check_ends=!false; // no obvious effect
if (g.size()==0) return ;
const Seg & sf=(*this)[g.front()-0];
const Seg & sb=(*this)[g.back()-0];
MM_SZ_LOOP(i,f,sz)
{
Partial & p=f[i];
const Seg & pf=(*this)[p.front()-0];
const Seg & pb=(*this)[p.back()-0];
// if the back of g mates to fron, they should concat
// otherwise need to revers g. 
IdxTy ac=see_add(sb,pf,pb,tri);
if (ac==1)  concatf(p,g,0); // p.push_front(segn);
else if (ac==2) concatf(p,g,1); //  p.push_back(segn);
if (ac!=0) return; 
// revers of above 
ac=see_add(sf,pf,pb,tri);
if (ac==1) concatf(p,g,2); //  p.push_front(segn);
else if (ac==2) concatf(p,g,3); // p.push_back(segn);
if (ac!=0) return ;
} // i 
f.push_back(g);
} // add_frag 3 param 



void concatf(Partial & d, const Partial & s, const IdxTy flags) const
{
//IdxTy szs=s.size();
//MM_ERR(MMPR2(__FUNCTION__,flags))
switch (flags)
{
//put s in front of d
case 0:{MM_SZ_LOOP(i,s,sz) { d.push_front(s[sz-1-i]); }     break; } 
//i put reversed s behind d
case 1:{MM_SZ_LOOP(i,s,sz) { d.push_back(s[sz-1-i]);  }    break; } 

//put reversed s in front of d
case 2:{MM_SZ_LOOP(i,s,sz) { d.push_front(s[i]);   }   break; } 
//i puts behind d
case 3:{MM_SZ_LOOP(i,s,sz) { d.push_back(s[i]);    }   break; } 

default:
MM_ERR(" bad code "<<MMPR(flags))

} // switch 


} // concatf
void condense_frags( Partials & f,  const Partials &fin, Sindex & tri ) const
{
//f=fin;
MM_ERR( " condensing, should be done incrementally lol "<<MMPR(fin.size()))
Partials fins=fin;
IdxTy iter=0;
f=fin;
while (true)
{
f.clear();
const IdxTy insz=fins.size();
MM_LOOP(ii,fins)
{
add_frag(f,(*ii),tri);
} // ii 
const IdxTy outsz=f.size();
MM_ERR(MMPR3(iter,insz,outsz))
//if (iter>10) break; 
if (insz==outsz) break; 
if (1==outsz) break; 
++iter;
fins=f;
} // true 

/*
MM_LOOP(ii,f)
{
const bool closed=circular((*ii),tri);
//if (closed) (*ii).front().closed(true; 
if (closed) (*this)[(*ii).back()].closed(true); 
(*this)[(*ii).back()].terminal(true); 
} // check circles 

*/

if (false) { 
MM_LOOP(ii,f)
{
Ss ss;
MM_LOOP(jj,(*ii))
{
ss<<" "<<(*jj); 
} // jj 
MM_ERR(" seq "<<MMPR(ss.str()))

} // ii 
} // false 
} // condense_frags

// put an oriented version into oriented 
// traversal should insure each grad points to the right
// or lower side is on the lefth 
void orient(const mjm_ff_trail & oriented, const IdxTy flags) const
{
// OBSOLETE se "fwd" and "fni" in the sorting method. 

} // orient 

IdxTy constrain_dists(const D & mn, const D & mx, const IdxTy flags)
{
IdxTy rc=0;
const D d2min=mn*mn;
const IdxTy trails=m_trail_heads.size();
if (trails==0) if (size()!=0)
{ MM_ERR(" no trail heads "<<MMPR2(trails,size())) }
IdxTy ndel=0;
TrailTy  & trail=*this;
MM_LOOP(ii,m_trail_heads)
{
IdxTy ntraildel=ndel;
const IdxTy ntrail=(*ii).first;
const IdxTy first=(*ii).second;
(*ii).second-=ndel;
auto jj=ii; ++jj;
const IdxTy end=(jj==m_trail_heads.end())?size():(*jj).second;
IdxTy jzed=1;
for(IdxTy j=(first+1); j<end; ++j)
{
// if segs are too close, remove and note to update indexe
D d2=trail[j].dist2(trail[j-jzed]);
MM_ERR(MMPR4(d2,d2min,mn,j-jzed))
if (d2<d2min) { ++ndel; ++jzed;  continue; }
if (ndel){  trail[j-ndel]=trail[j]; jzed=1; }

} // j 
// if this is closed, check first and last...
ntraildel=ndel-ntraildel;
MM_ERR(MMPR3(ntrail,ntraildel,ndel))

} // ii 
if (ndel) trail.resize(trail.size()-ndel); 
trail.index();
return rc;
} // constrain_dists
void clear_indexes()
{
m_st.clear();
m_sv.clear();
m_t1.clear();
m_t2.clear();
m_t1h.clear();
m_t2h.clear();

} // clear_indexes
// _trail derives from vector and it is a vector of segs.
// this may be confusing... 
Index m_st,m_sv,m_t1,m_t2;
HistoMap m_sth, m_svh;
HistoMap m_t1h, m_t2h;
IdxVal m_trail_flags;
IdxMap m_trail_heads; // starting location of each disjoint seg. 
Myt * m_mesh;
}; // _trail 

// functor or sparese matrix to load dof's from another
// mesh 
template <class Tr,class Myt>
class mjm_adapt_adopt_matrix
{


 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
typedef OsTy Os;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
// BAD should be in traits
enum { BAD=~0}; 
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;



typedef std::pair<IdxTy,D> Me;
typedef std::vector<Me>  Col;
typedef std::map<IdxTy, Col> Rows;
public:
mjm_adapt_adopt_matrix(): m_d_sz(0) {}
void resize(const IdxTy n) { m_d_sz=n; } 
void add(const IdxTy i, const IdxTy j, const D & v) 
{ if (v!=0) m_sparse[i].push_back(Me(j,v)); } 
typedef mjm_ff_dof<Tr,Myt> _dof;
_dof operator()(const _dof& s) const
{
_dof d;
d.resize(m_d_sz);
MM_SZ_LOOP(i, d.m_vals,szd)
{
const auto ii=m_sparse.find(i);
if (ii==m_sparse.end()) continue;
const auto & row=(*ii).second; 
//const auto & row=m_sparse[i];
D v=0;
MM_LOOP(jj,row) { v+=(*jj).second*s[(*jj).first]; }
d.m_vals[i]=v; 
} // ii 
return d;
} // operator()
StrTy dump() const
{
Ss ss;
IdxTy m1=0, m2=0;
MM_LOOP(ii,m_sparse)
{
if ((*ii).first>m1) m1=(*ii).first; 
MM_LOOP(jj,(*ii).second)
{
if ((*jj).first>m2) m2=(*jj).first;
} // jj 
} // ii 
ss<<MMPR4(m1,m2,m_d_sz,m_sparse.size());
return ss.str();
} // dump
IdxTy m_d_sz;
Rows m_sparse;

}; // mjm_adapt_adopt_matrix

#endif // MJM_FF_OBJECTS_H__ 
