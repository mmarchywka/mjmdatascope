#ifndef MJM_AFFINE_MESH_HOOD_H__
#define MJM_AFFINE_MESH_HOOD_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"
#include "mjm_canned_methods.h"
//#include "mjm_cli_ui.h"
#include "mjm_tokenized_collections.h"
#include "mjm_worm_blob.h"
#include "mjm_geo_contains.h"


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

// Fri 16 Jun 2023 07:32:28 AM EDT
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_affine_mesh_hood   
// QUICKCOMPILE  g++  -Wall -Wno-misleading-indentation  -std=gnu++11 -DTEST_MJM_AFFINE_MESH_HOOD -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_affine_mesh_hood.h  -o mjm_affine_mesh_hood.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_affine_mesh_hood("mjm_affine_mesh_hood" , "  ");

template <class Tr>
class mjm_affine_mesh_hood 
{
 typedef mjm_affine_mesh_hood Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

enum {BAD=~0 }; 
// TYPEDEF
typedef mjm_canned_methods Canned;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef std::vector<StrTy> Words;
typedef mjm_string_base_params<Tr> BaseParams;

typedef mjm_worm_blob<Tr,D> WORMBlob;
typedef mjm_geo_contains<Tr> GeoUtil;
// vector of points inthis lot 
//typedef std::pair<D,D> Point;

// hood and serial for a point
// NB serial can change if points are deleted.. 
typedef std::pair<IdxTy,IdxTy> HoodPt;
typedef std::set<HoodPt> HoodPtSet;


class _guard
{
public:
//const D & left, const D & right, const D& bot, const D& top
_guard() { Init(0,0,0,0,0); }
_guard( const D & left, const D & right, const D& bot , const D& top , const D& gf )
{ Init(left,right,bot,top,gf); }
void Init( const D & left, const D & right, const D& bot , const D& top , const D& gf )
{
ex=right-left;
ey=top-bot;
D guardx=gf*ex; 
D guardy=gf*ey; 
gl=left+guardx;
gr=right-guardx;
gt=top-guardy;
gb=bot+guardy;
} // Init

bool guarded(const D &x, const D & y ) const
{
const bool  iguarded=(x>gr)||(x<gl)||(y>gt)||(y<gb);
return iguarded;
}
D ex,ey,gl,gr,gb,gt;

}; // _guard
typedef _guard GuardType;


class _Point { public: 

_Point(): x(0),y(0),src(0) {}
_Point(const D & _x, const D & _y, const IdxTy _src)
	: x(_x),y(_y),src(_src) {}
template <class Tp> D dist2(const Tp & that) const 
//{ return x*x+that.x*that.x -2.0*x*that.x
{ return (x-that.x)*(x-that.x)+ (y-that.y)*(y-that.y); } 

D x,y; IdxTy src;

 }; // _Point 
typedef _Point Point;
class _hood_desc {
typedef std::vector<Point> Points;
public:
void add(const D & x, const D & y, const IdxTy src)
{ m_p.push_back(Point(x,y,src)); } // add
void add(const Point & p ) { m_p.push_back(p); } 

IdxTy size() const { return m_p.size(); } 
Points m_p;

}; // _hood_desc
typedef _hood_desc HoodDesc;
typedef std::map<IdxTy, HoodDesc> LotMap;
typedef typename LotMap::const_iterator Lii;
typedef typename LotMap::iterator Liinc;
typedef std::set<IdxTy> IdxSet;
// added for boxes, may be better with dense grid?
// TODO move to new class  
typedef std::map<IdxTy, IdxTy> BoxMap;

class _rule {
public:
_rule() : a(0), action(0) {}
const bool applies( const IdxTy tsrc) const 
{ const bool has=(src.count(tsrc)!=0);  
 return (has&&!exclude_set()) ||( !has&&exclude_set());  }
const bool exclude_set() const 
	{ return Bit(action, EXC_SET); }  
// the src set members are excluded from testing 
_rule & exclude_set(const bool x ) { twiddle(EXC_SET,x); return *this; }
// in case of proximity alert, remove the existing point
_rule & remove_hit(const bool x ) { twiddle(REMOVE_HIT,x); return *this; }
// only remove near collisions do not place probe point 
_rule & clean(const bool x ) { twiddle(CLEAN,x); return *this; }
void twiddle(const IdxTy n, const bool x) 
{ if (x) action|=(1<<n); else action &= ~(1<<n); }    
_rule & add(const IdxTy n) { src.insert(n); return *this;  }
_rule & clear_action() { action=0; return *this; } 
_rule &clear_src() { src.clear(); return *this; } 
const bool remove_hit() const 
	{ return Bit(action, REMOVE_HIT)|| Bit(action,CLEAN); }  
const bool place() const 
	{ return !Bit(action,CLEAN) ; }  
IdxSet src;
D a;
// EXC_SET : if set, the src set is those to exclude from testing 
// REMOVE_HIT : if set, remove the existing offender(s) and place this
// CLEAN : remove all offending pounts but do not place a new one.  
enum { EXC_SET=0, REMOVE_HIT=1, CLEAN=2  }; 
IdxTy action;
}; // _rule


class _delete_que
{

typedef std::set<IdxTy> DelSet;
typedef std::map<IdxTy, DelSet> DelMap;
public:
void add(const IdxTy hood, const IdxTy el) { m_del[hood].insert(el); } 
void clear() { m_del.clear(); } 
IdxTy size() const  { return  m_del.size(); } 
template <class Tv> IdxTy remove(Tv & x, const IdxTy hood) const
{
IdxTy n=0; 
const auto ii=m_del.find(hood);
if (ii==m_del.end()) { return 0; } 
const DelSet & s=(*ii).second;
IdxTy szy=x.size()-s.size();
//Tv y(szy);
IdxTy j=0;
MM_SZ_LOOP(i,x,xsz) 
	{ if ( !s.count(i)) { x[i-n]=x[i]; ++j; } else ++n; } // i 
x.resize(szy);
//x=y;
return n;
} // remove 

template <class Tm> IdxTy remove(Tm & m ) const 
{
IdxTy n=0;
MM_LOOP(ii,m_del)
{
const IdxTy hood=(*ii).first;
n+=remove(m[hood].m_p,hood);
} // ii 
return n;
} // remove 

DelMap m_del;

}; // _delete_que





// API

public:

typedef _delete_que delete_que_type;
typedef _rule rule_type;
typedef std::vector<rule_type> rule_list_type;
typedef WORMBlob points_io_type;

mjm_affine_mesh_hood() {Init();}
mjm_affine_mesh_hood(const StrTy & s, const IdxTy flags) {Init(s,flags);}
~mjm_affine_mesh_hood() {}
void clear() { m_map.clear(); m_pts=0; } 
IdxTy size() const { return m_pts;}
// delete empty hoods.. 
IdxTy hoods() const { return m_map.size(); } 
IdxTy hood(const D & x, const D & y ) const 
//{ return int((x-m_left)/m_dx) + int((y-m_bottom)/m_dy)*m_nx; }
{ return hoodx(x) + hoody(y)*m_nx; }
IdxTy nhood(const int n) const { const auto ii=m_map.find(n);
if (ii==m_map.end()) return 0; return (*ii).second.size(); } 
int  hoodxy(const int  nx, const int ny ) const { return nx + ny*m_nx; }
int hoodx(const D & x) const { return int((x-m_left)/D(m_dx)); }
int hoody(const D & y) const { return int((y-m_bottom)/D(m_dy)); }
D hoodxmin(int hx ) const { return hx*m_dx+m_left;  } // xmin for hood hx
D hoodymin(int hy ) const { return hy*m_dy+m_bottom;  } // xmin for hood hx
const D & left() const { return m_left; } 
const D & right() const { return m_right; } 
const D & bottom() const { return m_bottom; } 
const D & top() const { return m_top; } 
const D & height() const { return m_top-m_bottom; } 
const D & width() const { return m_right-m_left; } 
// from surface.idp, 
// left=1, top=2, right=3, boyt=4, surf=5
int edge(const D & x, const D & y, const IdxTy flags) const
{
if (x==m_left) return 1;
if (x==m_right) return 3;
if (y==m_top) return 2;
if (y==m_bottom) return 4;
return 0; 
} // edge
int edge_codes(const D & x, const D & y, const IdxTy flags) const
{
int rc=0;
if (x==m_left) rc|=(1<<1);
if (x==m_right) rc|=(1<<3);
if (y==m_top) rc|=(1<<2);
if (y==m_bottom) rc|=(1<<4);
return rc; 
} // edge_codes



void add(const D & x, const D & y, const IdxTy src )  
{ m_map[hood(x,y)].add(x,y,src); ++m_pts; } 

// set of all hooods containing line segment 
template <class Ty> IdxTy  isin(Ty & set, const D & x1,const D & y1, const D & x2, const D & y2, const IdxTy flags) const
{ return IsIn(set,x1,y1,x2,y2,flags); }

template <class Ty> IdxTy  fill(Ty & set,  const IdxTy flags=0) const
{ return Fill(set,flags); } 
template <class Ty> IdxTy  infer(const Ty & ps, const IdxTy flags)
 { return Infer( ps, flags); } 
template <class Ty> IdxTy  add(const Ty & ps, const IdxTy src, const IdxTy flags)
 { return Add( ps, src, flags); } 
template <class Ty> IdxTy  add(const Ty & ps, const IdxTy src, const rule_list_type & rl, const IdxTy flags)
 { return Add( ps, src, rl,flags); } 

//
template <class Ty> IdxTy  
	add_priority(const Ty & ps, const D & distbackground,const D & distprioriti, const IdxTy src, const IdxTy flags)
 { 
const bool over_all=Bit(flags,0); // otherwise just mask zero src 
const bool over_zero_space=Bit(flags,1); // otherwise just mask zero src 
rule_list_type rl;
rl.push_back(rule_type());
MM_ERR(MMPR4(distbackground,distprioriti, over_all,over_zero_space))
// distance is for first rule to clear out background pounts, src 0  
rl.back().a=distbackground;
if (over_zero_space)
{ // use first rule to nearby point from bg
 rl.back().exclude_set(false).clean(true).add(0);
rl.push_back(rule_type());
// this distance is between priority points 
rl.back().a= distprioriti; // .1*dist;
// now place if other priority points are not there.  
 rl.back().exclude_set(true).remove_hit(false).add(0);
}else 
// if over_all, remove all near collisions 
if (over_all) rl.back().exclude_set(true).remove_hit(true);
// oterwise only remove those in src zero the background points. 
else rl.back().exclude_set(false).remove_hit(true).add(0);
return Add( ps, src, rl,flags); } 

template <class Tc> void add_priority_curves(const Tc &cvs,const D & space,const D & gspace,const IdxTy base,const IdxTy flags )
{
//MM_DIE(MMPR2(space,gspace))
MM_SZ_LOOP(i, cvs,cvsz)
{
MM_SZ_LOOP(j,cvs[i],cvssz)
{
// d,src,flags
// if too close, remove background points but if too 
// close to existing curve just drop this point. 
// space is to background, second distance is between iso crve pts. 
//mh.add_priority(cvs[i][j].m_points,space,.1*space,i+1,2); // src and flags 
add_priority(cvs[i][j].m_points,space,gspace,base+i,flags); // src and flags 
} // j 
//break;
} // i 
} // add_priotiyy_curves



IdxTy delete_que( delete_que_type & dq, const IdxTy flags){ return  Delete(  dq,  flags); } 
// these need to be oriented but not sure which way lol. 
template <class Tp > IdxTy delete_tri(delete_que_type & dq, const Tp & p1, const Tp & p2, const Tp & p3 , const IdxTy flags)
 { return DeleteTri( dq,  p1, p2, p3 , flags); } 

IdxTy save_points(const StrTy & fn, const IdxTy flags)
{ return SavePoints(fn,flags); } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
// ADD-ons for boxes for layout. 
// Add a src point to all hoods that penetrate box but not on edge
void add_box(const D&x,const D&y, const D&w, const D&h, const IdxTy src)
{ AddBox(x,y,w,h,src); } 
template <class Tb> void add_box(const Tb&box, const IdxTy src)
{ AddBox(box.x,box.y,box.w,box.h,src); } 
// return true if any hoods pentrating box have points
IdxTy test_box(const D&x,const D&y, const D&w, const D&h)
{ return TestBox(x,y,w,h); } 
IdxTy test_box(const D&x,const D&y)
{ return TestBox(x,y); } 


private:
static bool Bit(const IdxTy f, const IdxTy b) { return  ((f>>b)&1)!=0; }
// should loop over map now 
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
void AddBox(const D&x,const D&y, const D&w, const D&h, const IdxTy src)
{
if (false)
{
static int cnt=0;
if (cnt<20) {  MM_ERR(MMPR4(x,y,w,h)<<MMPR3(m_dx,m_dy,src)) } 
++cnt;
}
D x0=x;
if (x0<m_left) x0=m_left;
D xf=x+w;
if (xf>m_right) xf=m_right;
D y0=y;
if (y0<m_bottom) y0=m_bottom;
D yf=y+h;
if (yf>m_top) yf=m_top;


for(D xx=x0; xx<(xf); xx+=m_dx)
{for(D yy=y0; yy<(yf); yy+=m_dy)
{
m_boxes[hood(xx,yy)]=src;
}}

} // AddBox
// should include a fail fast flag for bool return on first collision
IdxTy TestBox(const D&x,const D&y, const D&w, const D&h)
{
D x0=x;
if (x0<m_left) x0=m_left;
D xf=x+w;
if (xf>m_right) xf=m_right;
D y0=y;
if (y0<m_bottom) y0=m_bottom;
D yf=y+h;
if (yf>m_top) yf=m_top;

for(D xx=x0; xx<(xf); xx+=m_dx)
{for(D yy=y0; yy<(yf); yy+=m_dy)
{
const IdxTy code=TextBox(xx,yy);
if (code!=BAD) return code;
}}

return BAD;
} // TestBox

IdxTy TestBox(const D&x,const D&y)
{
const auto ii =m_boxes.find(hood(x,y));
if (ii==m_boxes.end()) return BAD;
return (*ii).second; 

} // TextBox

template <class Ty> 
IdxTy Add(const Ty & ps, const IdxTy src, const IdxTy flags)
{
MM_LOOP(ii,ps)
{
const auto & p=(*ii);
const D & x=p.x;
const D & y=p.y;
add(x,y,src); 
} // ii 
return 0; 
} // Add

template <class Tp > IdxTy DeleteTri(delete_que_type & dq, const Tp & p1, const Tp & p2, const Tp & p3 , const IdxTy flags)
{
GeoUtil gu;
HoodPtSet  hps;
std::set<IdxTy> hoods;
// add the hoods containing a piece of each seg
// these need not be oriented 
isin(hoods, p1.x,p1.y,p2.x,p2.y,0);
isin(hoods, p1.x,p1.y,p3.x,p3.y,0);
isin(hoods, p3.x,p3.y,p2.x,p2.y,0);
// get any empty area
fill(hoods);
// not a bit deal but can be cahced... 
GuardType gt( m_left, m_right, m_bottom , m_top , .02);
MM_LOOP(ii,hoods)
{
const IdxTy h=(*ii);
const auto jj=m_map.find(h);
// the point vector 
const auto & hd=(*jj).second.m_p;
MM_SZ_LOOP(k, hd,hdsz)
{
if (gt.guarded(hd[k].x,hd[k].y)) continue; 
if ( gu.contains(hd[k].x,hd[k].y,p1,p2,p3,0)) { dq.add(h,k); } //iuf
} // k  
} // ii 
return 0;
} // DeleteTri

IdxTy Delete( delete_que_type & dq, const IdxTy flags)
{
// del que not bigger than hood...
return dq.remove(m_map); 
} // Delete
template <class Ty> 
IdxTy Add(const Ty & ps, const IdxTy src, const rule_list_type & rl, const IdxTy flags)
{
const bool include_oob=Bit(flags,8); // include those out of bounds  
GuardType gt( m_left, m_right, m_bottom , m_top , .02);
#if 0 
D ex=m_right-m_left;
D ey=m_top-m_bottom;
D guard=.02*((ex<ey)?ex:ey); 
D gl=m_left+guard;;
D gr=m_right-guard;
D gt=m_top-guard;
D gb=m_bottom+guard;
#endif

MM_LOOP(ii,ps)
{
// TODO FIXME this is the etchfront point, just need a dman point class
//const Point & p=(*ii);
const auto  & p=(*ii);
if (!include_oob)  { 
const bool  iguarded=gt.guarded(p.x,p.y); // (p.x>gr)||(p.x<gl)||(p.y>gt)||(p.y<gb);
if (iguarded) continue; 
// actual oob points
if ((p.x<m_left)||(p.x>m_right)||(p.y<m_bottom)||(p.y>m_top))
{
continue; 
} // out 
} // oob

const IdxTy n=hood(p.x,p.y);
bool ok_to_place=true;
bool placed=false;
MM_LOOP(jj,rl)
{
const rule_type & r=(*jj);
ok_to_place=r.place();
IdxSet blocks;
// distance sequred to avoid sqrt later. 
D dlim2=r.a*r.a;
//IdxTy nfu=
// get all blocks that could contain a conflicting point 
// TODO this should distinguish zero from other sources

BlocksInHood(blocks,p,r.a,0);
MM_LOOP(kk,blocks)
{
// this is a const itor ... 
// really just want to modify and entry not the map?
Liinc ll=m_map.find(*kk);
if (ll==m_map.end()) continue;
// this should be ok its ref o vector not element 
// list of all points in the hood block kk
auto&  pv=(*ll).second.m_p;
// these are the existing points in this block 
MM_NCSZ_LOOP(m,pv,szm) // m_p is a vector, as deleting point doh 
{
const Point  px=pv[m]; // non-ref
// rearranged to print status before continue
if (!r.applies(px.src)) continue;
// if this point is too close to edge to remove
const bool  guarded=gt.guarded(px.x,px.y); // (px.x>gr)||(px.x<gl)||(px.y>gt)||(px.y<gb);
const D dpp2=px.dist2(p);
//if (dpp2<1e-18) MM_ERR(MMPR4(r.applies(px.src),dpp2,dlim2,px.src))
// point is close enough to care
if (dpp2<dlim2)
{
//MM_ERR(MMPR4(dpp2,dlim2,px.src,src))
//if(r.place()) { m_map[n].add(p); break;  }
// should be infrequent lol 
if (r.remove_hit()&&!guarded){
// TODO obviously coarse grids will slow down here lol
// better DS 
{ for(IdxTy md=m; md<szm-1; ++md) pv[md]=pv[md+1]; --szm; pv.resize(szm); }  
}
// for this conflict, don't place this candidate point 
else  { ok_to_place=false; } // break here? 
//if (!r.place()) {  ok_to_place=false; if (!r.remove_hit()) break;  }
} // distance 
} // mm lot desc

//(*ll).second.m_p=pv;
} // kk  
if (ok_to_place&&!placed) //  if(r.place()) 
	{ m_map[n].add(p.x,p.y,src); placed=true; }
} // jj , rule list 
//add(x,y,src); 
} // ii 
return 0; 
} // Add

// TODO FIXME this is just wasteful, better constrains possible 
template <class Tp> 
IdxTy BlocksInHood(IdxSet & blocks, const Tp  & p, const D & d, const IdxTy flags)
{
// presume point is in bounds or result messed up anyway.. 
//const IdxTy block=hood(p.x,p.y);
//const D dd=d*d;
const int nxx=(d/m_dx+.5);
const int nyy=(d/m_dy+.5);
const D & x=p.x;
const D & y=p.y;
const int  nx=hoodx(x);
const int ny=hoody(y);
const int nxmin=(nx>nxx)?(nx-nxx-1):0;
const int nxmax=(nx<(int(m_nx)-nxx-1))?(nx+nxx+1):int(m_nx);
const int nymin=(ny>nyy)?(ny-nyy-1):0;
const int nymax=(ny<(int(m_ny)-nyy-1))?(ny+nyy+1):int(m_ny);
if (false) { MM_ERR(MMPR4(nxx,nyy,nx,ny)<<MMPR4(nxmin,nxmax,nymin,nymax)) } 
for(int i=nxmin; i<=nxmax; ++i)
for(int j=nymin; j<=nymax; ++j)
blocks.insert(i+j*m_nx); 

return nx+int(m_nx)*ny;
} // BlocksInHood
template <class Ty> IdxTy  Fill(Ty & set,  const IdxTy flags) const 
{
const IdxTy inisz=set.size();
//mh.fill(hoods);
IdxTy rc=0;
// this gets the segs traversing could leave void
// in middlee, see if anything is missing...
std::map<IdxTy,IdxTy> map; 
MM_LOOP(ii,set) { map[(*ii)]=1; } 
// these are now in order but this assumes how the hood 
// works need to move to hood.. 
if (map.size()==0) return rc;
int last=(*map.begin()).first;
MM_LOOP(ii,map)
{
int h=(*ii).first;
if (h==BAD)
{
MM_ERR(" danger will robinson "<<MMPR4(h,last,set.size(),map.size()))
} // will die 
int d=h-last;
if (d<=1){ last=h;  continue; } 
// only fill horizonally 
int hy1= last/m_nx;
int hy2= h/m_nx;
if (hy2!=hy1) { last=h; continue;}
for(int g=(last+1); g<h; ++g ){ set.insert(g); }
last=h;
//if (d==m_nx ){ last=h; continue;} 
} // ii 
if (false) MM_ERR(MMPR2(inisz,set.size()))
return rc;
} // Fill

#if 0 
// add all block numbers into set that include a point of the
// line segment /A
template <class Ty> IdxTy  IsInOld(Ty & set, const D & x1,const D & y1, const D & x2, const D & y2, const IdxTy flags) const 
{ 
const IdxTy inisz=set.size();
// boundaries, 
//{ int((x-m_left)/m_dx) + int((y-m_bottom)/m_dy)*m_nx; }
//D hoodxmin(int h ) const { return hx*m_dx+m_left;  } // xmin for hood hx
int  h1x=hoodx(x1);
int  h2x=hoodx(x2);
int  h1y=hoody(y1);
int  h2y=hoody(y2);
// avoid problems with zero slopes. 
int  h1=hood(x1,y1);
int  h2=hood(x2,y2);
set.insert(h1);
set.insert(h2);
int hxmin=(h1x<h2x)?h1x:h2x;
int hxmax=(h1x<h2x)?h2x:h1x;
if ((h1x<0)|| (h2x<0)|| (h2y<0)|| (h1y<0)|| (h1<0)|| (h2<0))
{
MM_ERR(MMPR4(h1x,h2x,h1y,h2y)<<MMPR2(h1,h2)<<MMPR4(x1,y1,x2,y2))

} 
//D x=hx*m_dx+m_left; // xmin for hood hx
// equality includes zero diff , so avoid that. 
if (x2!=x1){
for (int  h=hxmin; h<=hxmax; ++h)
{
D x=hoodxmin(h);
D a=(x-x1)/(x2-x1);
D y=y1+a*(y2-y1);
set.insert(hood(x,y)); 
} // h
}
int hymin=(h1y<h2y)?h1y:h2y;
int hymax=(h1y<h2y)?h2y:h1y;
if (y1!=y2){
for (int  h=hymin; h<=hymax; ++h)
{
D y=hoodymin(h);
D a=(y-y1)/(y2-y1);
D x=x1+a*(x2-x1);
set.insert(hood(x,y)); 
} // h
}
if (false) MM_ERR(MMPR4(h1,h2,set.size(),inisz))
return 0; 
} // IsInOld
#endif

template <class Ty> IdxTy  IsIn(Ty & set, const D & x1,const D & y1, const D & x2, const D & y2, const IdxTy flags) const 
{ 
const IdxTy inisz=set.size();
// boundaries, 
//{ int((x-m_left)/m_dx) + int((y-m_bottom)/m_dy)*m_nx; }
//D hoodxmin(int h ) const { return hx*m_dx+m_left;  } // xmin for hood hx
int  h1x=hoodx(x1);
int  h2x=hoodx(x2);
int  h1y=hoody(y1);
int  h2y=hoody(y2);
// avoid problems with zero slopes. 
int  h1=hood(x1,y1);
int  h2=hood(x2,y2);
set.insert(h1);
set.insert(h2);
//int hxmin=(h1x<h2x)?h1x:h2x;
//int hxmax=(h1x<h2x)?h2x:h1x;
//int hymin=(h1y<h2y)?h1y:h2y;
//int hymax=(h1y<h2y)?h2y:h1y;

if ((h1x<0)|| (h2x<0)|| (h2y<0)|| (h1y<0)|| (h1<0)|| (h2<0))
{
MM_ERR(MMPR4(h1x,h2x,h1y,h2y)<<MMPR2(h1,h2)<<MMPR4(x1,y1,x2,y2))

} 
D a=0;
int dirx=1;
int diry=1;
int bx=0;
int by=0;
if (x2<x1){ bx=1;  dirx=-1; } 
if (x2==x1) dirx=0;
if (y2<y1){ by=1;  diry=-1; } 
if (y2==y1) diry=0;

int hx=h1x; // hxmin;
int hy=h1y; // hymin;
while (a<=1.0)
{
// this depends on ordering 
// hoodx : nx=(x-m_left)/m_dx;
// x=m_dx*nx+m_left; del(x)=m_dx
//  next bound is m_dx(nx+1)+m_left
//  next bound is m_dy(ny+1)+m_bottom
//D x=x1+a*(x2-x1);
//D y=y1+a*(y2-y1);
D xnext=m_dx*(hx+dirx+bx)+m_left;
D ynext=m_dy*(hy+diry+by)+m_bottom;
D axn=(x1!=x2)?(xnext-x1)/(x2-x1):(2);
D ayn=(y1!=y2)?(ynext-y1)/(y2-y1):(2);
if (axn<ayn) { hx+=dirx; a=axn;}
else if (axn>ayn) { hy+=diry; a=ayn;}
else 
{
a=axn;
if (a>1) break; 
MM_ERR(" 4 way cross "<<MMPR2(hx,hy))
set.insert(hoodxy(hx+dirx,hy)); 
set.insert(hoodxy(hx,hy+diry)); 
hx+=dirx; hy+=diry; 
} // 
if (a>1) break; 
set.insert(hoodxy(hx,hy)); 
} // true 
if (false) MM_ERR(MMPR4(h1,h2,set.size(),inisz))
return 0; 
} // IsIn 





// infer params from an existing point source 
template <class Ty> IdxTy  Infer(const Ty & ps, const IdxTy flags)
{
IdxTy rc=0;
const IdxTy sz=ps.size();
if (sz<2)
{
MM_ERR(" too few points "<<MMPR(sz))
return BAD;
}
m_left=ps[0].x;
m_right=ps[0].x;
m_bottom=ps[0].y;
m_top=ps[0].y;
MM_LOOP(ii,ps)
{
const auto & p=(*ii);
const D & x=p.x;
const D & y=p.y;
if ( x>m_right) m_right=x;
if ( x<m_left) m_left=x;
if ( y>m_top) m_top=y;
if ( y<m_bottom) m_bottom=y;
} // ii 
m_dx=(m_right-m_left)/D(m_nx);
m_dy=(m_top-m_bottom)/D(m_ny); 
MM_ERR(MMPR4(m_left,m_right,m_bottom,m_top)<<MMPR2(m_dx,m_dy))
return rc;
} // Infer

IdxTy SavePoints(const StrTy & fn, const IdxTy flags)
{
std::ofstream ofs(fn);
return SavePoints(ofs,flags); 

} // SavePoints
IdxTy SavePoints(OsTy & os, const IdxTy flags)
{ 
IdxTy nprec=12;
os<<std::setprecision(nprec);
IdxTy block=0; 
IdxTy serial=0;
MM_LOOP(ii,m_map)
{
MM_LOOP(jj,(*ii).second.m_p)
{
const auto & p=(*jj);
os<<p.x<<" "<<p.y<<" "<<p.src<<" "<<serial<<" "<<block<<CRLF;
++serial;
} // jj 
++block;
} // ii 

return 0; 
}  // SavePoints 

void Init(const StrTy & s, const IdxTy flags) 
{
BaseParams kvp(s);
Init();
kvp.get(m_dx,"dx");
kvp.get(m_dy,"dy");
kvp.get(m_nx,"nx");
kvp.get(m_ny,"ny");
kvp.get(m_left,"left");
m_right=m_left+m_dx*m_nx; // +1 ? 
kvp.get(m_bottom,"bottom");
// FIXME added for boxes may mess up grid
m_top=m_bottom+m_dy*m_ny;
MM_ONCE(dump(),)
} // Init

void Init()
{
// TODO FIXME to avoid div, store also 1/m_dx etc just a good habit.. 
m_dx=1; m_dy=1;
m_nx=100;
m_ny=100;
m_left=0;
m_right=m_left+m_dx*m_nx; // +1 ? 
m_bottom=0;
m_top=0; // not reallyused... 
m_pts=0;
} // Init 
StrTy Dump(const IdxTy flags=0) 
{Ss ss;  
ss<<MMPR4(m_map.size(),m_dx,m_dy,m_nx);
ss<<MMPR4(m_ny,m_left,m_right,m_top);
ss<<MMPR2(m_bottom, m_pts);
return ss.str(); } // Dump

// MEMBERS
LotMap m_map;
D m_dx,m_dy;
IdxTy m_nx,m_ny;
D m_left,m_right,m_top,m_bottom;
IdxTy m_pts;
// added to accomodate hot areas not points,
// TODO should make a new class...
BoxMap m_boxes;
}; // mjm_affine_mesh_hood

//////////////////////////////////////////////

template <class Tr>
class mjm_affine_mesh_hood_map : public std::map<typename Tr::StrTy, mjm_affine_mesh_hood< Tr > >  
{
 typedef mjm_affine_mesh_hood_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_affine_mesh_hood< Tr> >   Super;
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
mjm_affine_mesh_hood_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
static bool Bit(const IdxTy f, const IdxTy b)  { return  ((f>>b)&1)!=0; }
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

}; // mjm_affine_mesh_hood_map




////////////////////////////////////////////
#ifdef  TEST_MJM_AFFINE_MESH_HOOD
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
typedef tester_< mjm_affine_mesh_hood <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_AFFINE_MESH_HOOD "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_affine_mesh_hood<Tr>  Myt;
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

#endif // MJM_AFFINE_MESH_HOOD_H__ 
