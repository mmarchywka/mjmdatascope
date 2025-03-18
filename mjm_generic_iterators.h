#ifndef MJM_GENERICITERATORS_H__
#define MJM_GENERICITERATORS_H__

#include "mjm_globals.h"
#include "mjm_instruments.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*
Similar to block matrix, a general hierarchial way to iterator reasonably
efficiently over a series of dimenstions  but may not beat hard coded
or other approaches 
 this could be a LOT beter with template nesting level 
Taking functors or template function params would probably make
faster code but thisis easier to write and look at as long
as the ends of the loops are empty which may not be the case
for some summation orders.  
vi mjm_block_matrix.h 
 2330  g++ -DTEST_GENERIC_ITERATORS_MAIN__ -Wall -I.. -x c++ mjm_generic_iterators.h 
 2331  ./a.out

*/
// only for later compatibility
template <int DIM> 
class mjm_generic_itor 
{
typedef mjm_generic_itor Myt;
public:
class Tr {
public:
//typedef Ty D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::Sint Sint;
typedef typename Tr::StrTy StrTy;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::LocTy LocTy;


typedef LocTy location_type;
mjm_generic_itor( const LocTy & dims)
: m_dims(dims),m_sizes(elements(m_dims)),m_top(m_dims.size()-1)//,m_ptr(valloc())
,m_cursor(m_dims.size()),m_done(m_sizes[0]==0)
{}

//mjm_generic_itor( const LocTy & dims)
//: m_dims(dims),m_sizes(elements(m_dims)),m_top(m_dims.size()-1)//,m_ptr(valloc())
//,m_cursor(m_dims.size()),m_done(m_sizes[0]==0)
//{ }
 // could be faster with prealloc etc 
mjm_generic_itor( const IdxTy n)
{ LocTy dims; dims.push_back(n); init(dims); }
// fastest vary on the RIHGT
mjm_generic_itor( const IdxTy n, const IdxTy m)
{ LocTy dims; dims.push_back(n); dims.push_back(m); init(dims); }

mjm_generic_itor( const IdxTy n, const IdxTy m,const IdxTy p)
{ LocTy dims; dims.push_back(n); dims.push_back(m); dims.push_back(p); init(dims); }

mjm_generic_itor( const IdxTy n, const IdxTy m,const IdxTy p, const IdxTy q)
{ LocTy dims;
 dims.push_back(n); dims.push_back(m); dims.push_back(p);dims.push_back(q);
 init(dims); }

IdxTy size() const { return m_top+1; } 

void init( const LocTy & dims)
{
m_dims=dims;
m_sizes=elements(m_dims);
m_top=m_dims.size()-1;
m_cursor=LocTy(m_dims.size());
m_done=(m_sizes[0]==0);

}



//~mjm_block_matrix() { delete [] m_ptr; }
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }

const IdxTy & operator[]( const IdxTy i) const {return m_cursor[i];} 
// aaaarrrgh non const ref but need easy offset 
IdxTy & operator[]( const IdxTy i)  {return m_cursor[i];} 
// this would be better ... 
IdxTy & twiddle( const IdxTy i)  {return m_cursor[i];} 

Myt & operator++( )  {inc(); return *this;} 
const LocTy & cursor() const { return m_cursor; } 


Myt & set_dim(const IdxTy idx, const IdxTy v) { m_dims[idx]=v; return *this; } 

IdxTy inc()
{
IdxTy bumps=0;
IdxTy j=m_top;
++m_cursor[j];
while (m_cursor[j]>=m_dims[j])
{ m_cursor[j]=0; ++bumps;
if (j==0) { m_done=true; break;} 
--j;
++m_cursor[j];
}
return bumps;
}

// requires different init- actually just memset 
// not tested
IdxTy dec()
{
IdxTy bumps=0;
IdxTy j=m_top;
//++m_cursor[j];
//while (m_cursor[j]>=m_dims[j])
while (m_cursor[j]=0)
//{ m_cursor[j]=0; ++bumps;
{
 m_cursor[j]=m_dims[j]-1; ++bumps;
if (j==0) { m_done=true; break;} 
--j;
}
--m_cursor[j];
return bumps;
}

void flush(const IdxTy j)
{
if (m_top<1) { inc(); return; } 
inc(m_top-j);
for(IdxTy i=0; i<j; ++i) m_cursor[m_top-i]=0;

}


// keep separate for later optimization
// havent made non const call by value in a long time... 
// this is dumb since j starts at zero not top 
IdxTy inc( IdxTy j)
{
IdxTy bumps=0;
++m_cursor[j];
while (m_cursor[j]>=m_dims[j])
{ m_cursor[j]=0; ++bumps;
if (j==0) { m_done=true; break;} 
--j;
++m_cursor[j];
}
return bumps;
}


/*

void resize( const LocTy & dims)
{
//	delete [] m_ptr; 
 	m_dims=dims;
	m_sizes=elements(m_dims);
	m_top=(m_dims.size()-1);
//	m_ptr=(valloc());
}
void reform( const LocTy & dims)
{
//	delete [] m_ptr; 
 	m_dims=dims;
	const IdxTy oldsz=m_sizes[0];
	m_sizes=elements(m_dims);
	const IdxTy _top=(m_dims.size()-1);
//	m_ptr=(valloc());
//	if (m_top!=_top) MM_ERR(" reofrming changes size from "<<m_top<<" t " <<_top)
	m_top=_top;
	if (oldsz!=m_sizes[0]) 
		MM_ERR(" reofrming changes size from "<<oldsz<<" t " <<m_sizes[0])
}


*/


const LocTy&  sizes() const { return m_sizes; }

static LocTy loc(const IdxTy * p,const IdxTy sz)
{
LocTy x;
// const IdxTy * pptr=p;
//while (pptr!=0){ x.push_back(*pptr); ++pptr;}
for (IdxTy i=0; i<sz; ++i) x.push_back(p[i]);
return x;
}

//private:
LocTy elements(const LocTy & in) const
{
	const IdxTy sz=in.size();
//	 return 1 lement of value zero so others can get zero 
	if (sz==0) return LocTy(1);
	LocTy szv(sz);
	szv[sz-1]=in[sz-1];
	for (Sint i=sz-2; i>=0; --i) szv[i]=szv[i+1]*in[i];	
return szv;
}
const IdxTy pos(const LocTy & idx, const LocTy & szs) const
{
if (idx.size()<1) return 0;
IdxTy p=idx[m_top];
for (IdxTy i=0; i<m_top; ++i) p+=idx[i]*szs[i+1];
return p;
}
// this need only cary from the fieldincremented 
Sint inc( LocTy & idx, const IdxTy del) const { return inc(idx,del,m_top); } 
Sint inc( LocTy & idx, const IdxTy del, const IdxTy imin) const
{
if (idx.size()<1) return -1;
idx[imin]+=del;
//carry(idx);
carry(idx,imin);
return assert_range(idx,__LINE__);
}

void carry( LocTy & idx) const
{
IdxTy sz=idx.size();
// could fail on entry, not exit 
assert_range(idx,__LINE__);
// not suer what to do with incomplete index array 
if (sz==0 ) return; // fing unsigned fudd 
for (IdxTy i=m_top; i>0; --i)
{
// dangerous with unsigned lol i also zero 
while (idx[i]>=m_dims[i]) { ++idx[i-1]; idx[i]-=m_dims[i]; } 

}
assert_range(idx,__LINE__);

}
void carry( LocTy & idx, const IdxTy imin) const
{
IdxTy sz=idx.size();
// could fail on entry, not exit 
assert_range(idx,__LINE__);
// not suer what to do with incomplete index array 
if (sz==0 ) return; // fing unsigned fudd 
for (IdxTy i=imin; i>0; --i)
{
bool carry=0;
// dangerous with unsigned lol i also zero 
while (idx[i]>=m_dims[i]) { carry=true; ++idx[i-1]; idx[i]-=m_dims[i]; } 
if (!carry ) break; 
}
assert_range(idx,__LINE__);

}




template <class Ta> Sint  assert_range( const LocTy & idx, const Ta & lable="") const
{
Sint rc=0;
IdxTy sz=idx.size();
if (sz!=m_sizes.size()){
 MM_ERR(lable<<" location error "<<m_sizes.size()<<" vs "<<sz)
rc=1; }
for (IdxTy i=0; i<sz; ++i)
{
// dangerous with unsigned lol i also zero 
if (idx[i]>=m_dims[i])
{
MM_ERR(lable<<" index error for "<<i<<" value "<<idx[i]<<" exceeds "<<m_dims[i])
 rc+= (1<<i);
} 

}

return rc;
} //assert_range




const IdxTy postwo(const LocTy & idx, const LocTy & szs) const
{
if (idx.size()<1) return 0;
IdxTy p=idx[0];
for (IdxTy i=1; i<idx.size(); ++i) p+=idx[i]<<szs[i-1];
return p;
}
/*
D * valloc() const 
{
if ( m_dims.size()==0) return 0;
MM_ERR(info_string())
MM_ERR(" size is "<<m_sizes[0])
D * p= new D[m_sizes[0]]; 
::memset(p,0,sizeof(D)*m_sizes[0]);
return p;
}
*/

StrTy info_string() const 
{
SsTy ss;
const IdxTy sz=m_dims.size();
for (IdxTy i=0; i<sz; ++i) ss<<"("<<i<<" "<<m_dims[i]<<" "<<m_sizes[i]<<")";

return ss.str();
}

StrTy string() const 
{
SsTy ss;
const IdxTy sz=m_cursor.size();
for (IdxTy i=0; i<sz; ++i) ss<<" "<<m_cursor[i];
return ss.str();
}

StrTy loc_string(const LocTy & loc, const bool include_paren=true) const 
{
assert_range(loc, "FUDD" );
SsTy ss;
const IdxTy sz=loc.size();
if (include_paren) { ss<<"( "; } 
for (IdxTy i=0; i<sz; ++i) ss<<loc[i]<<" ";
if (include_paren) { ss<<")";}
return ss.str();
}


// each index must be mod this entry
// index zero is mod m_dims[top] 1,m_dims[top-1] etc 
 LocTy m_dims;
 LocTy m_logs;
// these are the total sizes for multiply index -1
// - m_sizes[0] is the toal array size
// m_sizes top multpiples index 1 msizes_top*ind1+ind0 for example 
//  there is no const "1" on top however. 
 LocTy m_sizes;
IdxTy m_top;
//D * m_ptr;
LocTy m_cursor;
bool m_done;

} ; // mjm_block_matrix


class mjm_fixed_sum_itor 
{
typedef mjm_fixed_sum_itor Myt;
public:
class Tr {
public:
//typedef Ty D;
typedef unsigned int IdxTy;
typedef unsigned int MUInt;
typedef double  D;
typedef double  MUD;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::MUInt MUIntTy;
typedef typename Tr::MUD MUDTy;
typedef typename Tr::D D;
typedef typename Tr::Sint Sint;
typedef typename Tr::StrTy StrTy;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::SsTy Ss;
typedef typename Tr::LocTy LocTy;

typedef LocTy location_type;
typedef std::vector<LocTy> PartList;
typedef std::vector<IdxTy> ExcludeList;
//mjm_fixed_sum_itor(const IdxTy n, const MUIntTy & sum )
//: m_sum(sum),m_n(n) {Init(); } 
mjm_fixed_sum_itor(const IdxTy k, const MUIntTy & mu,const MUIntTy & nu )
: m_k(k),m_mu(mu),m_nu(nu),m_excludes(false),m_exclude_level(k+1)  {Init(); } 

mjm_fixed_sum_itor(const IdxTy k, const MUIntTy & mu,const MUIntTy & nu,const IdxTy flags )
: m_k(k),m_mu(mu),m_nu(nu),m_excludes(false),m_exclude_level(k+1)  {Init(flags); } 


mjm_fixed_sum_itor(const IdxTy k, const MUIntTy & mu,const MUIntTy & nu, const ExcludeList & el  )
: m_k(k),m_mu(mu),m_nu(nu),m_excludes(true),m_exclude_level(el)  {Init(); } 


const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
const IdxTy & operator[]( const IdxTy i) const {return m_cursor[i];} 
IdxTy & operator[]( const IdxTy i)  {return m_cursor[i];} 
// this would be better ... 
IdxTy & twiddle( const IdxTy i)  {return m_cursor[i];} 

Myt & operator++( )  {Inc(); return *this;} 
const LocTy & cursor() const { return m_cursor; } 
const IdxTy size() const { return m_cursor.size(); } 
const IdxTy psize() const { return m_parts.size(); } 
const IdxTy ptr() const { return m_ptr; } 
const IdxTy k() const { return m_k; } 

Myt & inc( )  {Inc(); return *this;} 
// the selection of p0 in ref 1 is not compatible with this 
//MUIntTy p0() {if (m_k==0) return 0;  return m_cursor[size()-1]; } 
//MUIntTy p0() {if (m_k==0) return 0;  return m_cursor[0]; } 
//MUIntTy p0() {if (m_k==0) return 0;  return nonzed(); } 
MUIntTy p0() {if (m_k==0) return 0;  return deficit(); } 
IdxTy nonzed()
{
const IdxTy sz=size();

IdxTy p=0;
for(IdxTy i=0; i<sz; ++i) if (m_cursor[i]==0) ++p; 
return p;
}
IdxTy deficit(const bool _p0=true)
{
const IdxTy sz=size();

IdxTy p=0;
// for the p0 calc, this is the eqn in papaer 
if (!_p0) {for(IdxTy i=1; i<sz; ++i) p+=i*m_cursor[i];  }
else { for(IdxTy i=1; i<sz; ++i) p+=m_cursor[i]; }
return m_k-p;
}

//MUIntTy p0() { return m_cursor[0]; } 
// not const may cache 
MUDTy  multiplicity()
{
D mult=1;
const D m=D(m_mu)/D(m_nu);
if (m_k==0) 
{ // MM_MSG(MMPR3(m_done,mult,to_string(m_cursor)))

return mult ; }
// ref [1] needs p0 to go to zero, so this wrong
//const IdxTy max=m_k+1-m_cursor[0];
const IdxTy max=m_k-p0();
for(IdxTy r=0; r<max; ++r) mult*=(m+r);
//const IdxTy sz=size();
for(IdxTy s=1; s<=m_k; ++s) mult/=factorial(m_cursor[s]);
//MM_MSG(MMPR4(m_k,m_done,mult,to_string(m_cursor)))
return mult; 
}
// for tables  and specialization, rihgt now kluges lol 
D factorial( const IdxTy & n) 
{IdxTy j=n;  D v=1; while (j>1) {v=v*j; --j; }return v; }
D power( const D & b,const IdxTy & e) {
D v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; } 

 return v; }
template <class Tv> Tv power_t( const Tv & b,const IdxTy & e) {
Tv v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; } 

 return v; }




private:
static bool Bit(const IdxTy b, const IdxTy f) { return ((1<<b)&f)!=0; } 
void Init(const IdxTy flags=0)
{
const bool use_ml=Bit(0,flags);
SetToFirst();
SetDone();
if (use_ml) MakeList();
else MakeList2();
}
void SetToFirst()
{
SetSize(m_cursor,m_k+1);
/// the whole idea here is that these sum to zero 
//m_cursor[0]=m_k+0;
// instead of diffusing the delta at one, try to diffuse from 
// a single term of m_k 
//if (size()>1) m_cursor[1]=m_k;
//else m_cursor[0]=0;
m_cursor[m_k]=1;
m_ptr=0;
//if(m_k==0) m_cursor[0]=1;

}
void SetDone()
{
// these is always at least one term 
m_done=false;
}


// the p_s is the number of terms used with exponent s
// so this was not quite right
// This means sigma p is k, but we want \sum s*p_s = k
// doh see right after eqn (51)
// and as per the reference [1], \sum p_s= m 
void IncWrongPs()
{
const IdxTy sz=size();
// the impulse at zero spreads out until each bin ==1.
// diffusion occurs at lowest >1,
if (m_k==0) { m_done=true; return; } 
if (m_cursor[sz-1]==1) { m_done=true; }
else
{
IdxTy pos=0;
while (m_cursor[pos]==1) { ++pos; }
// pos points to the SRC to flip into
do {
if (pos>=(sz-1)) {m_done=true; break;  }
m_cursor[pos+1]+=1;
m_cursor[pos]-=1;
++pos;
// now if the src is 1, the dest needs to be bumped
} while  (m_cursor[pos-1]==1) ;

}
//MM_MSG(MMPR2(m_done,to_string(m_cursor)))
} // wrong Inc


void IncNotComplete()
{
const IdxTy sz=size();
// initially there is an impulse at m_k of 1, now diffuse into
// sums of smaller pieces until p[1]==m_k
if (m_k==0) { m_done=true; return; } 
if (m_cursor[1]==m_k) { m_done=true;  return; }
IdxTy i=sz-1;
// these can not all be zero although cursor[0] needs to be updated
while (m_cursor[i]==0) { --i; }
// i now points to a non-zero thing to distribute downward.
++m_cursor[1];
++m_cursor[i-1];
--m_cursor[i];
} // Inc()
void Inc()
{
++m_ptr;
if (m_ptr>=m_parts.size()) { m_done=true; return; }
m_cursor=m_parts[m_ptr];

}


// don't know the supported features of the LocTy yet
void SetSize( LocTy & l, const IdxTy sz) 
{
// not sure the others support this just vector 
l = LocTy(sz);
}
void Integrate(LocTy & d, const LocTy & s, const IdxTy flags=0)
{
const IdxTy sz=s.size();
d=LocTy(s);
if (sz==0) return; 
d[0]=s[0];
for(IdxTy i=1; i<sz; ++i) d[i]=d[i-1]+s[i];
}
void IntegrateFrom(LocTy & d, const LocTy & s, const IdxTy pos, const IdxTy flags=0)
{
if (pos==0) {Integrate(d,s,flags); return; } 
const IdxTy sz=s.size();
if (sz==0) return; 
for(IdxTy i=pos; i<sz; ++i) d[i]=d[i-1]+s[i];
}
StrTy to_string(const LocTy & loc)
{
const IdxTy sz=loc.size();
Ss ss;
if (sz<1) return ss.str();
ss<<",["<<0<<"]="<<loc[0];
for (IdxTy i=1; i<sz; ++i) ss<<",["<<i<<"]="<<loc[i];
return ss.str();
}
template <class Ty > int  Constraint(Ty & loc)
{
IdxTy sum=0;
const IdxTy sz=loc.size();
// gthis really should fail fat but assure good hop number
for(IdxTy i=1; i<sz; ++i)  sum+=i*loc[i];
//if (m_k==sum) return 0;
return m_k-sum; 
}
// in the fp burmann code, this is SLOWER with checks

template <class Td, class Tc> 
void MakePartList(Td & parts, Tc & curse, const IdxTy level, const IdxTy total , const IdxTy target)
{
// level 0 does not contribute to constraint 
if (level==0) return; 
// remainder term 
if (level==1) 
{
curse[1]=target-total;
if (curse[1]!=0) if (m_exclude_level[1]!=0) return; 
parts.push_back(curse);
return;
}
// keep trying until constraint violated 
curse[level]=0;
IdxTy ntotal=total;
if (m_exclude_level[level]!=0)
{
if (ntotal<=target) MakePartList(parts,curse,level-1,ntotal,target);
return;
}
while (ntotal<=target)
{
MakePartList(parts,curse,level-1,ntotal,target);
++curse[level];
ntotal+=(level);
}

}
///////////////////////////////////////////////////////////////////////////////////////
template <class Td, class Tc> 
void MakePartListAll(Td & parts, Tc & curse, const IdxTy level, const IdxTy total , const IdxTy target)
{
// level 0 does not contribute to constraint 
if (level==0) return; 
// remainder term 
if (level==1) 
{
curse[1]=target-total;
//if (curse[1]!=0) if (m_exclude_level[1]!=0) return; 
parts.push_back(curse);
return;
}
// keep trying until constraint violated 
curse[level]=0;
IdxTy ntotal=total;
//if (false) //if (m_exclude_level[level]!=0)
//{
//if (ntotal<=target) MakePartList(parts,curse,level-1,ntotal,target);
//return;
//}
while (ntotal<=target)
{
MakePartList(parts,curse,level-1,ntotal,target);
++curse[level];
ntotal+=(level);
}

}


///////////////////////////////////////////////////////////////////////////////////////
// this worked at least before making excludes conditional 
void MakeList2()
{
m_parts.clear();
m_ptr=0;
const IdxTy szi=m_cursor.size();
if (szi<1) return; 
LocTy curse(szi);
//MM_MSG(MMPR2(m_k,szi))
// the fp version appears SLOWER with exludes but shold verify
// DECISIONS are slower than some calculations Danger Will Robinson 
if (m_excludes)  MakePartList(m_parts, curse, szi-1, 0, m_k);
else  MakePartListAll(m_parts, curse, szi-1, 0, m_k);
if (m_parts.size()!=0) m_cursor=m_parts[0];
//MM_MSG(" partition size "<<MMPR(m_parts.size()))

}
// incrementa l recusion
// start with valid one, p1=m_k.
// set deficit=p1, set p1=0
// increment p2, deficit-=2, if deficit<0, reset p2 go to p3 etc

// check all possible just with local limits
// very slow 
void MakeList()
{
// this is now ok, " do not make cursor[m_k]==1 as this is the init one"
// just change m_cursor lol 
//mjm_generic_itor( const LocTy & dims)
typedef mjm_generic_itor<2>  Itor;
const IdxTy szi=m_cursor.size();
LocTy maxs(szi);
for(IdxTy i=1; i<szi; ++i) { maxs[i]=1+(m_k/i);
//MM_MSG(MMPR2(i,maxs[i])) 
}
maxs[0]=1;
Itor ii(maxs);
MM_ONCE("testing faster partitions",)
while (ii.ok())
{
auto & c=ii.cursor();
const int br=Constraint(c);
if (br==0) {
	 m_parts.push_back(c);
// slower... 
//	ii.flush(1); continue; 
}
if (false) { 	Ss ss;
	if (c.size()!=0) ss<<"["<<0<<"]="<<c[0];
	for (IdxTy i=1; i<c.size(); ++i ) ss<<",["<<i<<"]="<<c[i];
	MM_MSG(MMPR2(br,ss.str()))
	}
// if br<0, could really skip a lot 
++ii;
}
m_ptr=0;
if (m_parts.size()!=0) m_cursor=m_parts[0];
}

// cursor needs to sum to this 
//MUIntTy m_sum; // "k"
MUIntTy m_k,m_mu,m_nu; // "k"
//IdxTy m_n; // "m" = mu/nu number of bins. 
LocTy m_cursor; // location
LocTy m_sum ;//sum of less or eq values  
LocTy m_ulimits,m_llimits; // limits if any 
IdxTy m_ptr;
PartList m_parts;
bool m_excludes; // this needs to be a flag int or something 
ExcludeList m_exclude_level;
bool m_done; 
}; // mjm_fixed_sum_itor


template <int DIM> class mjm_fixed_itor 
{
typedef mjm_fixed_itor Myt;
public:
class Tr {
public:
//typedef Ty D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
//typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::Sint Sint;
typedef typename Tr::StrTy StrTy;
typedef typename Tr::SsTy SsTy;
//typedef typename Tr::LocTy LocTy;


class Location
{
public:
void zero() { ::memset(m_dims,0,sizeof(IdxTy)*DIM);}

IdxTy m_dims[DIM];

}; 

Location m_dims;
Location m_cursor;

bool m_done;

}; //mjm_fixed_itor


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

class mjm_rectangle_itor
{
typedef mjm_rectangle_itor Myt;
public:
class Tr {
public:
//typedef Ty D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
typedef std::vector<Sint> LocTy;
}; // Tr

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::Sint Sint;
typedef typename Tr::StrTy StrTy;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::LocTy LocTy;
mjm_rectangle_itor(const Sint olx, const Sint ilx, const Sint iux, const Sint oux,
const Sint oly, const Sint ily, const Sint iuy, const Sint ouy)
{
Init(olx,ilx,iux,oux,oly,ily,iuy,ouy);
}
// alt, just increemtn or store etc. 
const int geo_idx() const { return  i+nxj; } 
const int x() const { return i; }
const int y() const { return j; }
bool ok() const { return m_ok; } 
Sint inc() { return Inc(); } 
private:
enum { OLX=0, ILX=1, IUX=2, OUX=3, OLY=4, ILY=5, IUY=6, OUY=7};
void Init(const Sint olx, const Sint ilx, const Sint iux, const Sint oux,
const Sint oly, const Sint ily, const Sint iuy, const Sint ouy)
{
IdxTy pos=0;
lims[pos]=olx; ++pos;
lims[pos]=ilx; ++pos;
lims[pos]=iux; ++pos;
lims[pos]=oux; ++pos;
lims[pos]=oly; ++pos;
lims[pos]=ily; ++pos;
lims[pos]=iuy; ++pos;
lims[pos]=ouy; ++pos;
nx=lims[3]-lims[0];
ny=lims[7]-lims[4];
shell_xl=lims[ILX]-lims[OLX];
shell_xr=lims[OUX]-lims[IUX];
shell_yl=lims[ILY]-lims[OLY];
shell_yu=lims[OUY]-lims[IUY];
m_pathological=(shell_xl<=0);
m_pathological|=(shell_xr<=0);
m_pathological|=(shell_yl<=0);
m_pathological|=(shell_yu<=0);
if (m_pathological) 
	{ MM_MSG(" iterator labelled as pathological may not be "<<MMPR4(shell_xl,shell_xr,shell_yl,shell_yu)) 
	 MM_ONCE(" iterator labelled as pathological may not be "<<MMPR4(shell_xl,shell_xr,shell_yl,shell_yu),) }
i=lims[0];
j=lims[4];
nxj=nx*j;
Ok();
}
void Ok()
{
// m_ok=!( ((i>=lims[3]) &&(j>=lims[7])));
// but if the shell width is zero, this is still not ok 
m_ok=!( ((j>=lims[7])));
m_ok&=!m_pathological;

}
Sint Inc()
{
// could be at end of the row or piece 
//enum { OLX=0, ILX=1, IUX=2, OUX=3, OLY=4, ILY=5, IUY=6, OUY=7};
//++i; if ((j>=lims[5])&&(j<lims[6])) if (i==lims[1]) i=lims[2];
++i; if ((j>=lims[ILY])&&(j<lims[IUY])) if (i==lims[ILX]) i=lims[IUX];
//if (i>=lims[3])
if (i>=lims[OUX])
{
//i=lims[0];
i=lims[OLX];
++j;
nxj+=nx;
//if (j>=lims[7]) { Ok(); return 2; } 
if (j>=lims[OUY]) { Ok(); return 2; } 
return 1;
}


return 0;
}
// does not init well 
Sint lims[8];
Sint i,j;
bool m_ok,m_pathological;
Sint nx,ny,nxj,shell_xl,shell_xr,shell_yl,shell_yu;
}; // mjm_rectangle_itor;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////






class shooting_iterator
{

/*
Binary or Liner predictive search over single parameter for shooting boudndary values 
from one side to the other. See mjm_shooting7.h for usage. Can count exceptions
in optional object.  
General binary search over ordered or monotonic functions. 
Known obscure user dpenendies:


*/


//typedef ficklanc Tr;
class Tr
{
public:

typedef unsigned int IdxTy;
typedef double D;
//typedef mjm_block_matrix<D> MyBlock;
//typedef mjm_sparse_matrix<D> MySparse;

typedef std::stringstream Ss;

}; // Tr

// this needs a StrTy def somewhere? Where wtf? 
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
//typedef Tr::MyBlock  MyBlock;
//typedef Tr::MySparse MySparse;
typedef Tr::Ss Ss;


public:
// most of this is not inited.. 
	shooting_iterator():m_exception_flags(0), m_iter(0),m_pass(false),m_strategy(1), mp_cm(0)  {}

void cm(CounterMap * pcm ) { mp_cm=pcm; } 
bool params( const D & tgt,const D & inmin, const D & inmax, const D & tmax, const D & tmin, const IdxTy maxiter)
{
m_min=inmin;
m_max=inmax;
m_tgt=tgt;
m_tol_max=tmax;
m_tol_min=tmin;
m_iter_max=maxiter;
Reset();
return true;
}
void strategy( const IdxTy s) { m_strategy=s; } 
enum {
ABORT_OOR= 1 // if the linear search wants a parameter out of given range, revert to binary 

};

bool result( const D & res)
{
	++m_iter;
	m_res=res;
	const D tol=res-m_tgt;
	m_guesses.push_back(val());
	m_errors.push_back(tol);
	if ((tol>m_tol_min) &&(tol<m_tol_max))
	{
		m_pass=true;
		OnExit();
		return true;
	}
	if (!ok())
	{
		DiagnoseFail(res);
		OnExit();
	 	return true;
	}
//	const D eend=(state(sz,m_u)-state(sz-1,m_u))/(m_fixed(sz,m_x)-m_fixed(sz-1,m_x));
//	const D estart=(state(1,m_u)-state(0,m_u))/(m_fixed(1,m_x)-m_fixed(0,m_x));
		if (tol>0) m_guess_max=m_guess;
		else m_guess_min=m_guess;
return false; 
}
D guess() { Divide(); return val(); } 
D val() const { return m_guess; } 
IdxTy iterations() const { return m_iter; } 
bool ok() { return (m_iter<m_iter_max)&&(m_guess_max>m_guess_min); } 


StrTy to_string()
{
Ss ss;
ss<<" guesses ";
ss<<m_guess<<" ";
ss<<m_guess_min<<" ";
ss<<m_guess_max<<" ";
ss<<" targets ";
ss<<m_tgt_min<<" ";
ss<<m_tgt_max<<" ";
ss<<" minmax ";
ss<<m_min<<" ";
ss<<m_max<<" ";
ss<<m_tol_min<<" ";
ss<<m_tol_max<<" ";
ss<<" iters ";
ss<<m_iter<<" ";
ss<<m_iter_max<<" ";

return ss.str();
}

private:
void DiagnoseFail(const D & res )
{
MM_ERR(name+" shooter failed with "<<res<<" for "<<to_string()) 
MM_ERR(name+" tol value "<<(res-m_tgt)<<" not within "<<m_tol_min<<" to "<<m_tol_max<<" with guess "<<m_guess)
}

// this is supposed to be lintear, it should be able to guess the right
// answer after one iter :) 
// could make this a virtual vft faster than switch lol 
void Divide() {

switch (m_strategy)
{

case 0: {Binary(); break; }
case 1: { Linear(); break; }
case 2: { Log(); break; }
default:  m_guess=.5*(m_guess_min+m_guess_max);


}  // switch 

 } // Divide

void Binary() {  m_guess=.5*(m_guess_min+m_guess_max); }
void Log() {  m_guess=exp(.5*(log(m_guess_min)+log(m_guess_max))); }
// lol this was about 1/2 the total time now insignificant 
// however, it also "fixed"  a spurious charge value towards end of array which
// may reflect undiganozed corruption problems. 
void Linear()
{
const IdxTy sz=m_guesses.size();
if (sz<2) { Binary(); return ; }
const IdxTy p1=sz-2;
const IdxTy p2=sz-1;

const D g1=m_guesses[p1];
const D g2=m_guesses[p2];
const D r1=m_errors[p1];
const D r2=m_errors[p2];
const D d=g2-g1; // this can be zero if range is too small 
if (d==0) {  MM_INCP_MSG(mp_cm,"sifailure" )    Binary(); return; }  // tolerance on zero?
const D a=(r2-r1)/(g2-g1);
if (a==0) { Binary(); return; }  // tolerance on zero?
const D b=r1-g1*a;
const D g=-b/a;
//const bool abort_oor=false;
const bool abort_oor=((m_exception_flags&ABORT_OOR)!=0);
if (abort_oor) if ((g>m_guess_max)||(g<m_guess_min)) { MM_INCP_MSG(mp_cm,"sioor") Binary(); return; }  // 
if ((g>m_guess_max)) {m_guess_max=g; MM_INCP_MSG(mp_cm, "simax" )  MM_ERR(" simax "<<g)  }  // 
if ((g<m_guess_min)) {m_guess_min=g; MM_INCP_MSG(mp_cm, "simin" )  MM_ERR(" simin "<<g)  }  // 

m_guess=g; 

}


void Reset()
{
m_iter=0;
m_guess_min=m_min;
m_guess_max=m_max;
m_tgt_max=m_tol_max-m_tgt;
m_tgt_min=m_tol_min-m_tgt;
m_pass=false;
m_res=0; // maybe bad idea, just to define it 
// ConsistencyCheck();
m_guesses.clear();
m_errors.clear();
}
// all the parameters should make sense 
bool ConsistencyCheck()
{ return true; }
/// update the perf counters etc
void OnExit()
{
if (m_strategy==1) if (m_iter!=3) 
{ 
Ss ss; ss<<"siter"<<m_iter;
MM_INCP_MSG(mp_cm, ss.str()) 
MM_MSG(" linear shooting required "<<m_iter<<" shots,  final tol="<< m_errors[m_errors.size()-1])}
}


StrTy name;
IdxTy  m_exception_flags;
// const, essentialy traits of the search 
D m_guess_min, m_guess_max, m_tgt_max, m_tgt_min;
IdxTy m_iter_max;
// state progress of search 
D m_min, m_max, m_tol_min, m_tol_max,m_guess, m_tgt;
IdxTy m_iter;
D m_res;
bool m_pass;
IdxTy m_strategy;
std::vector<D> m_guesses, m_errors;
CounterMap * mp_cm;
} ; // shooting_iterator



template <class Tv,class Tr > class vector_iterator
{ // make an iterator that makes a vecror look like a map with pair 
typedef vector_iterator<Tv,Tr> Myt;
//typedef  collections_traits::Tr Tr;
typedef typename Tr::IdxTy IdxTy;
typedef std::pair<IdxTy, const Tv &> Tpair;
typedef std::vector<Tv> Tvec;
public:
vector_iterator(   Tvec & v) :m_i(0), m_v(v)// ,m_ptr(0,m_v[0])
,m_end(v.size(),m_v[m_v.size()]) {}
vector_iterator(  Tvec & v, const IdxTy i) :m_i(i), m_v(v)
,m_end(v.size(),m_v[m_v.size()]) {}
Tpair  operator*() { return Tpair(m_i,m_v[m_i]); }
const Tpair  operator*() const { return Tpair(m_i,m_v[m_i]); }
Myt  end() const { return Myt(m_v, m_end.first); }
Myt & operator++() { ++m_i; return *this; }
bool operator==(const Myt & that) { return (&m_v[m_i]==&that.m_v[that.m_i]); }
bool operator!=(const Myt & that) { return !((*this)==that); }
int operator-(const Myt & that) {int a=m_i; int b=that.m_i; return a-b; }
Myt  operator=(const Myt & that ) { Myt x(m_v,m_i); return x; }
private:
IdxTy m_i;
Tvec & m_v;
Tpair m_end;
}; // vector_iterator

template <class Tv,class Tr > class const_vector_iterator
{
typedef const_vector_iterator<Tv,Tr> Myt;
//typedef  collections_traits::Tr Tr;
typedef typename Tr::IdxTy IdxTy;
typedef std::pair<IdxTy, const Tv &> Tpair;
typedef std::vector<Tv> Tvec;
public:
const_vector_iterator( const  Tvec & v) :m_i(0), m_v(v)//,m_ptr(0,m_v[0])
,m_end(v.size(),m_v[m_v.size()]) {}
const_vector_iterator(  const Tvec & v, const IdxTy i) :m_i(i), m_v(v) // ,m_ptr(i,m_v[i])
,m_end(v.size(),m_v[m_v.size()]) {}
Tpair  operator*() { return Tpair(m_i,m_v[m_i]) ;}
operator int()const  { return m_i; } 
operator IdxTy()const  { return m_i; } 
const Tpair  operator*()const  { return Tpair(m_i,m_v[m_i]) ;}
Myt  end() const { return Myt(m_v, m_end.first); }
Myt & operator++() { ++m_i;  return *this; }
//bool operator==(const Myt & that) { return (&m_v[m_i]==&that.m_v[that.m_i]); }
bool operator==(const Myt & that) { return (m_i==that.m_i); }
bool operator!=(const Myt & that) { return !((*this)==that); }
int operator-(const Myt & that) {int a=m_i; int b=that.m_i; return a-b; }
Myt  operator=(const Myt & that ) { Myt x(m_v,m_i); return x; }
private:
IdxTy m_i;
const Tvec & m_v;
Tpair m_end;
}; // vector_iterator

////////////////////////////////////////////////////////////////////////


class acgt_string_iterator
{
typedef acgt_string_iterator Myt;
protected:
typedef  unsigned int IdxTy;
typedef std::string StrTy;
//typedef std::vector<Tv> Tvec;
public:
acgt_string_iterator( const  StrTy & v) :m_i(0),m_end(~0), m_string(v),m_cout(false)//,m_ptr(0,m_v[0])
{}
const StrTy &   operator*() { return m_string ;}
//const Tpair  operator*()const  { return Tpair(m_i,m_v[m_i]) ;}
//Myt  end() const { return Myt(m_v, m_end.first); }
Myt & operator++() { incseq(m_string) ;  return *this; }
bool operator==(const Myt & that) { return (m_string==that.m_string); }
bool operator!=(const Myt & that) { return !((*this)==that); }
Myt  operator=(const Myt & that ) { Myt x(m_string); return x; }
operator bool() { return !m_cout;} 
protected:


void incseq(StrTy & seq)
{
const IdxTy sz=seq.length();
char c[sz+2];
c[sz+1]=0;
bool cin=true;
for (IdxTy i=0; i<seq.length(); ++i)
{
const char ci=seq.c_str()[i];
char cn=ci;
if (cin){
cin=false;
if (ci=='A') cn='C';
else if (ci=='C') cn='G';
else if (ci=='G') cn='T';
else if (ci=='T'){  cn='A'; cin=true; }
} // cin

c[i+1]=cn;
}
m_cout=cin;
if (cin)
{
c[0]='A';
seq= StrTy(c);
return;
}
seq= StrTy(c+1);
}

IdxTy m_i,m_end;
StrTy m_string;
bool m_cout;
}; // acgt_string_iterator


class tree_string_iterator : public acgt_string_iterator
{
typedef acgt_string_iterator Super;
typedef tree_string_iterator Myt;
typedef  Super::IdxTy IdxTy;
typedef Super::StrTy StrTy;
typedef std::vector<StrTy> Stems;
typedef std::vector<IdxTy> C;
typedef std::vector<C> Candidates;
//typedef std::vector<Tv> Tvec;
public:
tree_string_iterator( const  StrTy & v) :Super(v),m_active(bad())//,m_ptr(0,m_v[0])
{Init(); }
const IdxTy size() const { return m_stems.size(); } 
operator bool() { return (m_stems.size()!=0);} 
Myt & inc(const bool success) 
{
// gret so that worked now add to it
// can not eerase yet as we need the old ones lol. 
//m_stems.erase(m_stems.begin()+m_active);
//m_cans.erase(m_cans.begin()+m_active);
if (success) 
{
StrTy next=m_string+StrTy("A");
Add(next); // m_stems.push_back(next);
next=m_string+StrTy("C");
Add(next); // m_stems.push_back(next);
next=m_string+StrTy("G");
Add(next); // m_stems.push_back(next);
next=m_string+StrTy("T");
Add(next); // m_stems.push_back(next);
//incseq(m_stems[m_active]);
//m_string=m_stems[m_active];
//return (*this); 
}
// gret so that worked now add to it
m_stems.erase(m_stems.begin()+m_active);
m_cans.erase(m_cans.begin()+m_active);



// TODO FIXME actually this is just a stack, just operate on teh easiest one to delete doh 
if (m_active>=m_stems.size()) m_active=0; 
m_string=m_stems[m_active];
return *this; 
}


C & cans() //const
{
return m_cans[m_active]; 
}
Myt & inc(const std::vector<IdxTy> & v) 
{
return inc(v.size()!=0); 

}
// gret so that worked now add to it


static IdxTy bad() { return ~0; } 
protected:
void Init()
{
Add(m_string); // m_stems.push_back(m_string);
m_active=0; 
}
void Add(const StrTy & s)
{
m_stems.push_back(s);
if (m_active<m_cans.size()) m_cans.push_back(m_cans[m_active]);
else m_cans.push_back(C());

}

Stems m_stems;
Candidates m_cans;
IdxTy m_active;

};  // tree_string_iterator

////////////////////////////////////////////////////////////////////////////
// 2019-11-28 from mjm_bibtex_entry.h
template <class Key, class Value,int Zero=0>
class mjm_vector_map_iter
{
typedef mjm_vector_map_iter Myt;
typedef unsigned int IdxTy;
typedef std::vector<Value> Vec;
typedef std::map<Key, Vec> TgtM;
public:
mjm_vector_map_iter(TgtM & m ) : m_map(m) {
m_ii=m_map.begin();
m_i=0;
see_ok();
}
const IdxTy & pos() const { return m_i; }
IdxTy vsize() { return (*m_ii).second.size(); } 
void pos(const IdxTy p )  { m_i=p; }
const Key & key() const { return (*m_ii).first; } 
const Vec & vec() const { return (*m_ii).second; }

Value & operator*() { return(*m_ii).second[m_i]; }
Myt & operator++() { ++m_i; 
if (m_i>=(*m_ii).second.size()){  ++m_ii; m_i=0; see_ok(); }
return *this;
}
 operator bool () {  return m_ii!=m_map.end(); }
void see_ok()
{
if (Zero==0) return; 
if (m_ii==m_map.end()) return;
if ((*m_ii).second.size()==0)
{
if ((Zero&1)!=0) { MM_ERR(" zero lengtrh vector "<<(*m_ii).first); } 
if ((Zero&2)!=0) {  
while (vsize()==0) { ++m_ii; if (!(*this)) return; } 
}}

}

TgtM & m_map;
IdxTy m_i;
typename TgtM::iterator m_ii;

};
/*
class iter_ord
{
public:
iter_ord(TgtM & m, const OrdVec & o  ) : m_map(m), m_ord(o) {
m_sz=m_ord.size();
m_o=0;
if ( m_sz>0) m_ii=m_map.find(m_ord[0]);
m_i=0;
}
Tgt & operator*() { return(*m_ii).second[m_i]; }
iter_ord & operator++() { ++m_i;
if (m_i>=(*m_ii).second.size()){ ++m_o ;
if (m_o<m_sz)  m_ii=m_map.find(m_ord[m_o]);
 m_i=0; }
return *this;
}
 operator bool () {  return (m_o<m_sz); }
TgtM & m_map;
const OrdVec  & m_ord;
IdxTy m_i,m_o,m_sz;
typename TgtM::iterator m_ii;

};
*/











/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


#ifdef TEST_GENERIC_ITERATORS_MAIN__
int main(int argc, char ** argv)
{
typedef unsigned int  IdxTy;
typedef mjm_generic_itor<10> T;
typedef T::location_type L;
const T::IdxTy  szi[]={4,5,6,7};
const T::IdxTy  els=sizeof(szi)/sizeof(T::IdxTy);
MM_ERR("location")
L s=T::loc(szi,els);
MM_ERR("ctor")
T::IdxTy i=0;
T x(s);
MM_ERR("dtor")
while (x) { MM_ERR(x.loc_string(x.cursor())) ++x; }
x=T(2,3,4,5);
// usually fall through is the best
//while (x) 
// the switch break willnot exit the infinite loop, need (x) for that 
IdxTy bump=3;
while (true){
// I thought I fixed these examples? 
switch (bump)
{
case 5: MM_ERR(" should not get to 5")
// any partial sums need to be included here,
// this based on empty bottoms of loops 
case 4:{  MM_ERR("  4  co now invalid switch break fails  ")   return 0;  }
case 3: { MM_ERR(" 3 carroty out of 3  ")  } 
case 2: MM_ERR("  2 this is the highest one")
case 1: MM_ERR("  1 this is the highest one")
case 0:{  MM_ERR(" case zero ") break ; } 

} ; //switch 
 MM_ERR(x.loc_string(x.cursor()))
bump=x.inc();
 }
return 0;
}
#endif // main

#endif




