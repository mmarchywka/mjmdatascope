#ifndef MJM_PARTITERATORS_H__
#define MJM_PARTITERATORS_H__

#include "mjm_globals.h"
#include "mjm_instruments.h"
//#include "mjm_generic_iterators.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*
g++ -DTEST_PART_ITERATORS_MAIN__ -O0 -gdwarf-3 -std=c++11 -Wall -I.. -x c++ mjm_part_iterators.h  


 g++ -DTEST_PART_ITERATORS_MAIN__ -std=c++11 -Wall -I.. -x c++ mjm_part_iterators.h


*/


/*
2025-02-04 This appears to preserve the sum \sigma i*n[i]
and \sigma n[i]
This works as expected then for finding the number of ways to
split n things into m bins but you still need to permute
the orders etc. So the first output has ome bin of 5 
while the second has a bin of 4 and 1 of 1 adding to 5.
The exponents all add to 6. 
For example, the usage below \sigma i*n[i] = 5 and \sigma n[i]=6.



mjm_orbital_fft.h653  sin=omega=1;realinit=1 _k=2 _nc=6 _nd=0 _nw=2 flags=0
mjm_orbital_fft.h662  psum=6 bins=5 ii.size()=6 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^5 c[1]^0 c[2]^0 c[3]^0 c[4]^0 c[5]^1
mjm_orbital_fft.h671  ii.ptr()=0 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^4 c[1]^1 c[2]^0 c[3]^0 c[4]^1 c[5]^0
mjm_orbital_fft.h671  ii.ptr()=1 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^3 c[1]^2 c[2]^0 c[3]^1 c[4]^0 c[5]^0
mjm_orbital_fft.h671  ii.ptr()=2 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^4 c[1]^0 c[2]^1 c[3]^1 c[4]^0 c[5]^0
mjm_orbital_fft.h671  ii.ptr()=3 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^3 c[1]^1 c[2]^2 c[3]^0 c[4]^0 c[5]^0
mjm_orbital_fft.h671  ii.ptr()=4 ii.psize()=6
mjm_orbital_fft.h670  DumpGroup(gm)= c[0]^2 c[1]^3 c[2]^1 c[3]^0 c[4]^0 c[5]^0
mjm_orbital_fft.h671  ii.ptr()=5 ii.psize()=6

*/


class mjm_partition_itor 
{
typedef mjm_partition_itor Myt;
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
mjm_partition_itor(const IdxTy sum, const IdxTy bins, const IdxTy & maxmul )
: m_k(sum),m_n(bins),m_maxmul(maxmul),m_excludes(!true),m_exclude_level() {Init(); } 

mjm_partition_itor(const IdxTy sum,  const IdxTy & maxmul )
: m_k(sum),m_n(0),m_maxmul(maxmul),m_excludes(!true),m_exclude_level() {
MM_ERR(" dont use this ctor .... ") 
Init(); } 


//mjm_fixed_sum_itor(const IdxTy k, const MUIntTy & mu,const MUIntTy & nu, const ExcludeList & el  )
//: m_k(k),m_mu(mu),m_nu(nu),m_excludes(true),m_exclude_level(el)  {Init(); } 


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
const IdxTy n() const { return m_parts.size(); } 

const IdxTy k() const { return m_k; } 

Myt & inc( )  {Inc(); return *this;} 
MUIntTy p0() {if (m_k==0) return 0;  return deficit(); } 
IdxTy nonzed()
{
const IdxTy sz=size();

IdxTy p=0;
for(IdxTy i=0; i<sz; ++i) if (m_cursor[i]==0) ++p; 
return p;
}
IdxTy sum() const
{
const IdxTy sz=size();
IdxTy p=0;
for(IdxTy i=0; i<sz; ++i) p+=m_cursor[i];
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
const D m=m_maxmul; // D(m_mu)/D(m_nu);
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

StrTy loc_string(const LocTy & loc, const bool include_paren=true) const
{
//assert_range(loc, "FUDD" );
SsTy ss;
const IdxTy sz=loc.size();
if (include_paren) { ss<<"( "; }
for (IdxTy i=0; i<sz; ++i) ss<<loc[i]<<" ";
if (include_paren) { ss<<")";}
return ss.str();
}



private:

void Init()
{
SetToFirst();
SetDone();
//MakeList();
//MakeList2();
//MakeList3();
MakeList4();
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
if (m_excludes) if (curse[1]!=0) if (m_exclude_level[1]!=0) return; 
parts.push_back(curse);
return;
}
// keep trying until constraint violated 
curse[level]=0;
IdxTy ntotal=total;
if (m_excludes) if (m_exclude_level[level]!=0)
{
if (ntotal<=target) MakePartList(parts,curse,level-1,ntotal,target);
return;
}
while (ntotal<=target)
{
MakePartList(parts,curse,level-1,ntotal,target);
if (m_maxmul!=0) if (level>m_maxmul) return; 
++curse[level];
ntotal+=(level);
//ntotal+=1 ; // (level);
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
if (m_maxmul!=0) if (level>m_maxmul) return; 
++curse[level];
ntotal+=(level);
//ntotal+=1; // (level);
}

}

////////////////////////////////////////////////////
// On a quick test this seems to give plausible results
// but have not verfieid it is picking up all the terms as
// many are probably smalll 
// TODO FIXME need to at least verify it does not duplidate terms or produce invalds
void MakeList4()
{
m_parts.clear();
m_ptr=0;
const IdxTy szi=m_cursor.size();
if (szi<1) return; 
//MM_ERR("start  parting ways "<<MMPR4(szi,m_k,m_n,m_parts.size()))
LocTy curse(szi);
// r is the total amount to distribute, m is the number of bins
// index is amount and value is number of bins with that amount
// put all of them into one bin, zero in the others
const IdxTy r=m_k;
const IdxTy m=m_n;
if (m<1) return; // there needs to be at least one bin 
// for one string, all of the residue is on that one string 

bool done=false;
// changing curse[r]= to += would workd again but still need to push and exit 
// although right now the done loop would function ok 
curse[0]=m-1;
 curse[r]+=1; 
// this is not a special case if the loop checks for lead==0
/*
if (r!=0) { curse[r]=1; } 
else { curse[0]=m; 
m_parts.push_back(curse);
done=true;
}
*/

// make the above simpler and let loop work 
bool move_top=true;
bool move_zed=true;
IdxTy lead=r;
// MM_ERR(MMPR4(r,m,lead,loc_string(curse)))
while (!done)
{

LocTy curseold=curse;
m_parts.push_back(curse);
move_top=true;
move_zed=true;
// if the can not  move, move down the top
for (IdxTy i=1; i<lead; ++i)
{
// TODO FIXME  break missing on first run
// this was supposed to break and it oes change results slightly in large occupancies  
if (curse[i]>1) 
	{curse[i-1]+=1; curse[i+1]+=1; curse[i]-=2; move_zed=false; move_top=false; break; }
}
if (move_zed) { if (curse[0]!=0) {curse[0]-=1; curse[1]+=1; }
else done=true;
 } 
if (move_top&&(lead>0)) { curse[lead]-=1; --lead; curse[lead]+=1; }
// really this need to be if nothing changed, 
if (lead==0) done=true; 
// probably do a faster test here for large r... 
if (curse==curseold)  done=true; 
// MM_ERR(MMPR2(lead,loc_string(curse)))

if (done) break;

}


//MM_ERR(" parting ways "<<MMPR3(m_k,m_n,m_parts.size()))

//else  MakePartListAll(m_parts, curse, szi-1, 0, m_k);
if (m_parts.size()!=0) m_cursor=m_parts[0];
//MM_MSG(" partition size "<<MMPR(m_parts.size()))

} // MakeList4

/////////////////////////////////////////////////////
void MakeList3()
{
m_parts.clear();
m_ptr=0;
const IdxTy szi=m_cursor.size();
if (szi<1) return; 
MM_ERR("start  parting ways "<<MMPR3(m_k,m_n,m_parts.size()))
LocTy curse(szi);
// r is the total amount to distribute, m is the number of bins
// index is amount and value is number of bins with that amount
// put all of them into one bin, zero in the others
const IdxTy r=m_k;
const IdxTy m=m_n;
if (m<1) return; // there needs to be at least one bin 
// for one string, all of the residue is on that one string 
if (m==1)
{
curse[r]=1;
m_parts.push_back(curse);
if (m_parts.size()!=0) m_cursor=m_parts[0];
return; 
}

if (m==2) // again sum(v)=m, r=sum(v*i)
{
const IdxTy rhalf=r>>1;
for (IdxTy j=0; j<=rhalf; ++j)
{
LocTy xcurse(szi);
// the += works for the even case of a middle 
xcurse[r-j]=1; xcurse[j]+=1;
m_parts.push_back(xcurse);
}
if (m_parts.size()!=0) m_cursor=m_parts[0];

return; 
}



LocTy cmax(szi);
MM_SZ_LOOP(i,cmax,csz) { 
if (i==0) {cmax[i]=m; continue; }
const IdxTy rmax=r/i;
if (rmax>m) cmax[i]=m; else cmax[i]=rmax;
}
//curse[r]=1; curse[0]=m-1;
//if (r>m) {}
//else { curse[1]=m; curse[0]=m-r;}
// current values for the curse pointer  s
IdxTy cm=0; // curse[1];
IdxTy cr=0; // curse[1];
//m_parts.push_back(curse);
while (true)
{

if ((cm==m)&&(cr==r)){  m_parts.push_back(curse);
if (curse[r]==1) break;
if (r==0) break;

}
bool cin=true;
for (IdxTy j=0; j<szi; ++j)
{
if ( cin ) { ++curse[j]; ++cm; cr+=j; } 
if (curse[j]>cmax[j])
{
cin=true; cm-=curse[j]; cr-=curse[j]*j; curse[j]=0;

} else {cin=false; } 
if (!cin ) break; 
} // j 
if (cin ) break; 

}

MM_ERR(" parting ways "<<MMPR3(m_k,m_n,m_parts.size()))

//MM_ERR(" ASSFUDD "<<MMPR2(m_excludes,loc_string(curse)))
//MM_MSG(MMPR2(m_k,szi))
// the fp version appears SLOWER with exludes but shold verify
// DECISIONS are slower than some calculations Danger Will Robinson 
//if (m_excludes)  MakePartList(m_parts, curse, szi-1, 0, m_k);
//else  MakePartListAll(m_parts, curse, szi-1, 0, m_k);
if (m_parts.size()!=0) m_cursor=m_parts[0];
//MM_MSG(" partition size "<<MMPR(m_parts.size()))

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
//MM_ERR(" ASSFUDD "<<MMPR2(m_excludes,loc_string(curse)))
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

#if MJM_INCLUDE_SLOW_EXHAUSTIVE_JUNK
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

#endif


// cursor needs to sum to this 
//MUIntTy m_sum; // "k"
MUIntTy m_k,m_n,m_maxmul; //m_mu,m_nu; // "k"
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




/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


#ifdef TEST_PART_ITERATORS_MAIN__

int main(int argc, char ** argv)
{
//typedef mjm_partition_itor<10> T;
typedef mjm_partition_itor T;
typedef unsigned int IdxTy;
//typedef T::location_type L;
//const T::IdxTy  szi[]={4,5,6,7};
//const T::IdxTy  els=sizeof(szi)/sizeof(T::IdxTy);
MM_ERR("location")
//L s=T::loc(szi,els);
//T::IdxTy i=0;
const IdxTy n=atoi(argv[1]);
const IdxTy m=atoi(argv[2]);

MM_ERR("ctor"<<MMPR2(n,m))
T x(n,m);
MM_ERR("dtor")
//while (x) { MM_ERR(x.loc_string(x.cursor())) ++x; }
//x=T(2,3); // 4,5);
// usually fall through is the best
//while (x) 
// the switch break willnot exit the infinite loop, need (x) for that 
IdxTy bump=3;
while (x.ok()){
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
default:{MM_ERR(" bump "<<MMPR(bump)) } 
} ; //switch 
 MM_ERR(x.loc_string(x.cursor()))
bump=x.inc();
 }
return 0;
}
#endif // main

#endif




