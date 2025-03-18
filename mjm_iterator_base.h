#ifndef MJM_ITERATOR_BASE_H__
#define MJM_ITERATOR_BASE_H__

#include "mjm_globals.h"
#include "mjm_instruments.h"
//#include "mjm_generic_iterators.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*

g++ -DTEST_ITERATOR_TRI_MAIN__ -O0 -gdwarf-3 -std=c++11 -Wall -I.. -x c++ mjm_iterator_base.h 


g++ -DTEST_ITERATOR_BASE_MAIN__ -O0 -gdwarf-3 -std=c++11 -Wall -I.. -x c++ mjm_iterator_base.h  


 g++ -DTEST_PART_ITERATORS_MAIN__ -std=c++11 -Wall -I.. -x c++ mjm_part_iterators.h


*/

// for now go with simple and virtuals 
//template <class Timp> class mjm_iterator_base  : public Timp 
 class mjm_iterator_base  // : public Timp 
{
typedef mjm_iterator_base Myt;
//typedef Timp Super ;
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
mjm_iterator_base( )
: m_k(0),m_n(0),m_maxmul(0) {Init(); } 


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
const IdxTy k() const { return m_k; } 

Myt & inc()  {Inc(); return *this;} 
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

StrTy loc_string(const bool include_paren=true, const bool sparse=false) const
{ return loc_string(m_cursor,include_paren,sparse); } 
StrTy loc_string(const LocTy & loc, const bool include_paren=true, const bool sparse=false) const
{
//assert_range(loc, "FUDD" );
SsTy ss;
const IdxTy sz=loc.size();
if (include_paren) { ss<<"( "; }
if (!sparse) { for (IdxTy i=0; i<sz; ++i) ss<<loc[i]<<" "; } 
else { for (IdxTy i=0; i<sz; ++i) if ( loc[i]!=0) {ss<<"("<<i<<","<<loc[i]<<")"; } } 
if (include_paren) { ss<<")";}
return ss.str();
}



private:

protected:

//virtual void Init() //  =0; 
 void Init() //  =0; 
{
//MM_ERR(" in base class init probably from ctor lol")
m_done=true; // unless derivers fix 

}
virtual void  SetToFirst()
{
SetSize(m_cursor,m_k+1);
m_cursor[m_k]=1;
//m_ptr=0;

}

// For most of these they are computationally intensive
// and a virtual ++ is ok for now. 
virtual void Inc() =0;


void SetDone() { m_done=false; }

void SetSize( LocTy & l, const IdxTy sz) { l = LocTy(sz); }

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


// cursor needs to sum to this 
//MUIntTy m_sum; // "k"
MUIntTy m_k,m_n,m_maxmul; //m_mu,m_nu; // "k"
//IdxTy m_n; // "m" = mu/nu number of bins. 
LocTy m_cursor; // location
//LocTy m_sum ;//sum of less or eq values  
//IdxTy m_ptr;
//PartList m_parts;
//bool m_excludes; // this needs to be a flag int or something 
//ExcludeList m_exclude_level;
bool m_done; 
}; // mjm_fixed_sum_itor


/*
Iteratate over sets of lengths and frequenceis with cursor[i]
indificating grequencey of length i such that the lowerst non-zero
is of lenght lmin=m_idx_min and the total length, sum of f*weights,  is less than
 m_max with weights set to length or index. 
*/

 class mjm_arb_limits_itor  : public mjm_iterator_base  // : public Timp 
{
typedef mjm_arb_limits_itor Myt;
typedef  mjm_iterator_base Super;
typedef std::vector<IdxTy> LocTy;
//typedef Super::LocTy LocTy;

public:
mjm_arb_limits_itor() : Super(){}
const IdxTy n() const { return m_sum;}
const IdxTy cnt() const { return m_cnt;}
const IdxTy highest() const { return m_highest;}

// min length filling with cursor size up to L1 and fill total max mfill  
void string_partitions(const IdxTy lmin, const IdxTy L1, const IdxTy mfill)
{
const IdxTy sz=L1+1;
SetSize(m_cursor,sz);
SetSize(m_weights,sz);
m_idx_min=lmin;
m_sz=sz;
for(IdxTy i=0; i<sz; ++i) {
m_weights[i]=i;} 
m_highest=m_idx_min; // this can include some zeros 
if (mfill==0) m_max=L1; else m_max=mfill;
m_sum=0;
m_cnt=0;
m_min_space=0;
m_done=(m_sum>=m_max);
}

private:
void Init()
{
//MM_ERR(" in derived init probably wtf")

}
void Inc(const IdxTy i )
{ ++m_cursor[i]; ++m_cnt; m_sum+=m_weights[i];
//if (m_weights[i]<(m_sz-1)) m_min_space+=1;
 }

void Co(const IdxTy i )
{
m_sum-=m_weights[i]*m_cursor[i];
m_cnt-=m_cursor[i];
//if (m_weights[i]<(m_sz-1)) m_min_space-=m_cursor[i];
m_cursor[i]=0;

}
virtual void Inc()
{
IdxTy i=m_idx_min;
Inc(i);
//++m_cursor[i]; ++m_cnt; m_sum+=m_weights[i];
MM_ONCE(" not sure m_sz is right for spae but littl ematter ",)
m_min_space=m_cnt-1;
//if (m_weights[i]<(m_sz-1)) m_min_space+=1;
while ((m_sum+m_min_space)>m_max)
{
Co(i);
//m_sum-=m_weights[i]*m_cursor[i];
//m_cnt-=m_cursor[i];
//if (m_weights[i]<(m_sz-1)) m_min_space-=m_cursor[i];
//m_cursor[i]=0;

++i;
if (i>= m_sz) { m_done=true; return; } 
if (i>m_highest) m_highest=i;
Inc(i);
m_min_space=m_cnt-1;
//++m_cursor[i]; ++m_cnt; m_sum+=m_weights[i];

}

} // Inc

private:

LocTy m_limits;
LocTy m_weights;
//LocTy m_occupied;
IdxTy m_sum,m_max,m_cnt,m_idx_min;
IdxTy m_sz,m_highest;
IdxTy m_min_space;
} ; // mjm_arb_limits_itor


typedef mjm_arb_limits_itor mjm_arbitrary_itor;

/////////////////////////////////////////////////////////////////////



#if 0 
 class mjm_strict_triangle_itor // : public mjm_iterator_base  // : public Timp 
{
typedef mjm_strict_triangle_itor Myt;
///typedef  mjm_iterator_base Super;
typedef unsigned int IdxTy; 
typedef std::vector<IdxTy>  LocTy;
typedef std::stringstream Ss;
public:
mjm_strict_triangle_itor(const LocTy & lengths ) :  m_lengths(lengths)
// , m_sz(c.size()), m_depth(depth) 
{Init(); }
//mjm_triangle_itor(const IdxTy sz, const IdxTy depth ) :m_counts(0), m_sz(sz), m_depth(depth)  { Init(); }


const LocTy & cursor() const { return m_cursor; } 
const IdxTy& size() const { return m_sz;}
const IdxTy n() const { return m_sum+1;}
const IdxTy& lim() const { return m_ilim;}


//const IdxTy& card() const { return m_card;}
//const IdxTy& depth() const { return m_depth;}
//const IdxTy count() const { return 0;}
void reset() { Reset(); } 
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
void inc() { Inc(); } 
Myt & operator++( )  {Inc(); return *this;} 

StrTy to_string(const LocTy & c, const IdxTy flags=0) const
{
const bool print_all=((flags&1)==0);
const bool print_non_zed=!print_all;
Ss ss;
const IdxTy sz=c.size();
if ( print_all)
{
if ( sz>0) { ss<<"["<<0<<"]="<<c[0]; } 
for(IdxTy i=1;i<sz; ++i) { ss<<", ["<<i<<"]="<<c[i]; } 
}
if (print_non_zed)
{
bool add_comma=false; 
for(IdxTy i=0;i<sz; ++i) { 
	if (c[i]!=0) { if (add_comma) ss<<","; 
	ss<<" ["<<i<<"]="<<c[i]; add_comma=true; }  
} 

}

return ss.str();
}
StrTy to_string(const IdxTy flags=0 ) const
{
Ss ss;
ss<<MMPR(to_string(m_cursor,1));
ss<<MMPR(to_string(m_lengths,1));
//if (m_counts!=0) { ss<<MMPR(to_string(*m_counts)); } 
//ss<<MMPR4( m_lowest,m_highest, m_i,m_j);
//ss<<MMPR4(m_m, m_max, m_biggest,m_done);
ss<<MMPR4(m_sz,m_bumps,m_sum,m_done)<<MMPR2(m_max,m_top);
return ss.str();
}

void make_mult_array( int * p )
{
Reset();
if (m_done) return; 
IdxTy lsum[m_top];
//lsum[0]=m_lengths[0];
//for(int i=1;i<int(m_sz) ;++i) { lsum[i]=m_lengths[i]+lsum[i-1]; } 
m_bumps=m_top;
while (ok()) { 
int j=m_top-1;
int lsumo=(m_bumps==m_top)?0:m_lengths[m_cursor[m_top-m_bumps-1]];
MM_ERR(MMPR4(m_top,m_bumps,lsumo,to_string()))

for(int  i=(m_top-m_bumps); i<int(m_top); ++i)
{
const int lsumi= m_lengths[m_cursor[i]];
lsumo+= lsumi; // m_lengths[m_cursor[i]];
lsum[i]=lsumo;
if (lsumi>m_ilim) { j=i; break; } 
//--p[lsumi];
const IdxTy loc=lsumo+i;
MM_ERR(MMPR4(i,lsum[i],m_cursor[i],m_top)<<MMPR2(loc,m_ilim))
if (loc>m_ilim)
{ // reset the reset of the cursor, 
j=i;

break; 
}
const bool parity=((i&1)==0);
if (parity) { ++p[loc]; } else { --p[loc]; } 
//MM_TRACE

} //i 
m_bumps=0;
while ( j>=0) {
MM_ERR(MMPR4(j,m_cursor[j],m_top,m_bumps))
 ++m_cursor[j]; ++m_bumps; 
if ( m_cursor[j]>m_max)
{
if (j==0) { m_done=true; return; }
//m_cursor[j]=m_cursor[j-1]+1; 
} // max 
--j;
} // j 
for(IdxTy k=j+1; k<m_top; ++k) { 
MM_ERR(MMPR4(k,j,m_top,m_sz))
if (k>0) { m_cursor[k]=m_cursor[k-1]+1;  } 
if ( m_cursor[k]>=m_sz ) {m_top=k ; break; } 
} 
} // tryue
//++(*this);

} // make_mult

// min length filling with cursor size up to L1 and fill total max mfill  

private:
void Init()
{
m_sz=m_lengths.size();
m_max=m_sz-1;
m_cnt=m_sz;
m_top=m_sz-1;
m_sum=0;
MM_LOOP(ii,m_lengths) { m_sum+=(*ii); } 
m_ilim=(m_sum>>1);
Reset();
}
void Reset()
{
m_cursor= LocTy(m_sz);
for(IdxTy i=0; i<m_sz; ++i) m_cursor[i]=i;
m_done=(m_sum==0);

}
private:
void Inc()
{
IdxTy i=m_top;
while (m_cursor[i]==m_max) { if ( i==0) { m_done=true; return; } --i; }
m_bumps=0;
const IdxTy pos=1+m_cursor[i];
IdxTy j=0;
while (i!=m_top) 
{if ( pos+j==m_max) { m_top=i; break; } 
++m_bumps;  m_cursor[i]=pos+j; ++i;  ++j; }
}


//const  LocTy*   m_counts;
const  LocTy&  m_lengths;
LocTy  m_cursor;
 IdxTy m_sz,m_sum,m_ilim,m_cnt;
IdxTy m_max,m_top, m_bumps;
bool m_done;

} ; // mjm_strict_triangle_itor

#endif


//////////////////////////////////////////////////////////////////



 class mjm_triangle_itor // : public mjm_iterator_base  // : public Timp 
{
typedef mjm_triangle_itor Myt;
///typedef  mjm_iterator_base Super;
typedef unsigned int IdxTy; 
typedef std::vector<IdxTy>  LocTy;
typedef std::stringstream Ss;
public:
mjm_triangle_itor(const LocTy & c, const IdxTy depth ) :  m_counts(&c), m_sz(c.size()), m_depth(depth) {Init(); }
mjm_triangle_itor(const IdxTy sz, const IdxTy depth ) :m_counts(0), m_sz(sz),
m_depth(depth)  {

Init(); }
const LocTy & cursor() const { return m_cursor; } 
const IdxTy& size() const { return m_sz;}

/*
The card calculation is wrong, but AFAICT it needs to
use the nested sum of integer thing, the example code used 10x4 
> x=10
> (1/24)*x^(4)+(1/4)*x^(3)+(11/24)*x^(2)+(1/4)*x^(1)
[1] 715
> 
./mjm_clustalw_etc.out -cmd "sum-nest 3 1 2" -quit
y1=+(1/2)*x^(2)+(1/2)*x^(1);

y2=+(1/6)*x^(3)+(1/2)*x^(2)+(1/3)*x^(1);

y3=+(1/24)*x^(4)+(1/4)*x^(3)+(11/24)*x^(2)+(1/4)*x^(1);

./a.out | grep h639 | wc
    715    7579   83799

*/

const IdxTy& card() const { return m_card;}
const IdxTy& depth() const { return m_depth;}
//const IdxTy count() const { return 0;}
void reset() { Reset(); } 
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
void inc() { Inc(); } 
Myt & operator++( )  {Inc(); return *this;} 

StrTy to_string(const LocTy & c, const IdxTy flags=0) const
{
const bool print_all=((flags&1)==0);
const bool print_non_zed=!print_all;
Ss ss;
const IdxTy sz=c.size();
if ( print_all)
{
if ( sz>0) { ss<<"["<<0<<"]="<<c[0]; } 
for(IdxTy i=1;i<sz; ++i) { ss<<", ["<<i<<"]="<<c[i]; } 
}
if (print_non_zed)
{
bool add_comma=false; 
for(IdxTy i=0;i<sz; ++i) { 
	if (c[i]!=0) { if (add_comma) ss<<","; 
	ss<<" ["<<i<<"]="<<c[i]; add_comma=true; }  
} 

}

return ss.str();
}
StrTy to_string(const IdxTy flags=0 ) const
{
Ss ss;
ss<<MMPR(to_string(m_cursor,1));
if (m_counts!=0) { ss<<MMPR(to_string(*m_counts)); } 
//ss<<MMPR4( m_lowest,m_highest, m_i,m_j);
//ss<<MMPR4(m_m, m_max, m_biggest,m_done);
ss<<MMPR4(m_sz,m_card,m_d,m_done)<<MMPR(m_max);
return ss.str();
}

// min length filling with cursor size up to L1 and fill total max mfill  

private:
void Init()
{
//if ( m_sz>0) m_card=1; else m_card=0;
// TODO this is clearly wrong but ok for now 
if (m_depth>0)
{
IdxTy df=1; // m_depth-1;
for (IdxTy i=m_sz-1; i>1; --i)  df=df*i;
m_card=df*m_depth;
} else m_card=0;
//MM_LOOP(ii,m_limits) m_card*=(*ii);
m_d=m_depth-1;
m_max=m_sz-1;
Reset();
}
void Reset()
{
m_cursor= LocTy(m_depth);
m_done=(m_card==0);

}
private:
void Inc()
{
IdxTy i=m_d;
while (m_cursor[i]==m_max) { if ( i==0) { m_done=true; return; } --i; }
const IdxTy pos=1+m_cursor[i];
while (i!=m_depth) { m_cursor[i]=pos; ++i; }
}


const  LocTy*   m_counts;
 LocTy  m_cursor;
 IdxTy m_card,m_sz, m_depth;
IdxTy m_max,m_d;
bool m_done;

} ; // mjm_triangle_itor
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
 class mjm_total_spaces_itor // : public mjm_iterator_base  // : public Timp 
{
typedef mjm_total_spaces_itor Myt;
///typedef  mjm_iterator_base Super;
typedef unsigned int IdxTy; 
typedef std::vector<IdxTy>  LocTy;
typedef std::stringstream Ss;
public:
mjm_total_spaces_itor(const LocTy & l, const IdxTy L ) 
:  m_l(l), m_L(L), m_ls(Sum(m_l)), m_ss(m_L-m_ls)
{Init(); }

const LocTy & spaces() const { return m_spaces; } 
//const LocTy & limits() const { return m_limits; } 
const IdxTy& size() const { return m_size;}
//const IdxTy& card() const { return m_card;}
const IdxTy& i() const { return m_i;}
//const IdxTy count() const { return 0;}
void reset() { Reset(); } 
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
void inc() { Inc(); } 
Myt & operator++( )  {Inc(); return *this;} 

StrTy to_string(const LocTy & c, const IdxTy flags=0) const
{
const bool print_all=((flags&1)==0);
const bool print_non_zed=!print_all;
Ss ss;
const IdxTy sz=c.size();
if ( print_all)
{
if ( sz>0) { ss<<"["<<0<<"]="<<c[0]; } 
for(IdxTy i=1;i<sz; ++i) { ss<<", ["<<i<<"]="<<c[i]; } 
}
if (print_non_zed)
{
bool add_comma=false; 
for(IdxTy i=0;i<sz; ++i) { 
	if (c[i]!=0) { if (add_comma) ss<<","; 
	ss<<" ["<<i<<"]="<<c[i]; add_comma=true; }  
} 

}

return ss.str();
}
StrTy to_string(const IdxTy flags=0 ) const
{
Ss ss;
ss<<MMPR(to_string(m_l,1));
ss<<MMPR(to_string(m_spaces,1));
IdxTy sum=0;
MM_LOOP(ii,m_spaces) sum+=(*ii);
//ss<<MMPR(to_string(m_limits));
//ss<<MMPR4( m_lowest,m_highest, m_i,m_j);
//ss<<MMPR4(m_m, m_max, m_biggest,m_done);
ss<<MMPR4(m_size,m_i,m_done,m_dec)<<MMPR4(m_dec,m_ls, m_ss,sum);
return ss.str();
}

// min length filling with cursor size up to L1 and fill total max mfill  

private:
void Init()
{
m_size=m_l.size();
Reset();
}
void Reset()
{
m_i=0;
m_spaces= LocTy(m_size+1);
//m_line_pos= LocTy(m_size);
m_spaces[0]=m_ss;
//m_line_pos.push_back(m_spaces[0]);
//for(IdxTy i=1; i<m_size; ++i) 
//	m_line_pos.push_back(m_line_pos.back()+m_l[i-1]+m_spaces[i]);
m_dec=0;
m_last=m_ss-1;
m_done=(m_ss==0);

}
private:
void Inc()
{
// MM_ERR(MMPR(to_string()))
if (m_spaces[m_size]==m_ss) { m_done=true; return; } 
if (0!=m_spaces[m_dec]) 
{ --m_spaces[m_dec]; ++m_spaces[m_dec+1];  return; }
++m_spaces[m_dec+2];
--m_spaces[m_dec+1];
IdxTy rem=0;
for(IdxTy i=0; i<=(m_dec+1); ++i) { rem+=m_spaces[i]; m_spaces[i]=0; }
m_spaces[0]=rem;
if (rem>0) m_dec=0; else ++m_dec;
}

IdxTy Sum( const LocTy& s) const
{ IdxTy sum=0; MM_LOOP(ii,s) { sum+=(*ii); } return sum; } 

const  LocTy&  m_l;
const IdxTy m_L;
const IdxTy m_ls, m_ss;
LocTy  m_spaces;
//LocTy  m_linepos;
 IdxTy m_i,m_size,m_dec,m_last;
bool m_done;

} ; // mjm_total_spaces_itor





////////////////////////////////////////////////////////////////
 class mjm_ragged_card_itor // : public mjm_iterator_base  // : public Timp 
{
typedef mjm_ragged_card_itor Myt;
///typedef  mjm_iterator_base Super;
typedef unsigned int IdxTy; 
typedef std::vector<IdxTy>  LocTy;
typedef std::stringstream Ss;
public:
mjm_ragged_card_itor(const LocTy & c) :  m_limits(c) {Init(); }
mjm_ragged_card_itor(const LocTy & c,const LocTy & nz)  {
MM_LOOP(ii,nz) { m_limits.push_back(c[*ii]); } 

Init(); }
const LocTy & cursor() const { return m_cursor; } 
const LocTy & limits() const { return m_limits; } 
const IdxTy& size() const { return m_size;}
const IdxTy& card() const { return m_card;}
const IdxTy& sum() const { return m_sum;}
//const IdxTy count() const { return 0;}
void reset() { Reset(); } 
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
void inc() { Inc(); } 
Myt & operator++( )  {Inc(); return *this;} 

StrTy to_string(const LocTy & c, const IdxTy flags=0) const
{
const bool print_all=((flags&1)==0);
const bool print_non_zed=!print_all;
Ss ss;
const IdxTy sz=c.size();
if ( print_all)
{
if ( sz>0) { ss<<"["<<0<<"]="<<c[0]; } 
for(IdxTy i=1;i<sz; ++i) { ss<<", ["<<i<<"]="<<c[i]; } 
}
if (print_non_zed)
{
bool add_comma=false; 
for(IdxTy i=0;i<sz; ++i) { 
	if (c[i]!=0) { if (add_comma) ss<<","; 
	ss<<" ["<<i<<"]="<<c[i]; add_comma=true; }  
} 

}

return ss.str();
}
StrTy to_string(const IdxTy flags=0 ) const
{
Ss ss;
ss<<MMPR(to_string(m_cursor,1));
ss<<MMPR(to_string(m_limits));
//ss<<MMPR4( m_lowest,m_highest, m_i,m_j);
//ss<<MMPR4(m_m, m_max, m_biggest,m_done);
ss<<MMPR4(m_size,m_card,m_sum,m_done);
return ss.str();
}

// min length filling with cursor size up to L1 and fill total max mfill  

private:
void Init()
{
m_size=m_limits.size();
if ( m_size>0) m_card=1; else m_card=0;
MM_LOOP(ii,m_limits) m_card*=(*ii);
Reset();
}
void Reset()
{
m_sum=0;
m_cursor= LocTy(m_limits.size());
m_done=(m_card==0);

}
private:
void Inc()
{
IdxTy i=0;
while (m_limits[i]==m_cursor[i])
{
// wtf transposed order?? 
m_sum-=m_cursor[i];
m_cursor[i]=0;
++i;
if (i==m_size) { m_done=true; return; } 
}

++m_cursor[i];
++m_sum;
}


 LocTy  m_limits;
 LocTy  m_cursor;
 IdxTy m_card,m_size,m_sum;
bool m_done;

} ; // mjm_ragged_card_itor




/////////////////////////////////////////////////////////////////////////////

 class mjm_cursor_pair_itor // : public mjm_iterator_base  // : public Timp 
{
typedef mjm_cursor_pair_itor Myt;
///typedef  mjm_iterator_base Super;
typedef unsigned int IdxTy; 
typedef std::vector<IdxTy>  LocTy;
typedef std::stringstream Ss;
public:
mjm_cursor_pair_itor(const LocTy & c, const IdxTy lowest, const IdxTy highest) 
:  m_cursor(c), m_lowest(lowest) , m_highest(highest) {Init(); }
//const IdxTy n() const { return m_sum;}
//const IdxTy cnt() const { return m_cnt;}
const IdxTy &lowest() const { return m_lowest;}
const IdxTy &highest() const { return m_highest;}
const IdxTy &biggest() const { return m_biggest;}
const IdxTy &i() const { return m_non_zed[m_i];}
const IdxTy &ni() const { return m_cursor[i()] ;}
const IdxTy &j() const { return m_non_zed[m_j];}
const IdxTy &nj() const { return m_cursor[j()] ;}
const IdxTy &m() const { return m_m ;}
const IdxTy &seq_i() const { return m_i ;}
const IdxTy &seq_j() const { return m_j ;}
const LocTy & non_zeds() const { return m_non_zed; } 
const LocTy & values() const { return m_nz_values; } 
const LocTy & cursor() const { return m_cursor; } 
const IdxTy size() const { return m_max;}
const IdxTy count() const { return (m_max*(m_max-1))>>1;}
const IdxTy multiplicity() const { if (m_i==m_j) { return (ni()*(nj()-1))>>1; }
else { return ni()*nj(); }
}
const IdxTy mm1() const { if (m_i==m_j) { return ((ni()-2)*(nj()-1))>>1; }
else { return (ni()-1)*nj(); }
}
void reset() { Reset(); } 
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
void inc() { Inc(); } 
Myt & operator++( )  {Inc(); return *this;} 

StrTy to_string(const LocTy & c, const IdxTy flags=0) const
{
const bool print_all=((flags&1)==0);
const bool print_non_zed=!print_all;
Ss ss;
const IdxTy sz=c.size();
if ( print_all)
{
if ( sz>0) { ss<<"["<<0<<"]="<<c[0]; } 
for(IdxTy i=1;i<sz; ++i) { ss<<", ["<<i<<"]="<<c[i]; } 
}
if (print_non_zed)
{
bool add_comma=false; 
for(IdxTy i=0;i<sz; ++i) { 
	if (c[i]!=0) { if (add_comma) ss<<","; 
	ss<<" ["<<i<<"]="<<c[i]; add_comma=true; }  
} 

}

return ss.str();
}
StrTy to_string(const IdxTy flags=0 ) const
{
Ss ss;
ss<<MMPR(to_string(m_cursor,1));
ss<<MMPR(to_string(m_non_zed));
ss<<MMPR4( m_lowest,m_highest, m_i,m_j);
ss<<MMPR4(m_m, m_max, m_biggest,m_done);
return ss.str();
}

// min length filling with cursor size up to L1 and fill total max mfill  

private:
void Init()
{
//MM_ERR(" in derived init probably wtf")
IdxTy m=0;
for (IdxTy i=m_lowest; i<=m_highest; ++i)
{
const IdxTy mi=m_cursor[i];
//if (mi!=0) {m+=mi; m_non_zed.push_back(i); }
if (mi!=0) {m+=mi; m_nz_values.push_back(mi); m_non_zed.push_back(i); }
}
m_m=m;
m_max=m_non_zed.size();
if (m_max>0) m_biggest=m_non_zed[m_max-1]; else m_biggest=0;
//m_done=(m_max<2);
//m_done=(m_m<1);
Reset();
}
void Reset()
{
m_i=0;
m_j=0;
m_done=(m_m<1);

}
private:
void Inc()
{

if (m_i==m_j)
{
++m_j;
if (m_j==m_max) { m_done=true; return ; } 
m_i=0;
}  else { ++m_i; } 
if (m_i>m_j) { MM_ERR(" FUDD "<<MMPR4(m_i,m_j,m_max,m_done)) } 
// I thought this was supposed to iterate over i<=j ??? 
if (false)
{++m_i;
if ( m_i==m_max) { ++m_j; 
if ( m_j==m_max) { m_done=true; return; } else m_i=0; }
} // false

}


const LocTy & m_cursor;
LocTy m_non_zed;
LocTy m_nz_values;
//LocTy m_occupied;
IdxTy m_lowest,m_highest;
IdxTy m_i,m_j,m_m, m_max,m_biggest;
bool m_done;

} ; // mjm_cursor_pair_itor







////////////////////////////////////////////////////////////////////////////



#ifdef TEST_COUNT_RAG_MAIN__
int main(int argc, char ** argv)
{
typedef unsigned int IdxTy;
typedef std::vector<IdxTy> LocTy;
LocTy c;
for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
typedef mjm_ragged_card_itor Mrc;
Mrc mrc(c);
//std::vector<IdxTy> c;
//c.push_back(
std::map<IdxTy, IdxTy> cnts;
for(mrc.reset(); mrc.ok(); ++mrc)
{
//MM_MSG(mti.to_string());
++cnts[mrc.sum()];
}
MM_LOOP(ii,cnts) { MM_MSG(MMPR2((*ii).first,(*ii).second)) } 
MM_LOOP(ii,c) { MM_MSG("C "<<MMPR((*ii))) } 
return 0; 
}
#endif


#ifdef TEST_ITERATOR_TRI_MAIN__
int main(int argc, char ** argv)
{
typedef mjm_triangle_itor Mti;
//std::vector<IdxTy> c;
//c.push_back(
Mti mti(10,4);

for(mti.reset(); mti.ok(); ++mti)
{
MM_MSG(mti.to_string());
}
return 0; 
}
#endif


#ifdef TEST_ITERATOR_BASE_MAIN__

int main(int argc, char ** argv)
{
//typedef mjm_partition_itor<10> T;
//typedef mjm_iterator_base T;
typedef mjm_arbitrary_itor T;
typedef unsigned int IdxTy;
//typedef T::location_type L;
//const T::IdxTy  szi[]={4,5,6,7};
//const T::IdxTy  els=sizeof(szi)/sizeof(T::IdxTy);
MM_ERR("location need 2 params ")
//L s=T::loc(szi,els);
//T::IdxTy i=0;
const IdxTy n=atoi(argv[1]);
const IdxTy m=atoi(argv[2]);
const IdxTy flags=(argc>3)?atoi(argv[3]):0;
const bool summary=((flags&1)!=0);
MM_ERR("ctor"<<MMPR3(n,m,flags))
//T x(n,m);
T x;
x.string_partitions(n,m);
MM_ERR("dtor")
IdxTy combos=0;
//while (x) { MM_ERR(x.loc_string(x.cursor())) ++x; }
//x=T(2,3); // 4,5);
// usually fall through is the best
//while (x) 
// the switch break willnot exit the infinite loop, need (x) for that 
//IdxTy bump=3;
//for (; x.ok(); bump=x.inc() ) {
for (; x.ok(); x.inc() ) {
//while (x.ok()){
// I thought I fixed these examples? 
/* switch (bump)
{
case 5: MM_ERR(" should not get to 5")
// any partial sums need to be included here,
// this based on empty bottoms of loops 
case 4:{  MM_ERR("  4  co now invalid switch break fails  ")   return 0;  }
case 3: { MM_ERR(" 3 carroty out of 3  ")  } 
case 2: MM_ERR("  2 this is the highest one")
case 1: {break; } MM_ERR("  1 this is the highest one")
case 0:{  // MM_ERR(" case zero ") 

break ; } 
default:{MM_ERR(" bump "<<MMPR(bump)) } 
} ; //switch 
*/
if (summary)
{
++combos;
if ((combos&((1<<16)-1))==0) 
{  MM_ERR(MMPR3(combos, x.cnt(),x.n())<<x.loc_string(x.cursor(),true,true))}
}
else {  MM_ERR(MMPR2(x.cnt(),x.n())<<x.loc_string(x.cursor(),true,true))}

//bump=x.inc();
 }
return 0;
}
#endif // main

#endif




