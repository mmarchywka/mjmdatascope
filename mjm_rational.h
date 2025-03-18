#ifndef MJM_RATIONAL_H__
#define MJM_RATIONAL_H__
#include "mjm_globals.h"

#include <vector>
#include <sstream>
#include <string>
#include <cmath> 
#include <math.h> 
#include <string.h> 
// 64 bit types 
#include <cstdint>


namespace mjm_int_ops
{

typedef uint64_t IntTy; 
double ratio(const IntTy & n, const IntTy & d)
{
//signed long int en,ed;
//const double mn=mpz_get_d_2exp(&en,n.m_n);
//const double md=mpz_get_d_2exp(&ed,d.m_n);
//signed long int de=en-ed;
//double x= mn/md;
//while (de>0) {x=x*2.0; --de; }
//while (de<0) {x=x/2.0; ++de; }
return double(n)/double(d);
}
}; // mjm_int_ops




// FIXME
// FIXME operators not intuitive on boundary cases like + and minus zero compare lol.
//  reduction or constant denominator options are needed. 
// prime factors and analysis would be nice to have
// dump some error vs denominator things etc. 

// FIXME right now the y dim of he grid layout explodes 
// as the thing is never reudced lol  
// this is a more general thing move somewhere....
// normally these are less than 100 but with maybe fine 
// deltas, to say 1e-6
// may be better map keys that doubles lol 
//template <class Tint=uint64_t, class Tbig= unsigned __int128 >

#ifdef MJM_RATIONAL_NAME
#warning  compiling mjm_rational as MJM_RATIONAL_NAME
#else
#define MJM_RATIONAL_NAME(args...)  mjm_rational(args)
#endif

#if MACRONAME
class MJM_RATIONAL_NAME 
#else
class  mjm_rational
#endif
{
#if MACRONAME
typedef MJM_RATIONAL_NAME Myt;  
#else
typedef mjm_rational Myt;
#endif

typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef unsigned int IdxTy;
public:
// while these may be used to defer division, speed is not a factor. 
// these are POSITIVE although n and d signs could beimportant...  
// 4 quadrant non-reduced by GCF,d=0 etc 
//typedef unsigned int IntTy; 
#ifdef MJM_RATIONAL_BIGINT 
typedef  MJM_RATIONAL_BIGINT  IntTy;
// this makes the integration part different frac of leg time 
enum {NEWDBITS=122};
//enum {NEWDBITS=30};
#else
enum {NEWDBITS=30};
typedef uint64_t IntTy; 
#endif
//typedef Tint IntTy; 

// nothing is supported for thi lol 
//typedef unsigned __int128 IntTy;
//typedef unsigned long long int MultTy;

#ifdef MJM_RATIONAL_BIGINT 
typedef  MJM_RATIONAL_BIGINT  MultTy;
#else
typedef unsigned __int128 MultTy;
//typedef Tbig  MultTy;
#endif

typedef  int SinTy; 
typedef double D;
// sign oflow etc etc 
typedef unsigned int FlagTy; 
// SIGN - set to 1 for negative numbers and not defined when n is zero
// NANBIT - set after div by zero 
// HOPITAL - in some cases limits may be yseful not sure yet
// OFLOW - it should be able to detect if overflow occurs when temps
//         are longer than n and d sizes. 
// REDUCE_FLAG - for indicating reduction should occur as a state not trait
enum FLAGS { SIGN=1, NANBIT=2, HOPITAL=4,OFLOW=8, REDUCE_FLAG=16, OF_FLOAT=(1<<16), OF_LATEX=(1<<17), OF_R=(1<<18) }; 
// with redundant representation, lots of choices but may be 
// better as template parameter... 
// REDUCE- whenver there is a change, remove GCF from n and d.
// DEFER - avoid any operations as long as possible
// CONSTD - try to keep denominators constant and propogate to operation results
// SANE_OPERATORS - equals etc do not distinguuish between plus and mnus zero for example.
enum NUMERICS { REDUCE=1, DEFER=2, CONSTD=4, SANE_OPERATORS=8,EXPBASE=10, BAD=~0}; 
// if m_d is zero, maybe something sane is put into hopital 
#if MACRONAME
 MJM_RATIONAL_NAME() // mjm_rational(): m_n(0), m_d(1),m_flags(0) {}
#else
  mjm_rational() // : m_n(0), m_d(1),m_flags(0) {}
#endif
 : m_n(0), m_d(1),m_flags(0) {}
// sign is normal convention 
#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
// MJM_RATIONAL_NAME // mjm_rational(const IntTy & n, const IntTy & d, const bool sign ): m_n(n), m_d(d) ,m_flags(sign?SIGN:0) {}
 (const IntTy & n, const IntTy & d, const bool sign ): m_n(n), m_d(d) ,m_flags(sign?SIGN:0) {}

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
// MJM_RATIONAL_NAME //mjm_rational(const IntTy & n, const IntTy & d): m_n((n<0)?(-n):n), m_d(d) ,m_flags((n<0)?SIGN:0) {}
 (const IntTy & n, const IntTy & d): m_n((n<0)?(-n):n), m_d(d) ,m_flags((n<0)?SIGN:0) {}

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
// MJM_RATIONAL_NAME //mjm_rational(const double & x): m_n(0), m_d(1),m_flags(0) {set(x); }
(const double & x): m_n(0), m_d(1),m_flags(0) {set(x); }

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
(const IntTy  & x): m_n(0), m_d(1),m_flags(0) 
{if (x<0) {m_flags|=SIGN; m_n=-x;} else m_n=x;  }

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
(const int  & x): m_n(0), m_d(1),m_flags(0) 
{if (x<0) {m_flags|=SIGN; m_n=-x;} else m_n=x;  }

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
(const unsigned int  & x): m_n(x), m_d(1),m_flags(0) 
{ }

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
(const StrTy & s): m_n(0), m_d(1),m_flags(0) 
{ 
//	MM_ERR(" string ctor not work "<<s )
//SsTy ss(s);
//ss>>m_n;
set(s);
// if this is a 
}

#if MACRONAME
 MJM_RATIONAL_NAME 
#else
  mjm_rational
#endif
(const StrTy & s, const IdxTy flags): m_n(0), m_d(1),m_flags(0) 
{ 
//	MM_ERR(" string ctor not work "<<s )
//SsTy ss(s);
//ss>>m_n;
set(s,flags);
// if this is a 
}


// flags added, want to preserver quirky behaviour for zero .. 
void set(const StrTy  & s, const IdxTy flags)
{
typedef unsigned int IdxTy;
//	Myt & r=(*this);
const char * p=s.c_str();
IdxTy pc=0;
IdxTy len=s.length();
IntTy n=0;
IntTy n2=1;
IntTy d=1;
IntTy d2=1;

bool num=true;
bool neg=false;
bool dec=false;
while (pc<len)
{
const char c=p[pc];
	if (c=='-') { neg=!neg; ++pc; continue; } 
	if (c=='.') 
{
if (dec==true) break;
dec=true;
++pc; continue;
} // / 

	if (c=='/') 
{
if (num==false) break;
d=0;
num=false;
//neg=false;
dec=false;
++pc; continue;
} // / 

	if (c>='0') if (c<='9')
{
if ( num ) n=n*10+ int(c)-int('0'); 
if (num&&dec) d2=d2*10;
if ( !num ) d=d*10+ int(c)-int('0'); 
if (!num&&dec) n2=n2*10;
++pc; continue;
}  //digit
MM_ERR(" fall through "<<MMPR4(pc,len,c,num)<<MMPR2(neg,dec))
break;
} // pc

m_n=n*n2;
m_d=d*d2;
if (neg) m_n=-m_n;
MM_ERR(MMPR(__FUNCTION__)<<MMPR4(m_n,m_d,d,d2)<<MMPR2(s,flags))
simplify();
MM_ERR(MMPR(__FUNCTION__)<<MMPR4(m_n,m_d,d,d2)<<MMPR2(s,flags))
} // set 2 param 


void set(const StrTy  & s)
{
	Myt & r=(*this);
const char * p=s.c_str();
typedef unsigned int IdxTy;
IdxTy off=0;
 { MM_ERR(" wtf "<<MMPR4(r.m_n,r.m_d,s,off)) } 
while (p[off]!=0)
{
	if (p[off]=='+') break;
	if (p[off]=='x') break;
	if (p[off]=='-') break;
	if (p[off]>='0') if (p[off]<='9') break;

++off;
}
//MM_MSG(MMPR2(r.m_n,r.m_d))
// 2020-01-28 this then leaves whatever which may be 0/0
// so just set it to zero 
r.m_n=0;
r.m_d=1;
// wtf?? 2022-11-16? 
//if ( p[off]==0) return ;
// this only uses shorter things
r.m_n=int(strtol(p+off,0,0));
r.m_d=1;
while (p[off]!=0)
{
	if (p[off]=='/') { ++off; break; } 
	if (p[off]==',') { ++off; break; } 

++off;
}
MM_MSG(MMPR2(r.m_n,r.m_d))
if (r.m_d==0) { MM_ERR(" wtf "<<MMPR4(r.m_n,r.m_d,s,off)) } 
if ( p[off]==0){ MM_ERR(" no denominat "<<MMPR(r.m_n))   return ; } 
r.m_d=int(strtol(p+off,0,0));
if (r.m_d==0) { MM_ERR(" wtf "<<MMPR4(r.m_n,r.m_d,s,off)) } 
MM_MSG("final value "<< MMPR2(r.m_n,r.m_d))
}


// undefined if n or limit is zero 
const SinTy  sign() const  { return ((m_flags&SIGN)==0)?1:-1;} 
const bool  signbit() const  { return ((m_flags&SIGN)!=0);} 
const bool  reduce() const  { return ((m_flags&REDUCE_FLAG)!=0);} 

template <class Tfudd> const Myt & approx(Tfudd & val) const { val= Tfudd(sign())*Tfudd(m_n)/Tfudd(m_d); return *this; }
template <class Tf> Myt & approx(Tf & val) { val= Tf(sign())*Tf(m_n)/Tf(m_d); return *this; }
// TODO this obviously fails in many cases, but need somethign better
//D approx() const { return D(sign())*D(m_n)/D(m_d); } 
// simplify at least assures equals means equals which is not otherwise
// the case 
// the gmp has divsion for example. 
// this garbage fixed in D ctor, 
//D approx() const {Myt z=(*this); z.simplify();  return D(sign())*D(z.m_n)/D(z.m_d); } 
D approx() const {Myt z=(*this); z.simplify();  return double(z); } 
// return ndigits and leave remainder in the object
StrTy some_remainder(const IdxTy ndigits, const IdxTy flags)
{
std::stringstream  ss;
// exponent, limited to 1<<63 :) 
int  nmult=0;
IdxTy i=0;
bool have_period=false;
if (m_d==0) { ss<<" inf "; return ss.str();  } 
if (m_d<0) { m_d=-m_d;  m_n=-m_n; } 
if (m_n<0) { ss<<"-"; m_n=-m_n; } 
while (i<ndigits)
{
int nx=0; 
if (m_n==0) break; 
while (m_n<m_d) { m_n=m_n*10;  ++i;  ++nx; ++nmult; if (have_period) if (i>=ndigits) break;  }
if (!have_period) if (nx) { ss<<"."; have_period=true; } 
while (nx>1) { ss<<"0"; --nx; } 
IntTy q=m_n/m_d;
ss<<q;
if (!have_period) { ss<<"."; have_period=true; } 
m_n=m_n-q*m_d;
++i;
}
if (ss.str()=="") ss<<"0";
while (nmult>0) { m_d=m_d*10; --nmult; } 
simplify();
return ss.str();
} // some_remainder


#ifndef ALLOW_IMPLICIT_TO_DOUBLE
explicit 
#endif
operator double() const {
Myt x=(*this);
x.simplify();
// TODO fix the stupid int code with a global function 
MM_ONCE(" need to fix the double rat code for ints doh",) 
//if (false ) {return double(x.sign())*(x.m_n.ratio(x.m_d)); } 
if (true) { 
std::stringstream  ss;

 ss<<x.some_remainder(20,0); return atof(ss.str().c_str()); }
if (!true ) {return double(x.sign())*(mjm_int_ops::ratio(x.m_n,x.m_d)); } 

// MM_MSG(MMPR3( double(x.sign()),(x.m_d),double(x.m_d)))  
 return double(x.sign())*double(x.m_n)/double(x.m_d); } 
//explicit operator double() const { return  approx(); } 

explicit operator int() const {
Myt x=(*this);
x.simplify();
// TODO fix the stupid int code with a global function 
MM_ONCE(" need to fix the double rat code for ints doh",) 
//if (false ) {return double(x.sign())*(x.m_n.ratio(x.m_d)); } 
if (true ) {return int(x.sign())*(mjm_int_ops::ratio(x.m_n,x.m_d)); } 
// MM_MSG(MMPR3( double(x.sign()),(x.m_d),double(x.m_d)))  
// return double(x.sign())*double(x.m_n)/double(x.m_d); 

} 
//explicit operator double() const { return  approx(); } 


void addn(const IntTy dn)
{
if (!signbit()){ m_n+=dn; return; }
if (dn>m_n){ negate(); m_n=dn-m_n;}
else m_n-=dn;
}
// these are unsigned, need to check magnitudes. 
Myt & step(const IntTy del) {IntTy dn=del*m_d; addn(dn); return * this; }
 Myt & inc(const IntTy del) {IntTy dn=del;  addn(dn); return * this; }

// format the output string fraction as floats so division works lol 
void of_float() { m_flags|= OF_FLOAT; }
void of_latex() { m_flags|= OF_LATEX; }
void of_R(const bool p=true ) { m_flags=p?(m_flags|OF_R):(m_flags& (~OF_R)); }

//Myt & inc(const IntTy del) {m_n+=del; return *this; }
Myt & incd(const IntTy del) {m_d+=del; return *this; }
Myt & negate() {m_flags^=SIGN; return *this; }
Myt & dec(const IntTy del) 
//{if(m_n<del) {negate(); m_n=del-m_n;} else m_n-=del;  return *this; }
{negate(); addn(del); negate();   return *this; }

Myt & back(const IntTy del) 
{const IntTy p=del*m_d;  return dec(p); }

Myt & setflag(const FlagTy x) { m_flags|=x; return *this; } 
Myt & clearflag(const FlagTy x) { m_flags&=~x; return *this; } 
Myt & clearset(const FlagTy x, const bool & pol) 
{if(pol) setflag(x); else clearflag(x);  return *this; } 


Myt & set_reduce() { m_flags|=REDUCE_FLAG; return *this; } 
Myt & clear_reduce() { m_flags&=~REDUCE_FLAG; return *this; } 

// exp is now really an EXPONENT, not a factor 
//void dnexp(const D &  exp)
void dnexp(const int &  exp)
{
if (exp==0) return; 
D f=1;
int aexp=(exp<0)?(-exp):exp;
for (int x=0; x<aexp; ++x) f=f*D(EXPBASE);
//if(exp>0) m_n=m_n*exp;
if(exp>0) m_n=m_n*f;
//else m_d=m_d*(-exp);
else m_d=m_d*(f);
}
// will use this d is _x is in range. 
Myt & setd(const D & _x, const IntTy d, const bool normalx=true)
{
bool neg=(_x<0) ;
clearset(SIGN,neg);
int exp=0;//  D exp=0;
D x=0;
if ( _x!=0)
{
const D __x=(neg)?(-_x):_x;
// put xinto range for easy denominator scanning 
if ( normalx) normalize(x,exp,__x,0);
	//IntTy n=IntTy(x*d);
	// this allows me to use the operator memeber but need to think threw
	// 2017-06-21 adding MMPR support, could make ifdef 
	IntTy n=IntTy(d*x);
	D e1=x-D(n)/D(d);
	D e2=-x+D(n+1)/D(d);
	m_d=d;
	if (e1<e2) m_n=n; else m_n=n+1;
if ( normalx)	dnexp(exp);
} else { m_n=0; m_d=d;}
return *this;
}

Myt & setint(const int  & _x,const int d=1)
{

if (_x<0) m_n=(-_x)*d;
else  m_n=_x*d;

 m_d=d;
// eh, these are both unsigned wtf 
//if ((m_n*m_d)<0) setflag(SIGN); else clearflag(SIGN);
if ((_x*d)<0) setflag(SIGN); else clearflag(SIGN);
return *this;
}
Myt & setuint(const IdxTy  & _x,const IdxTy d=1)
{ m_n=_x*d; m_d=d; clearflag(SIGN); return *this; }

//recipopcvoal
Myt & setrint(const int  & _x,const int d=1)
{
m_n=1; m_d=_x;
if ((m_n*m_d)<0) setflag(SIGN); else clearflag(SIGN);
return *this;
}
Myt & setruint(const IdxTy  & _x,const IdxTy d=1)
{ m_n=1; m_d=_x; clearflag(SIGN); return *this; }





// these stupid assignments are SLOWWWWWW due to search for d....
Myt & set(const D & _x,const IntTy dmax=1000)
{
static std::map<D,Myt> cached;
MM_ONCE(" caching rats may change some other behaviours ... danger ",)
if ( cached.find(_x)!=cached.end())
{
(*this)=cached[_x];
}
else{
D e=0;
// const dmax=1000; // 2-3 sig figs
bool neg=(_x<0) ;
clearset(SIGN,neg);
// no warning???? wtf 
//const D x=(neg)?(-x):x;
const D x=(neg)?(-_x):_x;
int exp=0;//  D exp=0;
if (x!=0)
{
D xnorm=0;
// put xinto range for easy denominator scanning 
// although div 5 destroys exact crap 
normalize(xnorm,exp,x,0);
// dmax is really stupid ,ned something more intellegint
e=rationalize( m_n, m_d, xnorm, dmax);
dnexp(exp);
reduce_few(m_d,m_n);
}
else { m_n=0; m_d=1; } 
//if (e!=0)if ( e> 1e-19)  MM_ERR(" setting "<<_x<<" to "<<string()<<" e= "<<e<<" exp="<<exp)
//if (e!=0) if ((e>1e-12)||(m_d>1000))  MM_ERR(" setting "<<_x<<" to "<<string()<<" e= "<<e<<" exp="<<exp)
if ((e>1e-8)||(m_d>1000000))
{
if (_x==(2*M_PI))
{ MM_ONCE(" 2pi over threshold "<<MMPR4(_x,string(),e,exp), )} 
else if (_x==(M_PI))
{ MM_ONCE(" 1pi over threshold "<<MMPR4(_x,string(),e,exp), )} 
else
  MM_ERR(" setting "<<_x<<" to "<<string()<<" e= "<<e<<" exp="<<exp)

}
 } // cached
cached[_x]=*this;
return *this; 
}
// put between say 1-10
//void normalize(D& norm, D & exp, const D & raw, const FlagTy & flags)
void normalize(D& norm, int & exp, const D & raw, const FlagTy & flags)
{
if (raw!=0) if ((raw>max_float())||(raw<min_float())) { 
MM_ERR(" out of range "<<raw) 
return; 
}
norm=raw;
exp=0; //exp=1;
if (raw==0) return ; 
// this really does destroy the exactness with the 5
// although chaning impl can maybe fix that. 
//while (norm>1) { norm=norm/D(EXPBASE); exp=exp*EXPBASE; }
while (norm>1) { norm=norm/D(EXPBASE); ++exp; }
//if (exp!=1) return; // dont do both for now 
if (exp!=0) return; // dont do both for now 
//while (norm<.1) { norm=norm*D(EXPBASE); exp=exp*EXPBASE;if (exp>0) exp=-exp; }
while (norm<.1) { norm=norm*D(EXPBASE); --exp;if (exp>0) exp=-exp; }
} 


// conservative, could add an exponent lol 
static IntTy bits() { return sizeof(IntTy)*8-2;}
// maximum number of bits to use in denominator search. 
// avioids everything being zero error but meaningless. 
//static IntTy dbits() { return 32-2;}
static IntTy dbits() { return NEWDBITS;}
static IntTy mulbits() { return sizeof(MultTy)*8-2;}
static D max_float() { return D(IntTy(1)<<bits()); }
static D min_float() { return 1.0/D(IntTy(1)<<bits()); }

// # this only works for a,b greater than zero
//    # this will recurse a long time for a,b much different 
//    # and indeed stack oflow... 
//    def gcd(self,a,b):
static D gcd_similar(const IdxTy a, const IdxTy  b) 
{
        if ((b==1) || ( a==1)) return 1;
        if (b==0) return a;
        if ( a==0) return b;
        if (a>b) return gcd_similar(a-b,b);
        return gcd_similar(a,b-a);
}
static void remove_twos(IntTy & n, IntTy & d) 
{
while ((n!=0)&&(d!=0))
{
if ((n&1)!=0) return; 
if ((d&1)!=0) return; 
n>>=1; d>>=1;
} // while 
}

// apply operations suggested by flags to n and d to conform to criteria optionally
// if these came from operating on operands with denominators d1 and d2. 
// more generally this should look at entire operands and result which may have
// already modified one operand. 
static void conform(IntTy & n, IntTy & d, const FlagTy & flags, const IntTy  d1=0, const IntTy d2=0) 
{



}

static void remove_factor(IntTy & n, IntTy & d, const IntTy & f) 
{
while ((n!=0)&&(d!=0))
{
if ((n%f)!=0) return; 
if ((d%f)!=0) return; 
n=n/f; d=d/f;
} // while 
}

IdxTy is_2n(const IntTy & n) const
{
IntTy x=n&(n-1);

if (x==0) 
{
IntTy np=n;
IdxTy p=0;
// should not sign extend, positive... 
while (np!=0) { ++p; np>>=1; } 
// mjm 2017-11-03 need to check this doh TODO FIXME 
return p-1; 
}
return IdxTy(BAD);
}
void better( D & e, IntTy & dbest, IntTy & nbest, const IntTy & d, const D & x)
{
	//IntTy n=IntTy(x*d);
	// infinte precision crap 
	IntTy n=IntTy(d*x);
	D enew=x-D(n)/D(d);
	if (enew<=e) { dbest=d; nbest=n; e=enew; }
	enew=-x+D(n+1)/D(d);
	if (enew<=e) { dbest=d; nbest=n+1; e=enew; }
}

D  rationalize( IntTy & _n, IntTy &_d, const D & x, const IntTy dmax)
{
D e=x; // error 
IntTy dbest=1;
IntTy nbest=1;
// check a few round ones
// again keep DECREASING magnitured 
//for (IntTy i=8; i<28; ++i)
// with 64 bit, almost EVERYTHING is exact with d circa 1<<54 lol. 
for (IntTy i=(dbits()-2); i>7; --i)
better(e,dbest,nbest,(IntTy(1)<<i),x);

// probably just easier here 
IntTy dguess=1000;
for (IntTy i=0; i<4; ++i)
{
better(e,dbest,nbest,(dguess),x);
dguess=dguess*10; 
}

// not sure how to comb this yet, 
// this is brute force, it needs to check the small ones
// last since it is just a dumb magnitude compare... 
for (IntTy i=1; i<=dmax; ++i)
{
	IntTy d=dmax-i;
better(e,dbest,nbest,d,x);
/*	IntTy n=IntTy(x*d);
	D enew=x-D(n)/D(d);
	if (enew<=e) { dbest=d; nbest=n; e=enew; }
	enew=-x+D(n+1)/D(d);
	if (enew<=e) { dbest=d; nbest=n+1; e=enew; }
*/
} // i
_n=nbest;
_d=dbest;
return e;
}
void rationalizetol( IntTy & _n, IntTy &_d, const D & x
	, const IntTy dmax, const D & tol)
{
D e=x; // error 
IntTy dbest=1;
IntTy nbest=1;
// this is brute force, it needs to check the small ones
// last since it is just a dumb magnitude compare... 
for (IntTy i=1; i<=dmax; ++i)
{
	IntTy d=i;
better(e,dbest,nbest,d,x);
/*	IntTy n=IntTy(x*d);
	D enew=x-D(n)/D(d);
	if (enew<=e) { dbest=d; nbest=n; e=enew; }
	enew=-x+D(n+1)/D(d);
	if (enew<=e) { dbest=d; nbest=n+1; e=enew; }
*/
//	if (enew<tol) break; 
	if (e<tol) break; 

} // i
//        gcdc=self.gcd(self.n,self.d)
_n=nbest;
_d=dbest;

}



// unsafe for non const
const IntTy & d() const { return m_d; }
const IntTy & n() const { return m_n; }
const bool equals(const Myt & that)const {
// zero? lol // FUDD signs lol 
// TODO FIXME 2017-11-18 does sign check apply to zero lol?
if (true) if ( m_n!=0) if (that.m_n!=0) 
       		if (sign()!=that.sign()) return false;  
// TODO FIXME this looks like it works but is very slow... 
if (true) return (that.m_d*m_n)==(that.m_n*m_d);
const double x1=that.approx();
const double x2=approx();
const bool goal=(x1==x2);
if (false) if (sign()!=that.sign())
{
if (goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
 
return false; 
}

if (m_d==that.m_d){ const bool res=  (m_n==that.m_n) ;
if (goal!=res) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())

return res; 
}

MultTy p1=MultTy(m_n)*MultTy(that.m_d);
MultTy p2=MultTy(m_d)*MultTy(that.m_n);
if (p1==p2){

if (!goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
 return true; 
}

if (goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())

return false;
}
const bool equals (const D & thatf) const
{
Myt that;
//that.setd(thatf,m_d,true);
that.set(thatf);
return this->equals(that); 
}

Myt operator-() const
{
Myt x;
x=(*this);
x.negate();
return x; 
}
Myt & operator++() { Myt& x=(*this);  x.m_n=x.m_n+x.m_d; return (*this); }
Myt & operator--() { Myt &x=(*this);  x.m_n=x.m_n-x.m_d; return (*this); }

Myt & operator=(const IdxTy n) { Myt x; x.setuint(n); return set(x); }
Myt & operator=(const StrTy & s) { Myt x; x.set(s); return set(x); }
Myt & operator=(const int  n) { Myt x; x.setint(n); return set(x); }
bool  operator<=(const Myt & that) const  { 
return (*this).equals(that)||lt(that); }
bool  operator==(const Myt & that) const  { return (*this).equals(that); }
bool  operator!=(const Myt & that) const { return!( (*this).equals(that)); }

bool operator<(const Myt & that) const { return lt(that); } 
// eh, this also needs to check the sign bit... 
const bool lt(const Myt & that)const {
#ifndef SKIP_RAT_SANITY
const double x1=that.approx();
const double x2=approx();
const bool goal=(x1>x2);
#endif // SKIP_RAT_SANITY
const bool ts=that.signbit();
const bool s=signbit();
// this may not be right with zero 
if ( m_n!=0) if (that.m_n!=0)
       	if ( s&&!ts){
#ifndef SKIP_RAT_SANITY
if (!goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())

#endif // SKIP_RAT_SANITY
 return true; // this is negative, so less than

}

if ( m_n!=0) if (that.m_n!=0) if ( ts&&!s)
{
#ifndef SKIP_RAT_SANITY

if (goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
#endif // SKIP_RAT_SANITY
 return false; // this is postiive , so not less than
}

// signbits should be the same , 
if (m_d==that.m_d){ return  s?(m_n>that.m_n):(m_n<that.m_n) ;}
MultTy p1=MultTy(m_n)*MultTy(that.m_d);
MultTy p2=MultTy(m_d)*MultTy(that.m_n);
if (p1==p2){
#ifndef SKIP_RAT_SANITY
if (goal){
Myt xthis=*this;
Myt xthat=that;;
xthis.simplify();
xthat.simplify();
 //MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
 MM_ERR(" operator error "<<MMPR4(xthis,xthat,goal,(x1-x2))<<MMPR2(x1,x2))
}
#endif // SKIP_RAT_SANITY

 return false; 
}
if (p1<p2) {
#ifndef SKIP_RAT_SANITY
if ((goal==s)&&(that.m_n!=0)&&(m_n!=00))
{
Myt xthis=*this;
Myt xthat=that;;
xthis.simplify();
xthat.simplify();
 //MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
 MM_ERR(" operator error "<<MMPR4(xthis,xthat,goal,(x1-x2))<<MMPR2(x1,x2))


}
#endif // SKIP_RAT_SANITY
//if (goal==s) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())

return !s; 

}

#ifndef SKIP_RAT_SANITY
if (goal!=s) MM_ONCE(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign(),)

#endif // SKIP_RAT_SANITY
return s;
}

bool operator>(const Myt & that) const { return gt(that); } 
const bool gt(const Myt & that)const {
#ifndef SKIP_RAT_SANITY
const double x1=that.approx();
const double x2=approx();
const bool goal=(x1<x2);
#endif // SKIP_RAT_SANITY
const bool ts=that.signbit();
const bool s=signbit();
// this may not be right with zero 
if ( s&&!ts){ 
#ifndef SKIP_RAT_SANITY
if (goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
#endif // SKIP_RAT_SANITY

return !true;} // this is negative, so less than
if ( ts&&!s){
#ifndef SKIP_RAT_SANITY
if (!goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
#endif // SKIP_RAT_SANITY

 return !false;} // this is postiive , so not less than
// signbits should be the same , 
if (m_d==that.m_d){ 

return  s?(m_n<that.m_n):(m_n>that.m_n) ;}
MultTy p1=MultTy(m_n)*MultTy(that.m_d);
MultTy p2=MultTy(m_d)*MultTy(that.m_n);
if (p1==p2) {
#ifndef SKIP_RAT_SANITY
bin_op_assert(true,goal,that);
#endif // SKIP_RAT_SANITY
//if (goal) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())


return false; }
if (p1<p2) {
#ifndef SKIP_RAT_SANITY
bin_op_assert(s,goal,that);
#endif // SKIP_RAT_SANITY
//if (goal!=s) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())


return s; }
#ifndef SKIP_RAT_SANITY
bin_op_assert(!s,goal,that);
#endif // SKIP_RAT_SANITY
//if (goal==s) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign())
return !s;
}



template <class Ty> void bin_op_assert(const Ty & res, const Ty & goal, const Myt & that) const
{
// for the fp add test, this ALWAYS fails but is ok for binary
if (goal==res)  return;
const double x1=that.approx();
const double x2=approx();
if (false)	MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign()<<" goal "<<goal<<" "<<res)

}
// for indefinite length, these can be better. 
Myt  operator*( const Myt   that ) const 
	{ 
//2025-01-25 mikemath fck 
//MM_DIE(" myt operator ")

Myt x=*this; return x.mult(that); } 
// this could be the same as const int.... 
Myt  operator*( const IntTy &  that ) const 
	{ 
//MM_DIE(" int operator ")
Myt x=*this; return x.mult(Myt(that,IntTy(1))); } 
// 2020-01-26
Myt  operator*=( const IntTy &  that ) //  const 
	{ 

//MM_DIE(" int operator ")

 return mult(Myt(that,IntTy(1))); } 


Myt  operator*( const int n ) const 
	{
//MM_DIE(" int operator ")

 Myt x=*this; Myt that; that.setint(n); return x.mult(that); } 
Myt  operator*( const IdxTy  n ) const 
	{
//MM_DIE(" int operator ")

 Myt x=*this; Myt that; that.setuint(n); return x.mult(that); } 
Myt  operator*( const double & d ) const 
	{
//MM_DIE(" int operator ")

 Myt x=*this; Myt that; that.set(d); return x.mult(that); } 

Myt  operator/( const Myt   that ) const 
	{ Myt x=*this; return x.div(that); } 

Myt  operator/( const int n ) const 
	{ Myt x=*this; Myt that; that.setrint(n); return x.mult(that); } 
Myt  operator/( const IdxTy  n ) const 
	{ Myt x=*this; Myt that; that.setruint(n); return x.mult(that); } 
Myt  operator/( const double & d ) const 
	{ Myt x=*this; Myt that; that.set(1.0/d); return x.mult(that); } 



// made non ref in case this and that are the same, need to think lol 
Myt & mult(const Myt  that)
{
#ifndef SKIP_RAT_SANITY
const double x1=that.approx();
const double x2=approx();
const double goal=(x1*x2);
#endif // SKIP_RAT_SANITY
//const int ss=sign()*that.sign();
// these are unsigned,
m_n=m_n*that.m_n;
m_d=m_d*that.m_d;
if (that.sign()<0) negate();
reduce_few(m_d,m_n);

#ifndef SKIP_RAT_SANITY
const double res=approx();
bin_op_assert(res,goal,that);
#endif // SKIP_RAT_SANITY
// this is already in reduce_few() now doh 
//if (reduce()) simplify();
//else if (that.reduce()) simplify();
return *this;

}
Myt & div(const Myt  that)
{
//const double x1=that.approx();
//const double x2=approx();
//const double goal=(x1*x2);
//const int ss=sign()*that.sign();
// these are unsigned,
m_n=m_n*that.m_d;
m_d=m_d*that.m_n;
if (that.sign()<0) negate();
//reduce_few(m_d,m_n);
//const double res=approx();
//bin_op_assert(res,goal,that);
if (reduce()) simplify();
else if (that.reduce()) simplify();
return *this;
}




// make a copy in case this and that are the same and think. 
// I guess operator need not make extra copy though 
Myt  operator-( const Myt   that ) const
{ Myt x=*this; return x.subtract(that); } 


Myt  operator+( const Myt   that ) const
{ Myt x=*this; return x.add(that); } 

Myt&  operator+=( const Myt   that ) //const
{  return (*this).add(that); } 

Myt&  operator-=( const Myt   that ) //const
{  return (*this).subtract(that); } 


//Myt&  operator*=( const Myt   that ) //const
//{  return (*this).add(that); } 



Myt & add(const Myt that)
{
#ifndef SKIP_RAT_SANITY
const double x1=that.approx();
const double x2=approx();
const double goal=(x1+x2);
#endif // SKIP_RAT_SANITY
const int ss=sign()*that.sign();
// these are unsigned,
if (m_d==that.m_d){ 
if(ss<0) 
{
if (that.m_n>m_n) { negate(); m_n=that.m_n-m_n;
#ifndef SKIP_RAT_SANITY
const double res=approx();
bin_op_assert(res,goal,that);
#endif // SKIP_RAT_SANITY
 return *this; }
else m_n-=that.m_n; 
#ifndef SKIP_RAT_SANITY
const double res=approx();
bin_op_assert(res,goal,that);
#endif // SKIP_RAT_SANITY

return *this;
}
// infinote precisions again crp 
//m_n+=ss*that.m_n;
m_n+=that.m_n*ss;
#ifndef SKIP_RAT_SANITY
const double res=approx();
bin_op_assert(res,goal,that);
#endif // SKIP_RAT_SANITY
 return *this ;}
// arggghhh...
//m_n=that.m_d*m_n+that.m_n*m_d;
IntTy n1=that.m_d*m_n;
MultTy n1m=MultTy(that.m_d)*m_n;
IntTy n2=that.m_n*m_d;
MultTy n2m=MultTy(that.m_n)*m_d;


if (ss<0) 
	{ if ( n2>n1) { negate(); m_n=n2-n1; } else m_n=n1-n2; } 
	else m_n=n1+n2; 
//m_d=m_d*that.m_d;
IntTy m_dnew=m_d*that.m_d;
MultTy dd=MultTy(m_d)*that.m_d;
if ((dd!=m_dnew)||(n1m!=n1)||(n2m!=n2)) {
	MM_ERR(" addition error for rats "<<m_n<<" "<<that.m_n<<" old d"<< that.m_d<<" "<<m_d<<" new d"<<m_dnew) } 
// the grid layout appears to oflow quickly with odd fractions lol
m_d=m_dnew;
// probably a better approach is to find common factors first 
// and avoid potential oflow here. ZZ
reduce_few(m_d,m_n);

#ifndef SKIP_RAT_SANITY
const double res=approx();
bin_op_assert(res,goal,that);
#endif // SKIP_RAT_SANITY
//if (goal!=res) MM_ERR(" operator error "<<x1<<" "<<x2<<" n="<<m_n<<" d="<<m_d<<" "<<that.m_n<<" "<<that.m_d<<" signs "<<sign()<<" "<<that.sign()<<" goal "<<goal<<" "<<res)
return *this;
}

Myt & subtract(const Myt &that)
{
Myt x=that;
x.negate();
return add(x);

}
// obviously needs some thought, but just watn code
// to limp along for now 
template <class Ti,class Tj > IdxTy brute_force_reduce(Ti & d,Ti  & n, const Tj &  f)
{
//if ((d==0)||(n==0))
//	MM_ERR(" reducing with zero fwiw "<<n<<" "<<d)
if (d<Ti(f)) return 0 ;
if (n<Ti(f)) return 0; 
IdxTy cnt=0;
while (((d%f)==0)&&((n%f)==0))
{ d=d/f; n=n/f; ++cnt; }
return cnt;
}
template <class Ti> IdxTy brute_force_reduce2(Ti & d,Ti  & n)
{
//if ((d==0)||(n==0))
//	MM_ERR(" reducing with zero fwiw "<<n<<" "<<d)
if (d<2) return 0 ;
if (n<2) return 0 ; 
IdxTy cnt=0;
while (((d&1)==0)&&((n&1)==0)) { ++cnt; d>>=1; n>>=1; } 
return cnt;
}

// wtf, something like this exists https://oeis.org/A000040
template <class Ti> void reduce_few(Ti & d,Ti  & n)
{
if (true) { simplify(); return ; } 
// factors of 10 and 3 are common 
const IdxTy two=brute_force_reduce2(d,n);
const IdxTy three=brute_force_reduce(d,n,3);
const IdxTy five=brute_force_reduce(d,n,5);
IdxTy mers=0;
for (IdxTy i=3; i<10; ++i) mers+=brute_force_reduce(d,n,(1<<i)-1);
for (IdxTy i=3; i<10; ++i) mers+=brute_force_reduce(d,n,(1<<i)+1);
const IdxTy sum=mers+two+three+five;

if (false) if (sum!=0) 
{ MM_ERR(" reduced "<<two<<" "<<three<<" "<<five<<" "<<mers<<" d and n : "<<d<<" "<<n)
}
}

IdxTy had_10e(const IntTy & n, IntTy & foo ) const
{
if (n<0) return 0;
IntTy x=n;
IdxTy e=0;
while (x>0)
{
IntTy mod=x%(IntTy(10));
if (mod!=0) break;
x=x/IntTy(10); ++e;
}
foo=x;
return e;

}

StrTy string_denom() const
{
SsTy ss;
IdxTy p2=is_2n(m_d);
if ((p2>10) && (p2!=IdxTy(BAD))) ss<<"(1<<"<<p2<<")";
else { IntTy foo; IdxTy  ten=had_10e(m_d,foo);
if (ten>5)
{
IntTy mx=(foo);
if ((m_flags&OF_FLOAT)==0)  ss<<mx<<"e"<<ten;
else  ss<<(m_d)<<".";

}
else{
if ((m_flags&OF_FLOAT)==0)  ss<<m_d;
//else  ss<<double(m_d);
else  ss<<(m_d)<<".";
}
}
return ss.str();
}
StrTy string(const IdxTy flags=0) const
{
SsTy ss;
Myt x=*this;
if ((flags&1)!=0) x.simplify();
if ((x.m_flags&OF_R)!=0)  
{
ss<<"(";
if (x.signbit()) ss<<"-";
if ((x.m_flags&OF_FLOAT)==0)  ss<<x.m_n;
else ss<<(x.m_n)<<".";
if (x.m_d!=1) ss<<"/"<<x.string_denom();  
ss<<")";
return ss.str();
}	
if ((x.m_flags&OF_LATEX)!=0)  
{
ss<<"{";
if (x.m_d!=1) ss<<"\\frac{";
if (x.signbit()) ss<<"-";
if ((x.m_flags&OF_FLOAT)==0)  ss<<x.m_n;
else ss<<(x.m_n)<<".";
if (x.m_d!=1) 
{ 
ss<<"}{"<<x.string_denom()<<"}";  
} 
ss<<"}";
return ss.str();

}
// 2022-12-25
//if (x.signbit()) ss<<"-";
//ss<<m_n<<"/"<<m_d; 
if ((x.m_flags&OF_FLOAT)==0)  ss<<x.m_n;
//else ss<<double(m_n)<<".";
else ss<<(x.m_n)<<".";
//ss<<m_n<<"/"<<string_denom(); 
//if (m_d!=1) { 
// FUDD DOH 
if (x.m_d!=1) { 
//MM_MSG(" test md code "<<MMPR2(m_d,(m_d!=1)))
ss<<"/"<<x.string_denom();  } 
return ss.str();

}
// this needs to take either decmials floats or factions.
bool read(const StrTy & input)
{
//SsTy ss(input);
//const char * c=input.c_str();
//const IdxTy sz=strlen(c);


return false;
}


Myt & set(const Myt & that)
{
m_n=that.m_n;
m_d=that.m_d;
m_flags=that.m_flags;
return *this;
}
// need a way to remove gcf but not here... 


// https://stackoverflow.com/questions/28267611/using-euclid-algorithm-to-find-gcfgcd
IntTy gcdslo(IntTy  a, IntTy b) {
    while (a != b) {
        if (a > b) {
            a -= b;
        }
        else {
            b -= a;
        }
    }
    return a;
}

IntTy gcd(IntTy a, IntTy b) {
    return (b == 0) ? a : gcd(b, a % b);
}
void simplify()
{
const IntTy gcx=gcd(m_n,m_d);
if(gcx>1) 
{

//MM_MSG(" reducing "<<MMPR3(m_n,m_d,gcx))
m_n/=gcx;
m_d/=gcx;
}

}

const IntTy & get_n() const { return m_n; } 
const IntTy & get_d() const { return m_d; } 

// this will have lots of special forms- power of 2, 10 etc
// remove common factors, retain an exponent for dynamic range etc.
// see state flags and other things. 
//MEMBERS
IntTy m_n, m_d;
FlagTy m_flags; 

#if MACRONAME
// MJM_RATIONAL_NAME 
friend std::ostream&  operator<<(std::ostream& os, const MJM_RATIONAL_NAME & r) ; // mjm_rational & r);
#else
//  mjm_rational
friend std::ostream&  operator<<(std::ostream& os, const mjm_rational & r) ; // mjm_rational & r);
#endif
//friend std::ostream&  operator<<(std::ostream& os, const MJM_RATIONAL_NAME & r) ; // mjm_rational & r);

}; // mjm_rational
////////////////////////////////////////////////////////////////////////////

#if MACRONAME
std::ostream&  operator<<(std::ostream& os, const MJM_RATIONAL_NAME & r ) // mjm_rational & r)
#else
std::istream&  operator>>(std::istream& is,  mjm_rational & r ) // mjm_rational & r)
{
	std::string s;
	is>>s;
	r=s;

return is;
}

std::ostream&  operator<<(std::ostream& os, const mjm_rational & r ) // mjm_rational & r)
#endif

{
os<<r.string(1);
return os;
}

mjm_rational pow(const mjm_rational & a, const int  p)
{
mjm_rational x=int(1);
int n=(p<0)?(-p):p;
for(int  i=0; i<n; ++i) x=x*a;
if (p<0) return mjm_rational(1)/x;
return x;
}

// rather than casting and calling floats,
// want to try to implement but should not be
//used if code is correct 
mjm_rational sqrt(const mjm_rational & a)
{
MM_ERR(" sqrt not valid for mjm_rational "<<MMPR(a))
return a;
}


mjm_rational sin(const mjm_rational & a)
{
MM_ERR(" sin not valid for mjm_rational "<<MMPR(a))
return a;
}


mjm_rational cos(const mjm_rational & a)
{
MM_ERR(" cos not valid for mjm_rational "<<MMPR(a))
return a;
}





#endif

