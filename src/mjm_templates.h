#ifndef MJM_TEMPLATES_H__
#define MJM_TEMPLATES_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <algorithm>

#include "mjm_globals.h"

namespace mjm
{
// this is a header file, no typedef junk 
//typedef  mjm_generic_traits
typedef double DefaultFloat;
typedef unsigned int IdxTy;
template < class Ty> class sort_second
{ public: int operator()( const Ty & a, const Ty & b) const  
//{ if ( a.second<b.second) return -1; if ( a.second>b.second) return 1; return 0; } 
{ if ( a.second<b.second) return -1;  return 0; } 
};

// sort map m using first or second depending on i and return in d
template <class Ty, class Td> void rank_map(Td & d, const Ty& m, const IdxTy ii)
{
//typedef typename Ty::value_type Tv;
typedef typename std::pair<typename Ty::key_type,typename Ty::mapped_type> Tv;
typedef typename Ty::const_iterator  Ti;
std::vector<Tv> t;
for ( Ti i=m.begin(); i!=m.end(); ++i) { t.push_back(*i); }
std::sort(t.begin(),t.end(),sort_second<Tv>());
for ( IdxTy i=0; i<t.size(); ++i) d[t[i].first]=t.size()-1-i;

}




template <class Tn, class Td, class Ts> Td & cross( Td & d, const Ts & s)
{
const Tn x=s[1]*d[2]-s[2]*d[1];
const Tn y=s[2]*d[0]-s[0]*d[2];
const Tn z=s[0]*d[1]-s[1]*d[0];
d[0]=x;
d[1]=y;
d[2]=z;
return d;
}
template < class Td, class Ts> void cross( Td & d , const Ts & r, const Ts & s)
{
d[0]=s[1]*r[2]-s[2]*r[1];
d[1]=s[2]*r[0]-s[0]*r[2];
d[2]=s[0]*r[1]-s[1]*r[0];
}
template <  class Ts> DefaultFloat dot( const Ts & r, const Ts & s)
{
return s[0]*r[0]+s[1]*r[1]+s[2]*r[2];
}


template <  class Ts> DefaultFloat  distance2(  const Ts & r, const Ts & s)
{
// all this crap should be done in elementwise template functions...
DefaultFloat a= r[0]-s[0] ;
DefaultFloat b= r[1]-s[1] ;
DefaultFloat c= r[2]-s[2] ;

return a*a+b*b+c*c;
}





// actually templates per element would be great here... 
template < class Td, class Ts> void add3( Td & d , const Ts & r, const Ts & s)
{
d[0]=s[0]+r[0];
d[1]=s[1]+r[1];
d[2]=s[2]+r[2];
}
template < class Td, class Ts, class Tf > 
void add3( Td & d , const Ts & r, const Ts & s, const Tf & a, const Tf & b )
{
d[0]=a*s[0]+b*r[0];
d[1]=a*s[1]+b*r[1];
d[2]=a*s[2]+b*r[2];
}

//template <class Td , int n > Td sp(const Td & s ) { return sp<Td,n-1>(s); } 

//template <class Td , int n > Td sp(const Td & s ) { return sp<Td,n-1>(s); } 
template < class Tn, class Ts> void normalize(  Ts & r, const Tn & mag=1)
{
Tn & x=r[0];
Tn & y=r[1];
Tn & z=r[2];
const Tn mold = mag/::sqrt(x*x+y*y+z*z);
x=x*mold;
y=y*mold;
z=z*mold;

}


template <class Ts> void normalize(  Ts & r)
{ normalize<DefaultFloat,Ts>(r); } 

// flush liberal denormals to zed, truncate precision etc
template <class Ts, class Tf> void ddnormals(Ts & s, const Tf & m)
{ if ((s> -m)&&( s< m)) s=0; }
template <class Ts, class Tf> void ddnormal(Ts & s, const Tf & m)
{
ddnormals(s[0],m);
ddnormals(s[1],m);
ddnormals(s[2],m);
}
template <class Ts, class Tf> void ddigits(Ts & s, const Tf & m)
{
std::stringstream ss;
ss.precision(m);
ss<<s[0]<<" "<<s[1]<<" "<<s[2];
ss>>s[0]>>s[1]>>s[2];
}






template <class Tv, class Os, class Tx=const char * > Os & dump(Os & os, Tv & v, const Tx & sep =" ")
{
os<<v[0]<<sep;
os<<v[1]<<sep;
os<<v[2];

}




// advance along directions f1 and f2 ( implied by f1 ) relative to 
// basis vectors d1 and d2. 
template <class Tn, class Td, class Ts> Td & advance(Td & d, const Ts & d1, const Ts & d2, const Tn & f1, const Tn & f2, const Tn & scale)
{
// needs a scale==1 specialization.. 
d[0]+=scale*(d1[0]*f1+d2[0]*f2);
d[1]+=scale*(d1[1]*f1+d2[1]*f2);
d[2]+=scale*(d1[2]*f1+d2[2]*f2);

return d;
}

template <class Tn, class Td, class Ts> Td & advance(Td & d, const Ts & d1, const Ts & d2, const Tn & phi, const Tn & scale)
{
return advance(d,d1,d2,::cos(phi),::sin(phi),scale); 
}

template <class Tn, class Td, class Ts> void theta(Td & d, const Ts & x1, const Ts & x2, const Tn & theta)
{
// note that d and x1 may be the same
Ts y=x1;
normalize(y);
const Tn a=::cos(theta);
const Tn b=::sin(theta);
d[0]=y[0]*a+x2[0]*b;
d[1]=y[1]*a+x2[1]*b;
d[2]=y[2]*a+x2[2]*b;


}


template <class Tc> StrTy toString(const Tc & x)
{
std::stringstream ss;
ss<<x;
return ss.str();

}




}; // mjm 


#endif // guard
