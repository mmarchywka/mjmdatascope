
#ifndef MJM_STRINGS_H__
#define MJM_STRINGS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <time.h>

// for testing, use string stream for radix conversion 
// until we get time to go back and write something good
#include <sstream>
#include <iomanip>
#include "mjm_globals.h"

// #include <sys/time.h>
// #include <sys/resource.h>

#include <time.h>
#include <memory.h>

class mjm_strings 
{
private:
typedef std::runtime_error ErTy;
typedef std::runtime_error CmdEr;

typedef unsigned int FlagTy;
typedef unsigned int IdxTy; 
typedef char ChTy;
typedef std::string StrTy;
//typedef unsigned long int UlongTy;
typedef time_t UlongTy;
typedef long long my64int;
public:

// https://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
static StrTy fancy_to_lower(const StrTy & y) 
{
StrTy x=y;
//std::string data = "Abc"; 
std::transform(x.begin(), x.end(), x.begin(), ::tolower);
return x; 
}
// find pattern p in s 
static bool fast_grep(const StrTy &s, const StrTy & p)
{
const char * pp=p.c_str();
const char * sp=s.c_str();
const IdxTy psz=p.length();
const IdxTy ssz=s.length();
if (ssz<psz) return false; 
const IdxTy sd=ssz-psz; 
for(IdxTy i=0; i<sd; ++i)
{
if (strncmp(sp+i,pp,psz)==0) return true;

}
return false;
} // fast_grep

static bool fast_grep(const StrTy &s, const std::vector<StrTy> & v)
{
const IdxTy sz=v.size();
IdxTy i=0;
MM_ILOOP(j,sz){
const StrTy & p=v[j];
const char * pp=p.c_str();
const char * sp=s.c_str();
const IdxTy psz=p.length();
const IdxTy ssz=s.length();
if (ssz<psz) return false; 
const IdxTy sd=ssz-psz; 
for(; i<sd; ++i)
{
if (strncmp(sp+i,pp,psz)==0) break; //  return true;

}
return false;
} //j 

return true;
} // fast_grep





static my64int str_to_time(const ChTy * s)
{
my64int x;
std::stringstream ss;
ss<<s;
ss>>x;
return x; 
}

static my64int xxx_str_to_time(const ChTy * s)
{
my64int x;
//sscanf(s,"%11d",&x);
sscanf(s,"%lld",&x);
return x; 
}



static UlongTy date_to_long(const ChTy * s)
{
tm x;
::memset(&x,0,sizeof(x));
// these format chars seem critical, need some way to pick or fix LOL
sscanf(s, "%04d-%02d-%02d %02d:%02d:%02d",
       &x.tm_year, &x.tm_mon,  &x.tm_mday,
       &x.tm_hour, &x.tm_min, &x.tm_sec);
x.tm_year-=1900; x.tm_mon-=1;

// osx()<<MM_MARK<<" "<<s<<" "<<x.tm_min<<" "<<x.tm_sec<<CRLF;
return mktime(&x);
}

static UlongTy xxx_date_to_long(const ChTy * s)
{
tm x;
::memset(&x,0,sizeof(x));
// these format chars seem critical, need some way to pick or fix LOL
sscanf(s, "%04d-%02d-%02d %02d:%02d:%02d",
       &x.tm_year, &x.tm_mon,  &x.tm_mday,
       &x.tm_hour, &x.tm_min, &x.tm_sec);
x.tm_year-=1900; x.tm_mon-=1;

// osx()<<MM_MARK<<" "<<s<<" "<<x.tm_min<<" "<<x.tm_sec<<CRLF;
return mktime(&x);
}


static void dedicated_char_replace(char * line, IdxTy len)
{
const ChTy cf='_';
ChTy * ep=line+len;
ChTy * pc=line;
while ( pc<ep)
{
ChTy & c=*pc;
if ( c=='\'') c=cf;
else if (c==',') c=cf;
else if (c=='\"') c=cf;
// somewhat dubious if these remain in query strings? 
else if (c=='&') c=cf;
++pc;

}


}

static void dedicated_char_replace_better(char * line, IdxTy len)
{
const ChTy cf='_';
ChTy * ep=line+len;
ChTy * pc=line;
while ( pc<ep)
{
ChTy & c=*pc;
//if ( c=='\'') c=cf;
if (c==',') c=cf;
//else if (c=='\"') c=cf;
// somewhat dubious if these remain in query strings? 
else if (c=='&') c=cf;
++pc;

}
}
// inplace url decode note also html
template < bool skip_to_question, bool do_html > static void url_decode(ChTy * s)
{
//IdxTy i=0;
ChTy * d=s;
if ( skip_to_question)
{
while ((*s!=0)&&(*s!='?')) { ++s;} 
if ( *s!=0) ++s;
}
//if ( *s==0)  { *d=0; return; }
while ( *s!=0)
{
if ( *s=='%')
{
++s;
ChTy cnew=(ChTy)aotip(s);
*(d++)=cnew;
}
else if (do_html&&(*s=='&'))
{
const ChTy * r=(++s);
while ((*s!=0)&&(*s!=';')) { ++s;} 
IdxTy sx=s-r;
ChTy dum[sx+2];
::strncpy(dum,r,sx);
dum[sx]=0;
if (*s!=0) ++s;
// ok, for newer code lists this needs to be map but ok for now.
// however, order should be  based on P 
if  ( ::strcmp(dum,"amp") == 0) *(d++)='&';
else if  ( ::strcmp(dum,"lt") == 0) *(d++)='<';
else if  ( ::strcmp(dum,"gt") == 0) *(d++)='>';
else if  ( ::strcmp(dum,"quot") == 0) *(d++)='"';
else { d+=sx; } // go ahdead and write whoel thingout. 
}
else *(d++)=*(s++);
// ++s;
} // *s; 

*d=0;
return;
}

static int aotip(ChTy * & s)
{
IdxTy i=0;
ChTy c=*s;
if (::isdigit(c)) {i=(i<<4)+(c-(int)'0'); ++s; }
else if (( c>='a')&&(c<='z'))  { i=(i<<4)+(c-(int)'a'); ++s; }
else if (( c>='A')&&(c<='Z'))  { i=(i<<4)+(c-(int)'A'); ++s; }
if (::isdigit(c)) {i=(i<<4)+(c-(int)'0'); ++s; }
else if (( c>='a')&&(c<='z'))  { i=(i<<4)+(c-(int)'a'); ++s; }
else if (( c>='A')&&(c<='Z'))  { i=(i<<4)+(c-(int)'A'); ++s; }
return i;

}


// map decode of query string
typedef std::map<StrTy, StrTy> map_type;
// recognizes only 2 real things, "&" and "="
// note that s is destroyed in this process
template < bool skip_to_question ,bool lower_key, bool do_url_decode> 
	static void query_decode(map_type & d,  ChTy * s)
{
if ( skip_to_question)
{
while ((*s!=0)&&(*s!='?')) { ++s;} 
if ( *s!=0) ++s;
}
//if ( *s==0)  { *d=0; return; }
while ( *s != 0)
{
ChTy * r=s;
while ((*s!=0)&&(*s!='=')) { ++s;} 
// if we have '=' then we have a key here, 
bool have_more=(*s!=0);
*s=0;
if ( lower_key) in_place_to_lower(r);
if ( do_url_decode) url_decode<false,false>(r); 
StrTy rkey=(r);
if ( !have_more) {d[rkey]=""; return ; }
++s;
r=s;
while ((*s!=0)&&(*s!='&')) { ++s;} 
have_more=(*s!=0);
*s=0;
if ( do_url_decode) url_decode<false,false>(r); 
d[rkey]=StrTy(r);
if (have_more) ++s;
}




}



// dedicated path parser, really for static images but now generalized
// slightly
// path element is a null terminatedlist of field names for the path
// components, a bit slow but matches qs decode approach better

template < bool skip_to_question ,bool lower_key, bool do_url_decode> 
	static void path_decode(map_type & d,  ChTy * s,const ChTy ** pathelements)
{
if ( skip_to_question)
{
while ((*s!=0)&&(*s!='?')) { ++s;} 
if ( *s!=0) ++s;
}
//if ( *s==0)  { *d=0; return; }
while (( *s != 0)&&(*pathelements!=0))
{
ChTy * r=s;
while ((*s!=0)&&(*s!='/')) { ++s;} 
// if we have '=' then we have a key here, 
bool have_more=(*s!=0);
*s=0;
if ( lower_key) in_place_to_lower(r);
if ( do_url_decode) url_decode<false,false>(r); 
StrTy rkey=(r);

/*if ( !have_more) {d[rkey]=""; return ; }
++s;
r=s;
while ((*s!=0)&&(*s!='/')) { ++s;} 
have_more=(*s!=0);
*s=0;
if ( do_url_decode) url_decode<false,false>(r); 
*/


d[StrTy(*(pathelements++))]=rkey;



if (have_more) ++s;
}




}

static void make_single(StrTy & f)
{//  Ss
std::stringstream  ss;
//MM_ERR(MMPR(f))
ss<<f;
const StrTy forig=f;
f=StrTy("");
StrTy g=f;
// FUDD 
while (!ss.eof() && ss.good()) {g="";  ss>>g;
//if (!ss.good()) break;
//if (ss.eof()) break; 
if ( f!="" ) f=f+StrTy("_");
f=f+g;
 }
//MM_ERR(MMPR3(forig,f,g))
if (f==StrTy("")) {f="blank"; return ; }
//if (true ) return;
const IdxTy sz=f.length();
if (sz==0) { f="blank"; return; }
char c[sz+1];
const char * p=f.c_str();
IdxTy pc=0;
for(IdxTy i=0; i<sz; ++i)
{
char d=p[i];
if (d==' ') d='_';
else if (d=='\t') d='_';
else if (d==',') d='_';

c[pc]=d; // p[i];
++pc;
}
c[pc]=0;
f=&c[0];

}
template <class Tin, class Tout> 
static StrTy common_prefix(Tout &d, const Tin & s)
{
const IdxTy sz=s.size();
// don't bother with just 1 although maybe a flag for that 
if (sz<2) return StrTy();
IdxTy sxmatch=0;
IdxTy szmin=s[0].size();
MM_ILOOP(i,sz) { if (s[i].size()<szmin) szmin=s[i].size();}
MM_ILOOP(i,szmin)
{

IdxTy j=1;
while (j<sz)
{
if (s[0].c_str()[i]!=s[j].c_str()[i]) break;
++j;
} // j
if (j!=sz) break;
++sxmatch; // actually this equals i, could just break loop etc. 
} // i 
MM_ILOOP(j,sz) { d.push_back(StrTy(s[j].c_str()+sxmatch)); } // j 
return s[0].substr(0,sxmatch);
} // common_prefix

static StrTy  to_lower(const StrTy &  s)
{
const IdxTy sz=s.length();
const char * sc=s.c_str();
char c[sz+1];
c[sz]=0;
for(IdxTy i=0; i<sz; ++i) c[i]=::tolower(sc[i]); 
return StrTy(c);
}

static void in_place_to_lower(ChTy * s)
{
while (*s!=0)
{
*s=(ChTy)(::tolower(*s));
++s;

} // while


}
// requote and escape strings
template <bool quote=false, bool esc_sp=true,char zed_replace=0> static IdxTy re_escape(ChTy * d, const ChTy * s,const IdxTy _maxlen)
{
if ( _maxlen<7) throw ErTy("maxlen in re_escape must be larger than 7");
const IdxTy maxlen=_maxlen-2;
IdxTy w=0;
if ( quote)d[w++]='"';
const IdxTy slen=::strlen(s);
if (zed_replace!=0)
{
if (slen==0)
{
if (quote) {d[0]='"'; d[1]=zed_replace; d[2]='"'; d[3]=0; return 3; }
else { d[0]=zed_replace; d[1]=0; return 1; } 
} //slen

}

for ( IdxTy i=0; i<slen; ++i)
{
const ChTy & c=s[i];
if ( esc_sp) if ( c==' ') d[w++]='\\';
d[w++]=c;
if ( w>=maxlen) throw ErTy(" field length too long during re-secape");
}
if ( quote) d[w++]='"';

d[w]=0; return w;


}
// inplace adapter, could make this a template on fptr
static void de_delim(ChTy *d, const ChTy * s)
{
::strcpy(d,s);
de_delim(d);

}

// replace dumb intra-field chars with non-delims

static void de_delim(ChTy * s)
{
while ( *s != 0)
{
ChTy & c=*s;
// this needs to fix delims and things that sql could object
if ( c==' ') c='_';
else if ( c=='\t') c='_';
++s;
} //  s!=0


}



// Strings with containing quotes and escaped quotes etc
// can now be de-escaped in negative definine size case
template <bool noquote=false > static void de_escape(ChTy * s, const IdxTy len)
{
IdxTy w=0;
bool now_esc=false;
for (IdxTy i=0; i<len; ++i)
{
ChTy & c = s[i];
// not sure if this is pointer aliasing as oberved
// by "-O3" ....
// branching probs are esc is low P
if ( noquote) if (c=='"') continue;

if ( c!='\\') s[w++]=c;
else 
if ( now_esc) { s[w++]=c; now_esc=false; }
else now_esc=true; // I supposed esc at eos should
// be considered a error but esc should pick up space
// unless eol on input file. 

}
s[w]=0;

}

static StrTy trim(const StrTy & s)
{
typedef std::stringstream Ss;
Ss ss;
ss<<s;
StrTy x;
while ( true) { StrTy y; ss>>y; if ( y=="") break; if (x!="") x=x+StrTy(" "); x=x+y; }
return x;
}


static StrTy trim_crlf(const StrTy & s)
{
typedef std::stringstream Ss;
const char * p=s.c_str();
const IdxTy len=s.length();
Ss ss;
for(IdxTy i=0; i<len; ++i) { if (p[i]!='\r')
if (p[i]!='\n') ss<<p[i]; } 
return ss.str();
}

static StrTy trim_int(const StrTy & s)
{
typedef std::stringstream Ss;
const char * p=s.c_str();
const IdxTy len=s.length();
Ss ss;
for(IdxTy i=0; i<len; ++i) { if (p[i]>='0')
if (p[i]<='9') ss<<p[i]; } 
return ss.str();
}

// https://stackoverflow.com/questions/180947/base64-decode-snippet-in-c

static std::string base64_encode(const std::string &in) {

    std::string out;

    int val = 0, valb = -6;
    for (unsigned char c : in) {
        val = (val << 8) + c;
        valb += 8;
        while (valb >= 0) {
            out.push_back("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[(val>>valb)&0x3F]);
            valb -= 6;
        }
    }
    if (valb>-6) out.push_back("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[((val<<8)>>(valb+8))&0x3F]);
    while (out.size()%4) out.push_back('=');
    return out;
}

template <class Ty> 
static IdxTy  base64_decode(Ty & d, const std::string &in) {

//    std::string out;

    std::vector<int> T(256,-1);
    for (int i=0; i<64; i++) T["ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[i]] = i;

    int val=0, valb=-8;
    for (unsigned char c : in) {
        if (T[c] == -1) break;
        val = (val << 6) + T[c];
        valb += 6;
        if (valb >= 0) {
            d.push_back(char((val>>valb)&0xFF));
            valb -= 8;
        }
    }
    return 0;
}





typedef  StrTy (* EscPtr)( const ChTy *, const IdxTy);
// make output an unambiguous single field
// escape delims or quote
static StrTy add_nothing(const ChTy *s, const IdxTy len)
{ return StrTy(s); } 
template < bool quote=false, bool dedelim=true, ChTy delim=' '> 
static 	StrTy add_escapes(const ChTy * s, const IdxTy len)
{
ChTy d[len*2+3];


return StrTy(d); 
}




// in place destructive parse as delimeters are replaced with zed
// in theory it would make sense for escape removal here
// but lets make it sep for now.
// multidelim will accept blank fieds denoted by multiple
// delims otherwise they are treated as singe whitespace
template <bool have_esc, char delim=' ',bool multidelim=false> static void parse_fields(ChTy * line, const IdxTy len, ChTy ** fields, IdxTy & field, 
const IdxTy maxfields )
{  

// escape logic
bool in_quote=false;
bool is_esc=false;

IdxTy prior=0;
IdxTy pc=0;
//remove extra junk
if ( !multidelim) { while (line[pc]==delim) { ++pc; } }
prior=pc;

while (pc<len)
{
char & c=line[pc++];
// this crap needs to be removed from loop in sane files
// insert escape logic here.
// rule: 
if ( have_esc)
{
// this stupid test occures on EACH char now... 
// if (pc>1) if (line[pc-2]=='\') continue; // no matter what ignored
if ( is_esc) { is_esc=false; continue; } 
// not already escaped then make it so 
if ( c=='\\') {is_esc=true; continue; }
// non escaped
// arguably quotes and escapes could be removed by us, but mess...
if (c=='"') { in_quote=!in_quote; continue; }
if (in_quote) continue; // otherwise let norm logic handle it. 
}
if (c==delim)
{
c=0;
fields[field++]=line+prior;
if ( field>=maxfields) return; 
//remove extra junk
if ( !multidelim)  { while (line[pc]==delim) { ++pc; } } 
prior=pc;

} // c== space
// escape logic alt
// cprior=c;

} //pc
// get trailing thing
//I guess a trailing delim should be a new blankc fiel?
// this came up in wendy's add 
if ( line[prior]!=0) fields[field++]=line+prior;
else if ( multidelim) fields[field++]=line+prior;

} //field count






static bool  ends_in(const ChTy * s, const IdxTy lens, const ChTy * p, const IdxTy lenp)
{
if (lenp<lens) return false;
const int cmp=strcmp(s,p+lenp-lens);
return cmp==0;
}

static bool  ends_in(const ChTy * s,  const ChTy * p)
{
const IdxTy lens=strlen(s);
const IdxTy lenp=strlen(p);

if (lenp<lens) return false;
const int cmp=strcmp(s,p+lenp-lens);
return cmp==0;
}





// find noth from end of char
static int count_back(ChTy * s, const ChTy c, IdxTy n,const IdxTy len )
{
// do not include zero, if zero just want whole string anyway
// doh, DO include zero as we want to remove the dot...
if ( len<1) return -1; 
for ( int j=len-1; j>=0; --j)
{
if ( s[j]==c) { --n; if (n<=0) return j; }
}
return -1;

}





// truncate s to the first occurence of any char in cs
static int trunc_to(ChTy * s, const ChTy *cs, const IdxTy ncs)
{
IdxTy i=0;
while ( s[i]!=0)
{
for(IdxTy j=0; j<ncs; ++j) 
 if ( s[i]==cs[j]) { s[i]=0; return i; }

++i;
}
return -1;

}









static int find_to(ChTy * s, const ChTy *cs, const IdxTy ncs)
{
IdxTy i=0;
while ( s[i]!=0)
{
for(IdxTy j=0; j<ncs; ++j) 
 if ( s[i]==cs[j]) {  return i; }

++i;
}
return -1;

}

// see is s ends with cs
static int ends_with(const ChTy * s, const ChTy *cs)
{
const int ncs=::strlen(cs);
const int ns=::strlen(s);
if ( ncs>ns) return -1;
// this returns zero if true doh??? 
return strcmp(cs,s+(ns-ncs));
}

// like java, find first occurence of cs in s
//static int indexOf(ChTy * s, const ChTy *cs)
static int indexOf(const ChTy * s, const ChTy *cs)
{
const int ncs=::strlen(cs);
const int  lim=::strlen(s)-ncs+1;
// eh, I take it int R uint doesn't work with int<0???
for ( int i=0; i<lim; ++i)
{

int j=0;
for(; j<ncs; ++j) 
 if ( s[i+j]!=cs[j]) { break; }
if ( j==ncs) return i; 
}
return -1;

}
// lens>lencs
static int indexOf(const ChTy * s, const int lens, const ChTy *cs, const int lencs)
{
const int ncs=lencs;
const int  lim=lens-ncs+1;
// eh, I take it int R uint doesn't work with int<0???
for ( int i=0; i<lim; ++i)
{

int j=0;
for(; j<ncs; ++j) 
 if ( s[i+j]!=cs[j]) { break; }
if ( j==ncs) return i; 
}
return -1;

}

// delsize 
// indexOf returns offset if found or -1 if not.
// same here up to user to sort out size.. 
static int contains(const ChTy * s, const int lens, const ChTy *cs, const int lencs)
{
if (lens<lencs) return indexOf(cs,lencs,s,lens);
if (lens>lencs) return indexOf(s,lens,cs,lencs);
return indexOf(s,lens,cs,lencs);
} // contains 



// return the first char past the string or zero

static int endOf(const ChTy * s, const ChTy *cs)
{
const int ncs=::strlen(cs);
const int  lim=::strlen(s)-ncs+1;
// eh, I take it int R uint doesn't work with int<0???
for ( int i=0; i<lim; ++i)
{

int j=0;
for(; j<ncs; ++j) 
 if ( s[i+j]!=cs[j]) { break; }
if ( j==ncs) return i+j; 
}
return 0;

}





// should fail at compile time
//duphus_agg() {}
// this assumes that we can do IN PLACE MOD on f
static void bin_geo_coord(ChTy * f, ChTy * fnew)
{
// be safe for now... 
// ChTy fnew[::strlen(f)<<1]; 
//typedef float LocTy;
// doh, this is to copy a MISTAKE in the java code
// and bin to 1 degree NOT .1 
int loc=(int)(((int)::atof(f))*10);
//int loc=(int)(::atof(f)*10);
if ( loc<0) { itoa(-loc,&fnew[1],10);  fnew[0]='-'; } 
else itoa((int)loc,fnew,10);
/*
if ( ::strlen(fnew)>::strlen(f))
{
os(3)<<MM_MARK<<" "<<fnew<<" is longer than input location "<<f<<CRLF;
throw std::runtime_error(" location information is bad see stream 3");
}

else ::strcpy(f,fnew);
*/

}

static StrTy ext_app(const char * p, const char * s)
{

        StrTy fo=StrTy(p)+StrTy(s);
if ( PEDANTIC_AND_CONDESCENDING)
osx()<<MM_MARK<<" appending to give "<<fo<<CRLF;
//      const char * fout=fo.c_str();
//      return fout;
        return fo;
}

// this is really stupid since the m_powers implies radix
//char * itoa(IdxTy v, char * s,IdxTy * pxx)
static char * itoa(IdxTy v, char * s,const IdxTy radix )
{
IdxTy loc=0;
IdxTy pow=1;
IdxTy powlast=1;
// IdxTy powidx=0;
while ( v>=pow) {powlast=pow; pow=pow*radix; }
//while ( v>=pow) {powlast=pow; ++powidx; pow=pxx[powidx]; }

while ( powlast!=0)
{
IdxTy i=v/powlast;
if (i<10) s[loc++]=(char)('0'+i);
else s[loc++]=(char)('a'+i-10);
v-=i*powlast;
// this is slow, the results should be saved somewhee...
// use m_powers
powlast/=radix;
}
if ( loc==0) s[loc++]='0';
s[loc]=0;
return s;

}
// tried to reconstruct the radix array code, need
// to run against other code and chec 
static char * itoa(IdxTy v, char * s,const IdxTy * pxx )
{
IdxTy loc=0;
IdxTy pow=1;
IdxTy powlast=1;
int  powidx=0;
//while ( v>=pow) {powlast=pow; pow=pow*radix; }
while ( v>=pow) {powlast=pow; ++powidx; pow=pxx[powidx]; }

while ( powidx>=0 )
{
IdxTy i=v/powlast;
if (i<10) s[loc++]=(char)('0'+i);
else s[loc++]=(char)('a'+i-10);
v-=i*powlast;
// this is slow, the results should be saved somewhee...
// use m_powers
--powidx;
powlast=pxx[powidx];
}
if ( loc==0) s[loc++]='0';
s[loc]=0;
return s;

}

// I guess here you could normalize with one
// division and then make the rest all multiplies
// or maybe even exhaustive or binary search rather than div?
// as written this will probably accumulate errors in radix, 
// need to recompute each time. 
static char * ftoa(float v, char * s,const IdxTy radix,const IdxTy digs )
{
if ( radix!=10) throw ErTy(" stl no support non-dec fp radix and I too lazy");

typedef std::stringstream ConvTy;
ConvTy ss;
// stupid manipulators no work on fp...

ss<<std::setprecision(digs)<< v; // for now, ignore radix etc.
::strcpy(s,ss.str().c_str());
// needs a length check really
return s;

}


template <class Ty> static char * stoa(const Ty & v, char * s )
{

typedef std::stringstream ConvTy;
ConvTy ss;
// stupid manipulators no work on fp...
// maybe this not work for long int???gg
ss<< v; // for now, ignore radix etc.
::strcpy(s,ss.str().c_str());
// needs a length check really
return s;

}





static char * not_work_ftoa(float v, char * s,const IdxTy radix,const IdxTy digs )
{
IdxTy loc=0;
IdxTy pow=1;
IdxTy powlast=1;
float dv=radix^digs;
dv=1.0/dv;
 IdxTy powidx=0;
while ( v>=pow) {++powidx; powlast=pow; pow=radix^powidx; }
//while ( v>=pow) {powlast=pow; ++powidx; pow=pxx[powidx]; }

while ( powidx>=(-digs))
//while ( powlast>=(dv))
{
IdxTy i=v/powlast;
if (i<10) s[loc++]=(char)('0'+i);
else s[loc++]=(char)('a'+i-10);
v-=i*powlast;
// this is slow, the results should be saved somewhee...
// use m_powers
//powlast/=radix;
--powidx;
powlast=radix^powidx;
}
if ( loc==0) s[loc++]='0';
s[loc]=0;
return s;

}



static int matches(const char * tgt, const char * test)
{
IdxTy tgtlen=strlen(tgt);
IdxTy testlen=strlen(test);
// sob- eh this is signed????

int poe=testlen-tgtlen;
// osx()<<MM_MARK<<" compare "<<poe<<" "<<tgtlen<<" "<<testlen<<" "<<tgt<<" - vs -  "<<test<<CRLF;
// we need to find all of tgt somewhere in test 
//return strcmp(tgt,test);

if( poe<0) return 1;
for ( int i=0; i<=poe; ++i)
{
if ( strncmp(tgt,test+i,tgtlen)==0) return 0;
}
return -1;

}

static StrTy exxx_t_app(const char * p, const char * s)
{

	StrTy fo=StrTy(p)+StrTy(s);	
if ( PEDANTIC_AND_CONDESCENDING) 
osx()<<MM_MARK<<" appending to give "<<fo<<CRLF; 
//	const char * fout=fo.c_str();
//	return fout;
	return fo;
}


static StrTy extraneous_ws(const StrTy & s)
{
//std::stringstream ss;
//Ss ss;
//ss<<s;
const IdxTy slen=s.length();
if (slen==0) return s; 
//while ( true) { StrTy y; ss>>y; if ( y=="") break; if (x!="") if ( y!=" ") x=x+StrTy(" "); x=x+y; } 
char * c=new char[slen+1];
memcpy(c,s.c_str(),slen);
c[slen]=0;
IdxTy start=0;
IdxTy send=slen-1;
while (c[start]==' ') ++start;
// make sure the send loop can terminate... 
if (start==slen) return StrTy();
while (c[send]==' ') --send;
c[send+1]=0;

StrTy x(c+start);
delete[]  c;
return x;
}

static StrTy ws_remove(const StrTy & s, IdxTy n=~0 )
{
IdxTy cnt=0;
const char * p=s.c_str();
const IdxTy sz=s.length();
char c[sz+1];
//int  level=0;
//IdxTy nsofar=0;
for (IdxTy i=0; i<sz; ++i) {  char  b=p[i];
//if (nsofar>=n) { c[cnt]=b; ++cnt; continue; } 
//if ( b=='[') { ++level; continue;}
//if ( b==']') {++nsofar;  --level; continue;}
//if ( level==0) { c[cnt]=b; ++cnt; }
if (b==' ') b='_';
c[cnt]=b; ++cnt;
}
c[cnt]=0;
return StrTy(c);
}






// BraceRemove 
static StrTy brace_remove(const StrTy & s, IdxTy n=~0 )
{
IdxTy cnt=0;
const char * p=s.c_str();
const IdxTy sz=s.length();
char c[sz+1];
int  level=0;
IdxTy nsofar=0;
for (IdxTy i=0; i<sz; ++i) { const char  b=p[i];
if (nsofar>=n) { c[cnt]=b; ++cnt; continue; } 
if ( b=='[') { ++level; continue;}
if ( b==']') {++nsofar;  --level; continue;}
if ( level==0) { c[cnt]=b; ++cnt; }
}
c[cnt]=0;
return StrTy(c);
}

static StrTy first_brace(const StrTy & s)
{
IdxTy cnt=0;
const char * p=s.c_str();
const IdxTy sz=s.length();
char c[sz+1];
int  level=0;
for (IdxTy i=0; i<sz; ++i) { const char & b=p[i];
if ( b=='[') { ++level; continue;}
// return first completed one and crry on with bad syntax 
if ( b==']') { --level;if ( level==0) break;  continue;}
if ( level>0) { c[cnt]=b; ++cnt; }
}
c[cnt]=0;
return StrTy(c);
}




static IdxTy vatou(const char * g, const IdxTy m )
{

	bool ok=true;
	if ( ::strlen(g)==0) ok=false;
	if (ok)
	{	const IdxTy lj=::strlen(g);
		for ( IdxTy j=0; j<lj; ++j)
			if ( !::isdigit(g[j])) { ok=false ; break; }
	}
	if ( !ok)
	{ 
	esx()<<MM_MARK<<" bad os number "<<g<<" for output file "<<" "<< CRLF;
	throw CmdEr(ext_app(" file spec not good ",g));
	}
// this relies on unsigned being >=0
	IdxTy idx= ::atoi(g);
if ( m==0) return idx;
	if (idx>m) 
	{
		esx()<<MM_MARK<<" bad os number "<<g<<" parsed to "<<idx
		<<" limit is "<<m<<CRLF;
		throw CmdEr(ext_app(" bad parse of input number  ",g));
	}
return idx;

}

#if 0 
// written fro mjm_snacks for variable subs in macro defs
// not used yet not sure what to do 
static StrTy replace(const char * s, const IdxTy first, const IdxTy last, const char * r)
{
this is fudded up 
const IdxTy fudd1=strlen(s);
const IdxTy fudd2=strlen(r);
const IdxTy flen=fudd1+fudd2-(last-first);
char c[fudd1+fudd2+1];
memcpy(c,s,first);
memcpy(c+first,r,last-first+1);
memcpy(c,s+last,fudd1-last);
c[d]=0;
return StrTy(c);

}
#endif


};  //myt


#endif

