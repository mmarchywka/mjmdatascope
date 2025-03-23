#ifndef MJM_INSTRUMENTS_H__
#define MJM_INSTRUMENTS_H__
#include <mjm_globals.h>
#include <mjm_read_buffer.h>
#include <map>
#include <algorithm>
#include <string>
#include <ios>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>
#include <limits>
#include <stdlib.h> // strtol
#include "mjm_format_text.h"
#include "mjm_misc_parse.h"
// add for allowing program variables to be passed easily 
// FUDD cirular fcking deifniiotns need another fck 
//#include "mjm_logic_base.h"

// https://stackoverflow.com/questions/13168833/c-command-line-interface
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>




// assured to compile away lol 
class mjm_dumb
{
public:
template <class Ty, class Tw, class Tz=std::string> static void OK(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{}
template <class Ty, class Tw, class Tz=std::string> 
		static void OKV(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{}

template <class Ty, class Tw, class Tz=std::string> static bool  nan(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{return false; }

template <class Ty, class Tw, class Tz=std::string, class Tl=std::string> static bool  denan(const  Ty & val , const Tw & warn, const Tz & warn2="", const Tl & warn3="") 
{return false; }
template <class Ty, class Tw, class Tz=std::string, class Tl=std::string> static bool  denanf(Ty & val , const Tw & warn, const Tz & warn2="", const Tl & warn3="") 
{return false; }

}; // mjm_dumb
// needs mjm_globals and cmath limits 
// http://stackoverflow.com/questions/2487653/avoiding-denormal-values-in-c
class mjm_numb
{

public:

typedef int CtTy;

static CtTy & exit_limit() { static CtTy x=1000; return x; }
static CtTy & count() { static CtTy x=0; return x; }
static bool suppress() { return count()>exit_limit(); } 
static CtTy & inconly() { ++count(); 
if (count()==exit_limit()) MM_ERR(" denorm warning limit reached, future errors supporessed") 

 return count(); } 
static CtTy & inc() { ++count(); 
if (count()>exit_limit())
{
MM_ERR(" exit on too many fp errors "<<count())
std::exit(666);
}
 return count(); }

template <class Ty> static std::string ntoh(const Ty & v )
{
std::stringstream ss;
const unsigned int sz=sizeof(Ty);
unsigned char * p= (unsigned char *)&v;
typedef unsigned int PrTy;
// this is ok for string, not numbers lol 
//for (unsigned int i=0; i<sz; ++i) {PrTy c=PrTy(*(p+i)); if( c<16) ss<<"0"; ss<<std::hex<<c; if ((i&3)==3) ss<<" "; }
for (unsigned int i=0; i<sz; ++i) {PrTy c=PrTy(*(p+(sz-1-i))); if( c<16) ss<<"0"; ss<<std::hex<<c; if ((i&3)==3) ss<<" "; }

return ss.str();

}

template <class Ty, class Tw, class Tz=std::string> static void OK(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{
if ( std::fpclassify( val ) == FP_SUBNORMAL )
{
MM_ERR( "subnormal "<<val<<" "<<warn<<" "<<warn2 )
}

}

template <class Ty, class Tw, class Tz=std::string> 
		static void OKV(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{
if ( val != 0 && ::fabsf( val ) < std::numeric_limits<Ty>::min() )
{
MM_ERR( "subnormal "<<val<<" "<<warn<<" "<<warn2 )

}

}
template <class Ty, class Tw, class Tz=std::string> 
		static void OKC(const  Ty & col , const unsigned int sz, const Tw & warn, const Tz & warn2="") 
{

for (unsigned int i=0; i<sz; ++i)
{
const auto & val=col[i];
if ( val != 0 && ::fabsf( val ) < std::numeric_limits<Ty>::min() )
{
MM_ERR( "subnormal "<<val<<" ["<<i<<"] "<<warn<<" "<<warn2 )
}
if ( isnan(val) )
{
MM_ERR( "nan "<<val<<" "<<ntoh(val)<<" idx="<<i<<" "<<warn<<" "<<warn2 )

}

}

}

template <class Ty, class Tw, class Tz=std::string> static bool  nan(const  Ty & val , const Tw & warn, const Tz & warn2="") 
{
if ( isnan(val) )
{
MM_ERR( "nan "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2 )
return true;
}
return false;
}


template <class Ty, class Tw, class Tz=std::string, class Tl=std::string> static bool  denan(const  Ty & val , const Tw & warn, const Tz & warn2="", const Tl & warn3="") 
{
if ( isnan(val) )
{
if (!suppress()) {  std::cerr<< "mjm_instruments "<<__LINE__<<"  nan "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2<<":"<<warn3 <<CRLF;  std::cerr.flush();
}
inc();

return true;
}
if ( std::fpclassify( val ) == FP_SUBNORMAL )
{
//MM_ERR( "subnormal "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2<<" "<<warn3 )
if (!suppress()) {  std::cerr<< "mjm_instruments "<<__LINE__<<"  sub "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2<<":"<<warn3 <<CRLF;  std::cerr.flush(); }
inconly();
return true;
}
return false;
}
template <class Ty > static bool dnz(Ty & val)
{
if ( std::fpclassify( val ) != FP_SUBNORMAL ) return false;
val=0;
return true; 

}

template <class Ty, class Tw, class Tz=std::string, class Tl=std::string> static bool  denanf(Ty & val , const Tw & warn, const Tz & warn2="", const Tl & warn3="") 
{
if ( isnan(val) )
{
MM_ERR( "nan "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2 )
val=0;
return true;
}
if ( std::fpclassify( val ) == FP_SUBNORMAL )
{
MM_ERR( "subnormal "<<val<<" "<<ntoh(val)<<" "<<warn<<" "<<warn2<<" "<<warn3 )
val=0;
return true;
}
return false;
}





}; // mjm_numb


class LineIterator
{

class TrDef{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::ostream OsTy;
typedef std::istream IsTy;
typedef std::ofstream Ofs;
typedef std::stringstream Ss;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef mjm_read_buffer<TrDef> Oflow;

protected:
typedef std::istream IsTy; 
typedef std::ifstream FisTy; 
typedef std::stringstream Ss; 
typedef unsigned int IdxTy;

typedef IdxTy Iv;
typedef mjm_var_parse<TrDef,char,Iv> CharTable;
public:
typedef std::vector<StrTy> Words;
enum { CODE_SPLIT_BASE=100};
enum { MAX_CODE_SPLIT=CODE_SPLIT_BASE+10};
enum { KVP_SPLIT_BASE=MAX_CODE_SPLIT+1};
enum { MAX_KVP_SPLIT=KVP_SPLIT_BASE+256};
enum { IDIOM_SPLIT_BASE=MAX_KVP_SPLIT+1};
StrTy SPLITS() const
{
Ss ss;
ss<<MMPR3(IDIOM_SPLIT_BASE,KVP_SPLIT_BASE,CODE_SPLIT_BASE);
return ss.str();
}
LineIterator():  m_f(0), m_close(false), m_serial(0),m_offset(0),m_sep(' '),m_splitter(0),m_prompt("mjm>"),m_disable_readline(false),m_save_last_line(false),m_use_alt_eol(false),
 m_remove_trailing_crlf(false),m_remove_white(false),
m_eol(0),m_eol_flags(0) {}
LineIterator(IsTy * f) : m_f(f), m_close(false), m_serial(0) ,m_offset(0), m_sep(' '), m_splitter(0), m_prompt("mjm>"),m_disable_readline(false),m_save_last_line(false),m_use_alt_eol(false),
 m_remove_trailing_crlf(false), m_remove_white(false),
m_eol(0),m_eol_flags(0) {}
LineIterator(const StrTy & nm) 
	: m_f(new FisTy(nm)), m_close(true), m_serial(0) , m_offset(0), m_sep(' '), m_splitter(0),m_prompt("mjm>") ,m_disable_readline(false),m_save_last_line(false),m_use_alt_eol(false),
 m_remove_trailing_crlf(false), m_remove_white(false),
m_eol(0),m_eol_flags(0)
{
//m_f= new FisTy(nm);
}

virtual ~LineIterator() { //if (m_close) { m_f->close();
clear_stack(); release(); } // 	if (m_close)	 delete m_f; }   
void source(const StrTy & fn,const bool die)
{
Source(fn,die); 

}
void set_prompt(const StrTy & s) { m_prompt=s; }
void set_prompt(const char *  s) { m_prompt=StrTy(s); }
IdxTy code_split_base() const { return CODE_SPLIT_BASE; } 
void set_split(const IdxTy & s, const char sep)
{
m_splitter=s;
m_sep=sep; 
}

const bool ok()const  {
// this could retrn false is stack empty
while (!ok(m_f)) { if (!pop_is()) break; } 
// may be true althogh new file should be checked 
return ok(m_f);  
} 
const bool good()const  {
// this could retrn false is stack empty
while (!good(m_f)) { if (!pop_is()) break; } 
// may be true althogh new file should be checked 
return good(m_f);  
} 

const bool ok(const IsTy * f) const {if (f==0) return false;  return  (f->good()&&!f->eof()) ; } 
const bool good(const IsTy * f) const {if (f==0) return false;  return  (f->good()) ; } 

const IdxTy size() const { return m_words.size()-m_offset; } 
const IdxTy size0() const { return m_words.size(); } 
const Words & words() const { return m_words; }
const StrTy & word(const IdxTy i ) const { return m_words[i+m_offset]; }
const double wordf(const IdxTy i ) const { return ::atof(m_words[i+m_offset].c_str()); }
const int  wordi(const IdxTy i ) const { return ::atoi(m_words[i+m_offset].c_str()); }
const unsigned int  wordu(const IdxTy i ) const { return ::atoi(m_words[i+m_offset].c_str()); }
// for hieratchail commands, just remove the current ones
//template <class Ty> Ty convert(const StrTy & x) const { return Ty(x); } 
//ZZtemplate <> 
void  convert(int & z, const StrTy & x) const { z=::atoi(x.c_str()); } 
void  convert(IdxTy & z,const  StrTy & x) const { z=::atoi(x.c_str()); } 
void  convert(double & z, const StrTy & x) const { z=::atof(x.c_str()); } 

template <class Ty>
bool get( Ty & d, const StrTy & nm) const
{
const IdxTy sz=size();
if (sz<2) return false; 
for (IdxTy i=0; i<(sz-1); ++i) if (word(i)==nm) { convert(d,word(i+1)); return true; } 

return false; 
}

void push()
{
// this can still get fed up if new line is obtained
if (m_offset>=size0()) return; 
m_offset+=1;
}
void pop()
{
if (m_offset>0) --m_offset;

}
bool nextok()
{
if (m_save_last_line){ m_save_last_line=false;  return true;}  // 2019-05-08
if (!ok()) return false;
next();
if (!true) return true; 
const bool asfck= good()||(m_line.length()!=0); 
// 2020-02-11 compiler complained irl apaprently due to commente dout below 
//const bool irl= (is_cin()&&readline_ok()) ; // ( m_f==&std::cin)
//MM_ONCE(" blanks doh "<<MMPR3(ok(),asfck,irl),)
//if ( irl) return true; 
//MM_ONCE(" blanks doh ",)
return asfck;
}

bool nextoks()
{
if (!ok()) return false;
if (m_save_last_line) m_save_last_line=false;
else next();
return true; 
}



bool is_cin() const { return  ( m_f==&std::cin); } 
bool readline_ok() const { return  ( !m_disable_readline); } 
void readline_ok(const bool x) { m_disable_readline=!x; } 
void save_last_line(const bool x) { m_save_last_line=x; } 
void use_alt_eol(const char c, const bool x=true) {m_eol=c;  m_use_alt_eol=x; } 
void cr_or_lf( const bool x=true) {if (x) m_eol_flags|=1; else m_eol_flags &=~1; } 
void remove_crlf(const bool b) {  m_remove_trailing_crlf=b; } 
void remove_white(const bool b) {  m_remove_white=b; } 
const StrTy & line() const { return m_line; } 
virtual void next()
{
if (is_cin()&&readline_ok()) // ( m_f==&std::cin)
{
	char * line = readline(m_prompt.c_str());
// this should allow for good() or !eof()  to return false preventing seg fault doh 
// even with a null return this says "good" wth????? 
    if(!line) {m_line="";  
if (m_f==&std::cin) m_f=0; 
//MM_ERR(" readline returns null danger will robinson "<<good()<<good(&std::cin)) 
return;   
} //break;
  remove_junk(line);
    if(*line) add_history(line);
    /* Do something with the line here */
	nextline(line);
    free(line);
	return;
}
enum { LLEN=1<<16, RWLEN=LLEN-3};
// 2021-02-16 stack overflow problem
//char line[LLEN];
char*  line = new char[LLEN];
// 
char * lp=&line[0];
Oflow * po=0;
while (true)
{
bool toolong=false;
if ((m_eol_flags&1)!=0)
{
IdxTy i=0;
for(; i<(RWLEN-1); ++i)
{
if (m_f->eof()) break;
if (!m_f->good()) break;
int c=m_f->get();
if (!m_f->good()) break;
if ((c=='\r')||(c=='\n')) break;
lp[i]=c;
}
if (i==(RWLEN-1)) toolong=true;
lp[i]=0;
}// cr or lf 
else if (m_use_alt_eol) m_f->getline(lp,RWLEN,m_eol);
else m_f->getline(lp,RWLEN);
if ( toolong ||(!(*m_f).eof() &&!(*m_f).good()))
{
MM_ERR(" line may bt too long not  truncating file, expand "<<MMPR(lp))
(*m_f).clear();
if ( po==0 ) po= new Oflow(20);
//(*po).append(lp,RWLEN-1);
(*po).append(lp);
continue; // skip the break;
}
else if (po) (*po).append(lp);
break;
} // true
if (po!=0) {(*po).cap();  lp=(*po).ptr(); } 
 remove_junk(lp);
//nextline(&line[0]);
nextline(lp);
delete  po;
delete[] line;
}
// TODO this should just look at raw line 
StrTy reassemble(const IdxTy i, const IdxTy f,const StrTy & sepnew=" ")
{   
StrTy din="" ; 
StrTy sep="";

StrTy sepx=sepnew; // " "; // FUDD StrTy(m_sep);
//for(IdxTy j=i; j<f; ++j) din=din+StrTy(" ")+li.word(j); from=din; 
for(IdxTy j=i; j<f; ++j) {din=din+sep+word(j);  sep=sepx; } 
return din;
//from=din; 
}


void remove_junk(char * line)
{
if (  m_remove_trailing_crlf ) 
{
IdxTy i=strlen(line);
while (i!=0)
{
--i;
bool crlf= ( line[i]=='\r')||(line[i]=='\n' ) ;
if (!crlf) break;
line[i]=0;

} // i!=0

} // if 
if ( m_remove_white)
{
IdxTy len=strlen(line);
IdxTy j=0;
bool f=false;
for(IdxTy i=0; i<=len; ++i)
{
const char ci=line[i];
if ( !f) { if ((ci!=' ')&& (ci!='\t')) f=true; } 
line[i-j]=ci; // line[i];
if (!f) ++j;
//if (!f) MM_ERR(MMPR3(i,j,line))
}

}
} // remove_junk

// use could call this and insert own lines without stream 
template <class Tline> void nextline ( const Tline & line) 
{
//Ss ss; ss<<line; StrTy n;
m_words.clear();
// can add quotes and comments and crap later 
// danger will robinson this does not update n if there is nothing there, need to check
//while (ok(&ss)) { ss>>n; m_words.push_back(n); }
//MM_ERR(" ASSFUCK "<<MMPR(m_splitter))
switch (m_splitter)
{
case 1: { parse_full(m_words,line,m_sep); break; } 
case 2: { parse_pdf_columns(m_words,line,m_sep); break; } 
case 3: { parse_full(m_words,line,m_sep,1); break; }  // preserve quotes
// note this is fullM apparently for mail :
case 4: { parse_fullm(m_words,line,m_sep,1); break; }  //sep not used  
case 5: { parse_fullm(m_words,line,m_sep,7); break; }  //sep not used  include ;  
// take consecutive delimeters as one. 
case 6: { parse_full(m_words,line,m_sep,2); break; } 
//void tokenize(Tv & words, const char * s, const IdxTy flags=0 )
// IIRC, break up into logical char groups that will often
// give deisred result or easily re-created. 
case 7: { tokenize(m_words,line,64); break; } 
case 8: { tokenize(m_words,line,0); break; } 
// just pass esc and quote etc 
case 9: { tokenize(m_words,line,25); break; } 
// IIRC, originally intended for TooBib, a few key words followed
// by literal unescaped code 
case CODE_SPLIT_BASE: 
case CODE_SPLIT_BASE+1: 
case CODE_SPLIT_BASE+2: 
case CODE_SPLIT_BASE+3: 
case CODE_SPLIT_BASE+4: 
case CODE_SPLIT_BASE+5: 
case CODE_SPLIT_BASE+6: 
case CODE_SPLIT_BASE+7: 
{ parse_partial(m_words,line,m_sep,  m_splitter-CODE_SPLIT_BASE, 0 ); break; } 
// should have = and then nothing parsed on rhs. 
case KVP_SPLIT_BASE: 
{ parse_kvp(m_words,line,m_sep,  m_splitter-KVP_SPLIT_BASE, 0 ); break; } 
// idiosncratic partial parses 
case IDIOM_SPLIT_BASE: 
{ parse_idiom(m_words,line,m_sep,  m_splitter-KVP_SPLIT_BASE, 0 ); break; } 


default:
{                
parse_sep(line);
//Ss ss; ss<<line; StrTy n;
//   while (ok(&ss)) {n="";  ss>>n; m_words.push_back(n); } } 
}
}

++m_serial; 
m_line=StrTy(line);
}



template <class Tline> void parse_sep( const Tline & line)
{
const IdxTy sz=strlen(line); // line.length();
const char * s=line; // .c_str();
char c[sz+1];
IdxTy p=0;
IdxTy q=0;
while (p<=sz)
{
c[p]=s[p];
if (c[p]==m_sep) {c[p]=0;}
if (c[p]==0) {  m_words.push_back(StrTy(c+q)); ++p; q=p; continue; } 
++p; 
}

}



bool Bit(const IdxTy x,const IdxTy b) { return ((1<<b)&x)!=0; } 
bool AllMask(const IdxTy x,const IdxTy m) { return ((m&x)==m); } 
bool AnyMask(const IdxTy x,const IdxTy m) { return ((m&x)!=0); } 

template <class Tv> 
void parse_full(Tv & words, const char * s,const char sep, const IdxTy flags=0 )
{
const IdxTy sz=strlen(s);
bool preserve_quotes=Bit(flags,0);
bool deb=false;
//MM_ERR(MMPR3(sz,sep,s))
char c[2*sz+1];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
bool esc=false;
bool quote=false;
const bool no_concat_seps=!Bit(flags,1);
//MM_ERR(" parser "<<MMPR4(no_concat_seps,preserve_quotes,flags,sep))
while (s[pc] != 0 )
{
const char & cpc=s[pc];
// do not pass the quote to esc to dest 
if (!esc) if (cpc=='\\')  { esc=true; ++pc; continue; }  // else esc=false;
if (!esc) if (cpc=='"') {deb=!true; quote=!quote; ++pc ; if (preserve_quotes){  c[dp]=cpc; ++dp;} continue;  } 
const bool br=(cpc==sep)&&(!esc)&&(!quote);
if (br){  c[dp]=0; 
if ((c[start]!=0)||no_concat_seps ) words.push_back(c+start); 
++dp; start=dp;    } else { c[dp]=cpc; ++dp; } 
++pc;
esc=false; // fck 2018-09-04 FUDD 
}
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

if (start!=dp){ c[dp]=0;  words.push_back(c+start);  } 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

}//  parse_full

//////////////////////////////////////////////

// parse only the first few words then pass remainder literally. Useful
// for code with a lot of quotes and crp
// pwords can be zero and entire line is returned as one word.  
template <class Tv> 
void parse_partial(Tv & words, const char * s,const char sep, const IdxTy pwords, const IdxTy flags=0 )
{
const IdxTy sz=strlen(s);
bool preserve_quotes=Bit(flags,0);
bool deb=false;
MM_ERR(MMPR3(sz,sep,s))
const IdxTy mxsz=10*sz+1;
char c[mxsz];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
bool esc=false;
bool quote=false;
const bool no_concat_seps=!Bit(flags,1);
IdxTy nwords=0;
//MM_ERR(" parser "<<MMPR4(no_concat_seps,preserve_quotes,flags,sep))
if ( pwords==0) { words.push_back(s); return; } 
while (s[pc] != 0 )
{
const char & cpc=s[pc];
// do not pass the quote to esc to dest 
if (!esc) if (cpc=='\\')  { esc=true; ++pc; continue; }  // else esc=false;
if (!esc) if (cpc=='"') {deb=!true; quote=!quote; ++pc ; if (preserve_quotes){  c[dp]=cpc; ++dp;} continue;  } 
const bool br=(cpc==sep)&&(!esc)&&(!quote);
if (br){  c[dp]=0; 
if ((c[start]!=0)||no_concat_seps ){  words.push_back(c+start);  ++nwords; } 
++dp; start=dp;    
if ( pwords<=nwords) { words.push_back(s+pc+1); return; } 
} else { c[dp]=cpc; ++dp; } 
++pc;
esc=false; // fck 2018-09-04 FUDD 
}
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

if (start!=dp){ c[dp]=0;  words.push_back(c+start);  ++nwords; } 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

}//  parse_partial
/////////////////////////////////////////
// parse_idiom : logical units like words or <foo> bar </foo>
template <class Tv> 
IdxTy parse_idiom(Tv & words, const char * s,const char sep, const IdxTy param, const IdxTy flags=0 )
{
IdxTy rc=0;
//const IdxTy sz=strlen(s);
switch (param)
{

default:
return parse_idiom_html(words,s,sep,param,flags); 
} // switch param 
return rc;

} // parse_idiom

class _parse_idiom_reg
{
public:
_parse_idiom_reg() { Init(); }
void inc() { ++pc; } 
void Init() { pc=dp=sp=0; 
 esc=quote=found_value=bad_hex_value=false;
} 
IdxTy pc,dp,sp;
bool esc; bool quote; bool found_value; bool bad_hex_value;

}; // _parse_idiom_reg
typedef _parse_idiom_reg parse_idiom_reg;
// now need to know output encoding however... 
//HTML Entity (decimal) &#8308;
//HTML Entity (hex) &#x2074;
// this is only PARSING, not fixing... 
template <class Tv> 
IdxTy parse_idiom_html(Tv & w, const char * s,const char sep, const IdxTy param, const IdxTy flags=0 )
{
IdxTy rc=0;
//IdxTy quote=0;
IdxTy gtlt=0;
IdxTy amp=0;
IdxTy pct=0;
// probably save a lot of time just passing this... 
const IdxTy sz=strlen(s);
const IdxTy mxsz=10*sz+1;
char *c = new char[mxsz];
parse_idiom_reg r;
bool esc=false;
bool escn=false;
while (s[r.pc] != 0 )
{
const char cc=s[r.pc];
esc=escn;
escn=false;
if (esc) {}
if (cc=='\\')
{
escn=true;
++r.pc; // 2022-07-06 this needs to be in other places doh ... 
continue;
} 
if (amp)
{
 if (cc==';'){ 
 c[r.dp]=0; if (r.dp){ w.push_back(StrTy(c)); r.dp=0; }
w.push_back(";"); 
--amp;
++r.pc; // 2022-07-06 this needs to be in other places doh ... 
continue;
}
}
if (pct)
{
c[r.dp]=cc; ++r.dp;
 if (r.dp==2){ --pct; 
 c[r.dp]=0; if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
continue;
}// writte thecomplete thing out 
}// pct 
else if (cc=='%')
{
// TODO needs to deal with unicode etc 
//IdxTy v=0;
 c[r.dp]=0; if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
++pct;
//c[r.dp]=char(v&255);
continue;
}
else if (cc=='&')
{ // just ignore for now...  PARSING 
c[r.dp]=0; if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
w.push_back(StrTy("&")); ++r.pc;
++amp;
continue;
}
else if (cc==' ') 
{
c[r.dp]=0; if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
w.push_back(StrTy(" "));
++r.pc;
continue;
} 
else if (cc=='<') 
{
++gtlt;
c[r.dp]=0;
if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
}

c[r.dp]=cc; ++r.dp;
if (cc=='>') 
{
--gtlt;
c[r.dp]=0;
if (r.dp){ w.push_back(StrTy(c)); r.dp=0; } 
}
r.inc();
} // s
c[r.dp]=0;
if (r.dp) w.push_back(StrTy(c));
delete [] c; 
return rc;
} // parse_idiom



// IIRC lines are for hierarchial key with "=" and then a value 
// return value matters ... 
template <class Tv> 
IdxTy parse_kvp(Tv & words, const char * s,const char sep, const IdxTy param, const IdxTy flags=0 )
{
const IdxTy sz=strlen(s);
bool preserve_quotes=Bit(flags,0);
const bool no_concat_seps=!Bit(flags,1);
const bool silent_err=Bit(flags,2);
const bool add_null_value=!Bit(flags,3);
bool deb=false;
//MM_ERR(MMPR3(sz,sep,s))
const IdxTy mxsz=10*sz+1;
// stack space crap again... 
// char c[mxsz];
char *c = new char[mxsz];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
bool esc=false;
bool quote=false;
bool found_value=false;
bool bad_hex_value=false;
//MM_ERR(" parser "<<MMPR4(no_concat_seps,preserve_quotes,flags,sep))
while (s[pc] != 0 )
{
	//typedef char Ch;
	const char & cpc=s[pc];
	// do not pass the quote to esc to dest 
	if (!esc) if (cpc=='\\')  { 
		esc=true; ++pc; continue; 
		}  // else esc=false;
	if (!esc) if (cpc=='"') {
		deb=!true; quote=!quote; ++pc ; 
		if (preserve_quotes){  
			c[dp]=cpc; ++dp;
			if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
			} 
		continue;  
		}
	if (false) if (!esc) if (!quote) if (cpc=='=') {
		c[dp]=0; 
		if ((c[start]!=0)||no_concat_seps ) words.push_back(c+start); 
			if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
		++dp; start=dp;    
		// anything remaining is the value
		words.push_back(c+start); 
		found_value=true;
		break;
		} // have = rest of line is value  
	if (esc){
		if (cpc=='x') { // hex value, may also use u and o. 
			if (dp>=(mxsz-2)) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
			IdxTy xs=numb( c, dp, s, pc, 2);
			if (xs!=0) { esc=false; break; } 
			esc=false;
			continue;
			} // cpc==x
		} // handling escp stuff 
	const bool br=(cpc==sep)&&(!esc)&&(!quote);
	if (br){  
		c[dp]=0; 
		if ((c[start]!=0)||no_concat_seps )
{	 const StrTy wnew=StrTy(c+start); 
	if (wnew!="=") 	words.push_back(c+start); 
	else
	{
		++dp; start=dp;   ++pc; 
			if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
		while (s[dp]!=0) { c[pc]=s[dp]; ++dp; ++pc;

			if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
 } 
		c[pc]=0;
		words.push_back(c+start); 
found_value=true;
		++dp; start=dp;
		break;
	} // else  equals 
} // word 
		++dp; start=dp;    
			if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
		} 
	else 
		{ 
		if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
c[dp]=cpc; ++dp; } 
		if (dp>=mxsz) { MM_ERR(" size too bib "<<MMPR3(pc,dp,mxsz)) } 
	++pc;
	esc=false; // fck 2018-09-04 FUDD 
} // while
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 
// could be an equal left dangling...
if (start!=dp){ c[dp]=0;  
	 const StrTy wnew=StrTy(c+start); 

if (wnew != "=") words.push_back(c+start);  
else
{
words.push_back(StrTy());
found_value=true;
}
} 

delete [] c; 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 
if (esc||quote||!found_value||bad_hex_value) {
	if (!silent_err) {  MM_ERR(" parse errors "<<MMPR4(esc,quote,found_value,bad_hex_value)<<MMPR(s))  } 
	if (!found_value) {
		if(add_null_value) words.push_back(StrTy()); 
		 return ~0; 
	}
	return 1; 
	} 
return 0;
}//  parse_kvp

IdxTy numb( char * c, IdxTy & dp, const char * s, IdxTy &pc, const IdxTy n)
{
	typedef char Ch;
	Ch cx[n+1];
	IdxTy i=0;
	for(;i<n; ++i) { 
		++pc;
	 	cx[i]=s[pc];
		if ( cx[i]==0) { return 1; } //
		} // i 
	cx[n]=0;
	int  v;
	std::istringstream(StrTy(cx)) >> std::hex>>v;
	//ss>>std::hex >> v;
	c[dp]=Ch(v); ++dp;
	++pc;
	return 0;
	} // numb

//////////////////////////////////////////////

// number identifierstart identifiermore  scope+ scope-  white eol quote escaep punc alpha digit   
// blocks {([< >])} and quotes '" and escape \   etc    
// classes white alpha numeric punc eol \n \r

static bool whites(const char c) { return (c==' ')||(c=='\t'); }
static bool mail(const char c) { return (c==' ')||(c=='\t')||(c==';'); }
template <class Tv> 
void parse_fullm(Tv & words, const char * s,const char sep, const IdxTy flags=0 )
{
const IdxTy sz=strlen(s);
bool preserve_quotes=Bit(flags,0);
// 2020-02-11 compiler complained .. 
//bool all_white_space=Bit(flags,1); // this is not really doing anything, func ptr 
bool mail_delims=Bit(flags,2); // white space and emicolons 
bool save_all=Bit(flags,2);
auto sepf=(&LineIterator::whites);
if (mail_delims)  sepf=(&LineIterator::mail);
bool deb=false;
//MM_ERR(MMPR3(sz,sep,s))
const IdxTy mxsz=10*sz+1;
char c[mxsz];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
bool esc=false;
bool quote=false;
while (s[pc] != 0 )
{
const char & cpc=s[pc];
// do not pass the quote to esc to dest 
if (!esc) if (cpc=='\\')  { esc=true; ++pc; continue; }  // else esc=false;
if (!esc) if (cpc=='"') {deb=!true; quote=!quote; ++pc ; if (preserve_quotes){  c[dp]=cpc; ++dp;} continue;  } 
const bool br=((*sepf)(cpc))&&(!esc)&&(!quote);
if (br){  c[dp]=0; const bool save_word=(start!=dp)||save_all;  if (save_word) { words.push_back(c+start);} ++dp; start=dp;    } else { c[dp]=cpc; ++dp; } 
++pc;
esc=false; // fck 2018-09-04 FUDD 
}
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 
if (esc||quote) { MM_ERR(" stray quote or escape "<<MMPR2(esc,quote)<<MMPR(s)) } 
if (start!=dp){ c[dp]=0;  words.push_back(c+start);  } 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

}//  parse_fullm





//  2014  pdftotext -layout xxx.pdf
// columns are generall separpated by multiple spaces, need some care 
template <class Tv> 
void parse_pdf_columns(Tv & words, const char * s,const char sep )
{
const IdxTy sz=strlen(s);
bool deb=false;
//MM_ERR(MMPR3(sz,sep,s))
char c[2*sz+1];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
IdxTy nsep=0;
bool esc=false;
bool quote=false;
while (s[pc] != 0 )
{
const char & cpc=s[pc];
// do not pass the quote to esc to dest 
if (!esc) if (cpc=='\\')  { esc=true; ++pc; nsep=0; continue; }  // else esc=false;
if (!esc) if (cpc=='"') {deb=!true; quote=!quote; ++pc ; nsep=0; continue;  } 
const bool br=(cpc==sep)&&(!esc)&&(!quote);
//if (br){ ++nsep; c[dp]=cpc; ++dp;if (c[start]==sep) ++start;    if ( nsep>2) { if (start!=(dp-nsep)) { c[dp-nsep]=0;  if (c[start]!=sep) { words.push_back(c+start);} start=dp; } }      } else { nsep=0; c[dp]=cpc; ++dp; } 
//if (br){ ++nsep; c[dp]=cpc; ++dp;if (c[start]==sep) ++start;    if ( nsep>2) { if (start<(dp-nsep)) { c[dp-nsep]=0;  if (c[start]!=sep) {StrTy w= StrTy(c+start);  Ss ss; ss<<start<<sep; w=ss.str()+w;  words.push_back(w);} start=dp; } }      } else { if ( c[start]==sep) start=dp; nsep=0; c[dp]=cpc; ++dp; } 
if (br){ ++nsep; c[dp]=cpc; ++dp;if (c[start]==sep) ++start;    if ( nsep>2) { if (start<(dp-nsep)) { c[dp-nsep]=0;  if (c[start]!=sep) {
push_pdf_word(words,c,sep,start);
} start=dp; } }      } else { if ( c[start]==sep) start=dp; nsep=0; c[dp]=cpc; ++dp; } 
++pc;
esc=false; // fck 2018-09-04 FUDD 
}
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

if (start!=dp){ c[dp-nsep]=0; //  words.push_back(c+start);  } 
push_pdf_word(words,c,sep,start); } 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

}//  parse_pdf_columns
template <class Tv >
void push_pdf_word(Tv & words, const char * c, const char sep, IdxTy start)
{
StrTy w= StrTy(c+start);  
Ss ss; ss<<start<<sep; w=ss.str()+w;  
words.push_back(w); 

}
////////////////////////////////////////////////////////


template <class Tv> 
void tokenize(Tv & words, const char * s, const IdxTy flags=0 )
{
const IdxTy sz=strlen(s);
bool preserve_quotes=Bit(flags,0);
const bool no_concat_seps=!Bit(flags,1);
const bool check_esc=Bit(flags,3);
const bool check_quote=Bit(flags,4);
const bool remove_crlf=!Bit(flags,5);
const bool skip_non_ascii=Bit(flags,6);
const bool debug=Bit(flags,8);
if (debug) 
{
MM_ERR(MMPR4(preserve_quotes,no_concat_seps,check_esc,debug)<<MMPR(skip_non_ascii)<<MMPR(flags)<<MMPR2(sz,s))
}
bool deb=false;
//MM_ERR(MMPR3(sz,sep,s))
char c[2*sz+1];
IdxTy pc=0;
IdxTy dp=0;
IdxTy  start=0;
bool esc=false;
bool quote=false;
//MM_ERR(" parser "<<MMPR4(no_concat_seps,preserve_quotes,flags,sep))
Iv cclass=0;
Iv cclass_old=0;
while (s[pc] != 0 )
{
// change order for skipping 
//cclass_old=cclass;
//const char & cpc=s[pc];
//
const char & cpc=s[pc];
if (skip_non_ascii) if ((cpc&128)!=0)  { ++pc; continue; } 
cclass_old=cclass;
//

cclass=m_char_table.lut(cpc);
/*
enum { BAD=~0, CARD=(1<<BITS),LC =(1<<0), DIGIT=(1<<1),
WHITE=(1<<2), CTRL=(1<<3), UC=(1<<4), PUNC=(1<<5),EOL=(1<<6),START=(1<<7),
STOP=(1<<8), OPER=(1<<9),UNUSED=(1<<10), UNUSEDB=10 } ;

// more terms than needed but overkill here is of no consequence. 
enum { ASCII=LC|DIGIT|WHITE|CTRL|UC|PUNC|EOL|START|STOP|OPER } ;
enum { ALPHA=LC|UC } ;
enum { ALPHANUM=LC|UC|DIGIT } ;
enum { PRINT=LC|DIGIT|WHITE|UC|PUNC|EOL|START|STOP|OPER } ;


*/
enum { PUNC=CharTable::PUNC,PRINT=CharTable::PRINT,
WHITE=CharTable::WHITE, 
ALPHANUM=CharTable::ALPHANUM, 
ALPHA=CharTable::ALPHA, 
DIGIT=CharTable::DIGIT, 
EOL=CharTable::EOL,
OPER=CharTable::OPER };
const bool is_eol=AnyMask(cclass,EOL);
const bool xistion=(cclass_old!=0)&&(
AnyMask(cclass,PUNC)
||AnyMask(cclass,OPER)
||!AnyMask(cclass,PRINT)
||(!AnyMask(cclass,WHITE)&&AnyMask(cclass_old,WHITE))
||(AnyMask(cclass,WHITE)&&!AnyMask(cclass_old,WHITE))
||(!AnyMask(cclass,DIGIT)&&AnyMask(cclass_old,DIGIT))
||(AnyMask(cclass,DIGIT)&&!AnyMask(cclass_old,DIGIT))
||(!AnyMask(cclass,ALPHA)&&AnyMask(cclass_old,ALPHA))
||(AnyMask(cclass,ALPHA)&&!AnyMask(cclass_old,ALPHA))
);

// do not pass the quote to esc to dest 
if (check_esc&&!esc) if (cpc=='\\')  { esc=true; ++pc; continue; }  // else esc=false;
if (check_quote&&!esc) if (cpc=='"') {
	deb=!true; 
	quote=!quote; 
	++pc ; 
	if (preserve_quotes)
		{  c[dp]=cpc; ++dp;} 
	continue;  
} 
//const bool br=(cpc==sep)&&(!esc)&&(!quote);
const bool br=(xistion)&&(!esc)&&(!quote);
if (br){  
	c[dp]=0; 
	// really pointless here now 
//	if ((c[start]!=0)||no_concat_seps ) 
		words.push_back(c+start); 
	++dp;  // dp=0???? 
	start=dp;
// put the sep into a new stinrg
// just remove the "else" below 
//	c[dp]=cpc; ++dp;    
} 
//else  - // ALWAYS keep the "sep" char(s)
{ if (!(is_eol&&remove_crlf)) {  c[dp]=cpc; ++dp;} } 
++pc;
esc=false; // fck 2018-09-04 FUDD 
} // s[pc]!=0
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

if (start!=dp){ c[dp]=0;  words.push_back(c+start);  } 
if (deb) { MM_ERR(MMPR4(pc,dp,start,words.size())) } 

}// tokenize 




/////////////////////////////////////////////////////////

IdxTy lineno() const { return m_serial; } 

// reparse the line with different rules, can not set must call
// after each line... 
// implement quoting or maybe line concat etc. 
void reparse(const IdxTy flags)
{



}
StrTy rest_of_line( const IdxTy skip) const 
{
Ss ss;
ss<<m_line;
StrTy n;
for(IdxTy i=0; i<skip; ++i) ss>>n;
return ss.str();
}

// the ok() call is not const
StrTy dump()  const
{
Ss ss;
ss<<" line "<<m_serial;
ss<<" words "<<size();
ss<<" ok "<<ok();
for (IdxTy i=0; i<size(); ++i) ss<<" word "<<i<<" "<<word(i); 

return ss.str(); 

}
private:
void release() const // yes, the ok() fking method fks all this sshot
{
//if (m_close) (*m_f).close();
if (m_close) delete m_f;
m_f=0;
m_close=false;
}
class stream_wrap
{
public:
stream_wrap(): m_f(0),m_close(false) {}
stream_wrap(IsTy * is, bool c ): m_f(is),m_close(c) {}
~stream_wrap() {} // we do not close here... 
IsTy * m_f;
bool m_close;
}; // stream_wrap
typedef std::vector<stream_wrap> StreamStack;
void Source(const StrTy &fn,const bool die)
{ 
MM_ERR(" Sourcing "<<fn<<" "<<die)
if (die)
	set(new FisTy(fn), true);
else	push(new FisTy(fn),true); 

}
	
public:
void set(IsTy * is, const bool close) { release(); m_f=is; m_close=close; } 

void push( char ** args, const IdxTy argc)
{
Ss* ssp= new Ss();
for(IdxTy i=0; i<argc; ++i) (*ssp)<<args[i]<<CRLF; 

push(ssp,true);

}


void push(IsTy * is, const bool close) { push_is(); m_f=is; m_close=close; } 
void push_is()
{
m_stream_stack.push_back(stream_wrap(m_f,m_close));
m_f=0;
}
bool pop_is() const 
{
//MM_ERR(" check pop ")
if (m_stream_stack.size()==0) return false;
//MM_ERR(" popping an old stream")
set(m_stream_stack.back());
m_stream_stack.pop_back(); // (stream_wrap(m_f,m_close));
return true; 
}
void clear_stack()
{ while (m_stream_stack.size()!=0) { pop_is(); } }
void set(const stream_wrap & sw) const 
{
release();
m_f=sw.m_f;
m_close=sw.m_close;
}
private:
mutable IsTy * m_f;
mutable bool m_close;
mutable StreamStack m_stream_stack;
Words m_words;
StrTy m_line;
IdxTy m_serial; 
IdxTy m_offset;
char m_sep;
IdxTy m_splitter;
StrTy m_prompt;
bool m_disable_readline;
bool m_save_last_line;
bool m_use_alt_eol;
bool m_remove_trailing_crlf;
bool m_remove_white;
char m_eol;
IdxTy m_eol_flags;

CharTable m_char_table;

}; // LineIterator




class ReadWriteMap
{
public:
typedef std::map<StrTy, StrTy> Map;
typedef std::vector<Map> MapStack;
ReadWriteMap( const StrTy & nm) {load(nm); }
ReadWriteMap( ) { }
// there was some good reason not to do this doh
//const Map & map() const { return m_map; } 
StrTy & operator[](const StrTy & k) { return m_map[k]; } 
Map::iterator begin() { return m_map.begin(); } 
Map::const_iterator begin()const  { return m_map.begin(); } 
Map::iterator end() { return m_map.end(); } 
Map::const_iterator end() const  { return m_map.end(); } 

Map::const_iterator find(const StrTy & k) const  { return m_map.find(k); } 
//const StrTy & operator[](const StrTy & k) const { return m_map[k]; } 
void clear() {Clear(); } 
void load( const StrTy & nm)
{
typedef std::ifstream IsTy; 
typedef std::stringstream Ss; 
typedef unsigned int IdxTy;
IsTy f;
m_nm=nm;
f.open(nm);
MM_ERR(" loading configuration file " <<nm)
IdxTy lines=0;
// this is not robust for crap in file lol 
while (f.good()&&!f.eof())
{
enum { LLEN=1<<16, RWLEN=LLEN-3};
char line[LLEN];
f.getline(line,RWLEN);
if ( !f.eof() &&!f.good())
{
MM_ERR(" line may be too long and missing rest of fie")
}
const IdxTy len=strlen(line); 
for(IdxTy i=0; i<len; ++i) if (line[i]=='=') { line[i]=' '; break; } 
Ss ss;
ss<<line;
StrTy n,v;
ss>>n; 
ss>>v; // while (ss.good()&&!ss.eof()) { StrTy x=""; ss>>x; v=v+StrTy(" ")+x; } 
const IdxTy vpos=n.length()+1;
if (vpos<len) v=StrTy(line+vpos);
m_map[n]=v;
++lines;
}
MM_ERR(" loaded "<<lines<<" lines from configuration file " <<nm)

f.close();
}

// this needs to work the samas the file laoder
void set_line(const StrTy & line)
{
typedef std::stringstream Ss; 
Ss ss;
ss<<line;
StrTy n,v;
ss>>n; ss>>v;
m_map[n]=v;

}

// this should return bool on setting 
void get(const StrTy & nm, StrTy & v)
{ if (m_map.find(nm)!=m_map.end()) v=StrTy(m_map[nm]); }

void get(const StrTy & nm, double & v)
{ if (m_map.find(nm)!=m_map.end()) v=::atof(m_map[nm].c_str()); }

void get(const StrTy & nm, int & v)
//{ if (m_map.find(nm)!=m_map.end()) v=::atoi(m_map[nm].c_str()); }
{ if (m_map.find(nm)!=m_map.end()) v=myatoi(m_map[nm].c_str()); }

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); }
static int myatoi(const char * c) {const bool inv=(c[0]=='~') ; int x= ::strtol(c,0,0); return inv?(~x):x;  }

// really should check... 
bool  get(const StrTy & nm, unsigned int & v)
//{ if (m_map.find(nm)==m_map.end()) return false; v=(unsigned int)(::atoi(m_map[nm].c_str())); return true; }
{ if (m_map.find(nm)==m_map.end()) return false; v=(unsigned int)(myatoi(m_map[nm].c_str())); return true; }


template <class Ty> 
const StrTy  get_string(const StrTy & nm, const Ty & def ) const
{ 
	if (m_map.find(nm)==m_map.end()) return StrTy(def);
	const  StrTy s=(*(m_map.find(nm))).second; 

return s;  }
const double    get_double(const StrTy & nm, const double & def ) const
{ 
	if (m_map.find(nm)==m_map.end()) return def;
	const  StrTy s=(*(m_map.find(nm))).second; 

return ::atof(s.c_str());  }

int  get_int(const StrTy & nm, const int & def) const 
{ 
	if (m_map.find(nm)==m_map.end()) return def;
	const  StrTy s=(*(m_map.find(nm))).second; 

	return ::strtol(s.c_str(),0,0); 
//return ::atoi(s.c_str()); 

}

unsigned int  get_uint(const StrTy & nm, const unsigned int & def) const 
{ 
	if (m_map.find(nm)==m_map.end()) return def;
	const  StrTy s=(*(m_map.find(nm))).second; 

	return (unsigned int ) ::strtol(s.c_str(),0,0); 
//	return ::atoi(s.c_str());  


}



// the string in the param list generates a char[] NOT a char * WT FUDD????
///template <class Tfck> bool get_bool(const Tfck[]   nm, const bool def) { return get_bool(StrTy(nm),def); }
// "fck" should put char * in the fcking list fck 
template <class Tfck> bool get_bool(const Tfck   nm[], const bool def) { return get_bool(StrTy(nm),def); }
bool get_bool(const StrTy & nm, const bool def) const 
{ 
	if (m_map.find(nm)==m_map.end()) return def;
	const  StrTy s=(*(m_map.find(nm))).second; 
	if (s=="true") return true;
	if (s=="false") return false;
	if (s=="1") return true;
	if (s=="0") return false;
	if (s=="yes") return true;
	if (s=="no") return false;
	MM_ERR( "returning false for bool "<<s<<" for key "<<nm<<" in "<<m_nm)
	return false;

}
static bool str_to_bool(const StrTy & s) 
{
	const char * p=s.c_str();
	const char c=p[0];
	if (c=='1') return true;
	if (c=='t') return true;
	if (c=='T') return true;
	if (c=='y') return true;
	if (c=='Y') return true;
	if (c=='+') return true;

	if (c==0) return false;
	if (c=='N') return !true;
	if (c=='n') return !true;
	if (c=='F') return !true;
	if (c=='f') return !true;
	if (c=='0') return !true;
	if (c=='-') return !true;
//	if (s=="true") return true;
//	if (s=="false") return false;
//	if (s=="1") return true;
//	if (s=="0") return false;
//	if (s=="yes") return true;
//	if (s=="no") return false;
	MM_ERR( "returning false for bool "<<s)
	return false;
}

// for things like rational class
template <class Ty> bool  get_read(const StrTy & nm, Ty & v)
{ if (m_map.find(nm)==m_map.end()) return false; v.read(m_map[nm]); return true; }



// doubles may not be exact, verify ... 
template <class Ty> bool set(const StrTy & nm, const Ty & v, const bool always=true)
{
auto ii=m_map.find(nm);
if ( !always) if ( ii!=m_map.end()) return false; 
std::stringstream ss;
ss.precision(64);
ss<<v;
m_map[nm]=ss.str();
return true;
}
void push() { m_saved.push_back(m_map); }
void pop() { if ( m_saved.size()!=0) m_map=m_saved.back(); m_saved.pop_back(); }
enum { SAME_MAP=1, NO_MAP=-1, DIFFER=0 };
int same() 
{ 
	if (m_saved.size()==0) return NO_MAP; 
	if (m_map==m_saved.back()) return SAME_MAP ; else return DIFFER; 
}
template <class Tx> void get_all(Tx & keys, Tx & values) const 
{
for(auto ii=m_map.begin(); ii!=m_map.end(); ++ii)
{
keys.push_back((*ii).first);
values.push_back((*ii).second);
}

}
Map & map() { return m_map; }
const Map & map()const  { return m_map; }

void Clear() { m_map.clear(); m_saved.clear(); } 
// Map & map()  { return m_map; }

StrTy m_nm;
Map m_map;

MapStack m_saved;

} ; // ReadWriteMap

#define MM_INC_MSG(cm,msg) cm.inc_location( __FILE__,__LINE__,msg);
#define MM_INCP_MSG(pcm,msg) { if (pcm!=0) { (*pcm).inc_location( __FILE__,__LINE__,msg);}}


class CounterMap
{


public:
typedef unsigned int IdxTy;
typedef double R;
typedef long long TimeTy;
typedef std::stringstream Ss;
typedef std::string StrTy;
typedef std::map<StrTy,IdxTy> MapTy;
typedef std::map<int,R> ArrayTy;
typedef std::map<StrTy,R> RMapTy;
typedef std::map<StrTy,TimeTy> TMapTy;
typedef std::map<StrTy,ArrayTy> ArrayMapTy;

CounterMap(): m_zed(markabs()) {mark("mytzed"); }
CounterMap(const CounterMap & that ): m_zed(that.m_zed) {}
CounterMap & inc(const StrTy & nm) 
{ if ( m_map.find(nm)!=m_map.end()) ++m_map[nm] ; else m_map[nm]=1;  return *this; }

CounterMap & inc(const StrTy & nm, IdxTy del) 
{ if ( m_map.find(nm)!=m_map.end()) m_map[nm]+=del ; else m_map[nm]=del;  return *this; }

const IdxTy   get_count(const StrTy & nm) const 
{ auto ii=m_map.find(nm);
if (ii==m_map.end()) return 0;   return (*ii).second; }

const R   get_time(const StrTy & nm) const  // in SECONDS 
{ auto ii=m_tmap.find(nm);
if (ii==m_tmap.end()) return 0;   return R((*ii).second)/R(1000); }

// for counting error message hits
template <class Tx, class Ty, class Tz > CounterMap & inc_location(const Tx & file, const Ty & line, const Tz & msg) 
{ 
Ss ss;
ss<<file; ss<<line; ss<<":"<<msg;
const StrTy nm=ss.str();
// others may look up the message without needing line info  hence inc msg. 
if ( m_map.find(nm)!=m_map.end()) ++m_map[nm] ; else m_map[nm]=1;  inc(msg,1); return *this; }



// add to array nm at index idx
CounterMap &add(const StrTy & nm, const int i, const  R & val) 
{ m_arrays[nm][i]+=val; return *this; }

CounterMap &add(const StrTy & nm, const  R & val) 
{ if ( m_rmap.find(nm)!=m_rmap.end()) m_rmap[nm]+=val ; else m_rmap[nm]=val;
inc(nm);
 return *this; }
CounterMap &  set(const StrTy & nm, const  R & val) 
{ m_rmap[nm]=val; inc(nm);return * this; } 
CounterMap &  mark(const StrTy & nm) 
{ m_tmap[nm]=mark(); return * this; } 

// need to distinguish intervals from absolute times.. 
CounterMap &  lap(const StrTy & total, const StrTy & m) 
{ m_tmap[total]=mark()-m_tmap[m]; return * this; } 

CounterMap &  cum(const StrTy & total, const StrTy & m) 
{ m_tmap[total]+=mark()-m_tmap[m]; return * this; } 

CounterMap &  cumdel(const StrTy & total, const StrTy & m) 
{ m_tmap[total]+=mark()-m_tmap[m]; m_tmap.erase(m); return * this; } 


CounterMap &  minmax(const StrTy & nm, const R & val, const bool min) 
{ 
if ( m_rmap.find(nm)==m_rmap.end()){ m_rmap[nm]=val; return *this; } 
if ( min) { if (val<m_rmap[nm]) m_rmap[nm]=val; }
else { if (val>m_rmap[nm]) m_rmap[nm]=val; }

return * this;

} 


// want some feature to stack times or notes or something. 
//CounterMap &  push(const StrTy & total, const StrTy & m) 
//{ m_tmap[total]+=mark()-m_tmap[m]; return * this; } 

CounterMap &  reset(const StrTy & nm) 
{ cif(m_map,nm); cif(m_tmap,nm); cif(m_rmap,nm); return * this; } 

CounterMap &  clear(const StrTy & nm) 
{ cif(m_map,nm); cif(m_tmap,nm); cif(m_rmap,nm); return * this; } 

template <class Ty> void cif(Ty &m, const StrTy & nm)
//{ if ( m.find(nm)!=m.end()) m[nm]=0;  }
{ if ( m.find(nm)!=m.end()) m.erase(nm);   }

//{ m_arrays[nm][i]+=val; }
StrTy array_to_string(const StrTy & nm)
{
	Ss x;
	ArrayTy & y=m_arrays[nm];
	IdxTy cnt=0; 
	for ( auto ii=y.begin(); ii!=y.end(); ++ii)
	{
		if (cnt!=0) x<<",";  ++cnt;
		x<<"["<<(*ii).first<<"]="<<" "<<(*ii).second; // <<",";
	}
	return x.str();
}
template<class Os> Os & dump(const StrTy & lbl, Os & os, const StrTy & sep="=", const StrTy & sep2=";" ) 
{
os<<lbl<<" counts:  "; 
auto ii=m_map.begin();
auto ie=m_map.end();
while (ii!=ie) { os<<(*ii).first<<sep<<(*ii).second<<sep2; ++ii; }
os<<std::endl;
os<<lbl<<" vars:  "; 
auto iir=m_rmap.begin();
auto ier=m_rmap.end();
while (iir!=ier)
{ os<<(*iir).first<<sep<<(*iir).second<<sep<<((*iir).second/m_map[(*iir).first])<<sep2; ++iir; }
os<<std::endl;
os<<lbl<<" times:  "; 
auto iit=m_tmap.begin();
auto iet=m_tmap.end();
while (iit!=iet) { os<<(*iit).first<<sep<<format_time((*iit).second)<<sep2; ++iit; }
os<<std::endl;
return os; 

}

StrTy time_table(const StrTy label)
{
lap("mytzed","mtzed");
std::stringstream s;
typedef mjm_format_text FoTy;
const StrTy sep=" ";
auto iit=m_tmap.begin();
auto iet=m_tmap.end();
// IdxTy max_sz_map(Ti  ii, const Te & ff)
IdxTy nmsz=FoTy::max_sz_map(iit,iet);
while (iit!=iet) 
{ 
	s<<label<<" "<<FoTy::left((*iit).first,nmsz)<<sep<<FoTy::align_on(format_time((*iit).second),StrTy("."),5,14,0)<<CRLF; 

++iit; }
//os<<std::endl;

return s.str();
} 


static StrTy format_time(const TimeTy & t) 
{
static const TimeTy hour=1000*60*60;
static const TimeTy minute=1000*60;
static const TimeTy second=1000;
std::stringstream s;
TimeTy tx=t;
bool start=false; 
if ((tx>=hour)||start) { TimeTy ty=tx/hour; s<<ty<<":"; tx-=ty*hour;start=true; } 
if ((tx>=minute)||start) { TimeTy ty=tx/minute; s<<ty<<":"; tx-=ty*minute;start=true; } 
if ((tx>=second)||start) { TimeTy ty=tx/second; s<<ty<<""; tx-=ty*second;start=true; } 
s<<"."<<tx; 
return s.str(); 
}


TimeTy markabs() const
{  struct timeval tp;
    gettimeofday(&tp, NULL);
    TimeTy  mslong = TimeTy(tp.tv_sec) * 1000L +TimeTy( tp.tv_usec / 1000); 
	return mslong;
}

TimeTy mark() const { return markabs()-m_zed; } 

StrTy format_entry(const StrTy & nm)  const
{ 
std::stringstream s;
s<<nm<<"=";
if ( m_map.find(nm)!=m_map.end()) 
	s<<(*(m_map.find(nm))).second; 
return s.str(); 
}
void clear()
{
m_map.clear(); m_rmap.clear(); m_arrays.clear(); 
}
MapTy m_map;
RMapTy m_rmap;
TMapTy m_tmap;
ArrayMapTy m_arrays;
// this creates a problem for copy ctor
// const 
TimeTy m_zed;
// this finally provoked a code sensitivity not revealed by optimization, dphi indicies bad
//double pad[100];
};

template <class Tag, class CodeKey> class mjm_code_info
{

public:
typedef Tag info_tag;
typedef CodeKey code_key;
class Tr{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::D D;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::StrTy StrTy;

typedef typename std::map<code_key,StrTy> CodeMap;

// the look up needs to use find to be const fix later or leak 
//static const StrTy &  code(const code_key & key) { return codes()[key]; }
static  StrTy &  code(const code_key & key) { return codes()[key]; }


 static CodeMap & codes() 
{
static CodeMap x=CodeMapInit();
return x;
}
// other classes put stuff here 
 static CodeMap  CodeMapInit() 
{
CodeMap x;
return x;
}

} ; // mjm_code_info


/*

Subrotines, variables etc handled by bash.
User interface with up arrow etc could be nice.
Needs a way to shift command words for hierarchial commands
Handling quotes etc would be nice too. 

*/

class CommandInterpretter : public LineIterator
{
typedef LineIterator Super;
typedef CommandInterpretter Myt;
typedef std::stringstream Ss; 
typedef unsigned int IdxTy;
typedef StrTy Word;
typedef StrTy Line;
public:
typedef std::vector<Word> Words;
typedef std::vector<Line>  Lines;
typedef Lines Subroutine;
typedef StrTy Symbol;
typedef StrTy Label;
typedef std::map<Symbol,Subroutine> Subroutines;
typedef std::map<Label,IdxTy> Labels;


enum { CODE_SPLIT_BASE=Super::CODE_SPLIT_BASE,MAX_CODE_SPLIT=Super::MAX_CODE_SPLIT};
enum { KVP_SPLIT_BASE=Super::KVP_SPLIT_BASE};
CommandInterpretter() : Super() {}
CommandInterpretter(IsTy * f) : Super(f) {}
// this is a file or maybe a pipe name 
CommandInterpretter(const StrTy & nm) : Super(nm) {}
enum {BAD=~0};
IdxTy bad() const { return BAD; } 
virtual void next() { Super::next(); }

static StrTy cmd_set(const LineIterator & li) { return li.word(1)+StrTy(" ")+li.word(2); }
StrTy cmd_set() const { return cmd_set(*this); } 
static bool cmd_ok(const LineIterator & li, const IdxTy sz, const IdxTy n)
//static bool cmd_ok( LineIterator & li, const IdxTy sz, const IdxTy n)
{
 if (sz<n) { MM_ERR("ERROR "<<li.dump())  return false; }
return true;

}
bool cmd_ok( const IdxTy n) { return cmd_ok(*this,size(),n) ; }

static IdxTy n(const LineIterator & li, const IdxTy sz, const IdxTy p)
{
IdxTy n=~0;
 if (cmd_ok(li,sz,p+1))  n=::atoi((li.word(p).c_str()));
return n;
}
IdxTy n( const IdxTy p) { return n(*this,size(),p) ; }
static double d(const LineIterator & li, const IdxTy sz, const IdxTy p)
{
double d =0;
 if (cmd_ok(li,sz,p+1))  d=::atof((li.word(p).c_str()));
return d;
}
double d( const IdxTy p) { return d(*this,size(),p) ; }

virtual const bool  help()
{
	
 if (!cmd_ok(2))  return false; 
 const StrTy & cmd=word(1);
	if (cmd=="-help") return true;
	if (cmd=="help") return true;
	if (cmd=="--help") return true;
return false; 
}
IdxTy find( const StrTy & cmd) const
{
for (IdxTy i=0; i<size(); ++i) if (word(i)==cmd) return i;
return BAD; 
}
void set(const StrTy & cmds, const IdxTy flags)
{
const bool parse = Bit(flags,0);
const bool semi_split = !Bit(flags,1);

if (flags==0) { Super::set(new Ss(cmds), true);  return; } 
//if (flags==1)
if (parse) 
{
Ss  * ss =new Ss();
for (IdxTy i=0; i<cmds.length(); ++ i)
{
const char c=cmds.c_str()[i];
if ((c==';')&&semi_split) (*ss)<<CRLF; else (*ss)<<c;
}
// apparently it takes the last line ok 
//(*ss)<<CRLF;
//MM_MSG(" fcking string is fck "<<(*ss).str())
//MM_ERR(" fcking string is fck "<<(*ss).str())
Super::set(ss, true); 
}
}
void source(const StrTy & fn, const bool die=false) { Super::source(fn,die); } 





private:
// needs an interactive interface for up arrow etc. Maybe a curses terminal lol 
Lines m_code;
// right now do all of this in bash
Labels m_labels;

Subroutines m_subs;

}; // CommandInterpretter
class CommandInterpretterVector
{
public:
typedef  CommandInterpretter Ci;
typedef std::vector<Ci> Civ;
typedef unsigned int IdxTy; // TODO FIXME only needed by mjm_snacks??
CommandInterpretterVector( const IdxTy n) : m_civ(n) {}
CommandInterpretterVector( )  {}
// give the line or stream to ci n and get the words from ci 0


Civ m_civ;

};  //  CommandInterpretterVector


class CommandInterpretterParam 
{
typedef unsigned int IdxTy; 
typedef std::string StrTy;
typedef std::ostream OsTy;
public:
typedef std::vector<StrTy> Words; 
//  TODO include order messed up somewhere, header lib mayu be bad  
CommandInterpretterParam( CommandInterpretter & li, const IdxTy min=2)
:m_li(li),help_only(false),m_os(&std::cout) // , m_variables(0)
{
MM_LOOP(ii,li.words()) { w.push_back((*ii)); }
while (w.size()<min) { w.push_back(StrTy()); } 
if (w.size()>1) p1=w[1];
if (w.size()>2) p2=w[2];
{ MM_SZ_LOOP(i,li.words(),liws) { if (i>1) w2.push_back(li.words()[i]); }}
}

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

const StrTy & cmd() { return w[0]; } 
const StrTy & line() { return m_li.line(); } 
const Words& words() const { return w; }
const Words& words2() const { return w2; }
const IdxTy size() const  {return w.size(); } 
const IdxTy params()const {if (w.size()==0) return 0;  return w.size()-1; } 
const StrTy  wif(const IdxTy i)  const 
{ if (i<w.size()) return w[i]; return StrTy(""); } 

const StrTy  eval(const IdxTy i) 
{ if (i<w.size()) 
{
StrTy x=w[i];
if (vars())
{ 
const char * c=x.c_str();
if (c[0]==0) return x; 
//if (c[0]=='$') {return (*m_varaibles).get(StrTy(c+1)); } 
}
return x;
}
return StrTy(""); } 
// if this is only a request for help
bool help() const { return help_only; }
// method can call this with help info for uniform handling 
bool help(const StrTy & hinfo) { MM_ERR(hinfo); return true; }
void set_help(const StrTy & hinfo) { m_help_info=hinfo; }


/*
StrTy dump_variables() const
{
Ss ss;
std::vector<StrTy> k,v;
if (! vars()) { return ss.str(); }  
//(*m_variables).get_all(k,v);
MM_SZ_LOOP(i,k,sz)
{
ss<<k<<"="<<v<<CRLF;
}

return ss.str();
}
*/

CommandInterpretter & li() { return m_li; }
bool vars() const { return false; } // m_variables !=0; }
CommandInterpretter & m_li;
bool help_only;
std::vector<StrTy> w; 
std::vector<StrTy> w2; 
StrTy p1,p2;
OsTy * m_os;
StrTy m_help_info;
//mjm_logic_base * m_variables;
}; // CommandInterpreterParam 






#endif

