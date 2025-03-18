#ifndef MJM_TAXONOMY_STUFF_H__
#define MJM_TAXONOMY_STUFF_H__

// #warning "made from mjm_collections 2021-01-25 not tested "
#include "mjm_collections.h"
#include "mjm_gen_viz.h"
#if 0
#include "mjm_globals.h"
#include "mjm_string_tokenizer.h"
// needed for treelayout apparently
#include "mjm_gen_viz.h"
#include "mjm_data_model_error_log.h"
#include "mjm_generic_iterators.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#endif

#define LUT_LOAD_FWD  0

// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#ifndef MM_DMEL
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
#endif

#define MM_DMELF(fm,code,x) f.DMel(code, f.m_ss<<MM_STR_LOC<<" "<<x); 

#if 0 
namespace collections_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel;
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // biom_hdf5_traits

#if 0

// define input and consistent inverse output formats for supporting
// data structures. This is csv with quote and escape supporting
// ivector elements with decorated integers or footnotes.  
template <class Tr=collections_traits::Tr > class form_def
{
typedef form_def<Tr>  Myt;
//typedef collections_traits::Tr Tr;
protected:
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;


public:
form_def(): m_sep(",") {}
// saving operations 
StrTy operator()(const StrTy & s)const  { return StrTy("\"")+s+StrTy("\""); } 
StrTy operator()(const IdxTy & i)const  {Ss s; s<<i; return s.str(); } 


StrTy operator()(const std::vector<IdxTy> & v)const  
{Ss s;s<<"\""; MM_LOOP(ii,v) {if (ii!=v.begin()) s<<sep();   s<<(*ii); } s<<"\""; return s.str(); } 


StrTy operator()(const std::vector<IdxTy> & v, const std::vector<StrTy> & footnotes)const  
{Ss s;s<<"\""; 
const IdxTy szf=footnotes.size();
MM_SZ_LOOP(i,v,sz)
//MM_LOOP(ii,v) 
{if (i!=0) s<<sep();   s<<v[i]; if (i<szf) { s<<footnotes[i];}  } 
s<<"\""; return s.str(); } 



StrTy operator()(const int & i)const  {Ss s; s<<i; return s.str(); } 
 // loading ops
Myt &  operator()(StrTy & d, const StrTy & s)  { d=s; return *this; } 
Myt &  operator()(unsigned int  & d, const StrTy & s)  { d=myatoi(s); return *this; } 


Myt &  operator()(std::vector< unsigned int>   & d, const StrTy & s)  
{
StrTy r=imbue_taboo(s);
  Ss ss(r); 
while ( true) {StrTy x=""; ss>>x; if (x=="") break;   d.push_back( myatoi(x)) ; } return *this; } 

Myt &  operator()(std::vector< unsigned int>   & d, std::vector<StrTy> &footnotes, const StrTy & s)  
{
StrTy r=imbue_taboo(s);
  Ss ss(r); 
while ( true) {StrTy x=""; ss>>x; if (x=="") break;   
const IdxTy sz=x.length();
const char * p=x.c_str();
//char y[sz+1];
StrTy note="";
for (IdxTy i=0; i<sz; ++i) 
{
	const char c=p[i];  
	if (!::isdigit(c)&&(c!='-') ) 
	{ 
		note=StrTy(p+i); break; 
	}  
} 
d.push_back( myatoi(x)); 
footnotes.push_back(note);
}
 return *this; } 



Myt &  operator()( int  & d, const StrTy & s)  { d=myatoi(s); return *this; } 
Myt &  operator()( double  & d, const StrTy & s)  { Ss ss(s); ss>>d; return *this; } 

const StrTy&  sep() const { return m_sep; } 

StrTy imbue_taboo(const StrTy & s) const
{const IdxTy sz=s.length(); char r[sz+1]; r[sz]=0;  const char sep=m_sep.c_str()[0]; 
for(IdxTy i=0; i<sz; ++i) { char c=s[i]; if (c==sep) r[i]=' '; else r[i]=c; }
return StrTy(r);
}


int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }


void setup(CommandInterpretter & li) {
li.set_split(1,m_sep.c_str()[0]);

}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}
StrTy m_sep;
Ss m_ss;

}; // form_def


#endif
















#if 0 
template <class Tv,class Tr > class vector_iterator
{ // make an iterator that makes a vecror look like a map with pair 
typedef vector_iterator<Tv,Tr> Myt; 
//typedef  collections_traits::Tr Tr;
typedef Tr::IdxTy IdxTy;
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
const Tpair  operator*()const  { return Tpair(m_i,m_v[m_i]) ;} 
Myt  end() const { return Myt(m_v, m_end.first); } 
Myt & operator++() { ++m_i;  return *this; } 
bool operator==(const Myt & that) { return (&m_v[m_i]==&that.m_v[that.m_i]); } 
bool operator!=(const Myt & that) { return !((*this)==that); } 
Myt  operator=(const Myt & that ) { Myt x(m_v,m_i); return x; } 
private:
IdxTy m_i;
const Tvec & m_v;
Tpair m_end;
}; // vector_iterator

#endif


#if 0

class null_string_tokenizer
{
public:
template <class Ty> const Ty&  operator()( const Ty &  n ) const { return n; }
}; // null_string_tokenizer


// convert strings to numbers for high frequency strings with minimal
// effort. Needs to take care of accidental adds etc.
class string_tokenizer
{
typedef string_tokenizer Myt;
typedef collections_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy, IdxTy> Fmap;
typedef std::vector<StrTy> Rmap;
public:
typedef Myt mapped_type;
string_tokenizer():BADINT("invalid") {}
const StrTy&  operator()( const IdxTy &  n ) const {if ( n<m_vec.size())  return m_vec[n]; return BADINT;  }
const IdxTy&  operator()( const char *  k  )  { return (*this)(StrTy(k)); } 
const IdxTy&  operator()( const StrTy & k  )  {
auto x= m_map.find(k); if (x!=m_map.end()) return (*x).second;
const IdxTy n=m_vec.size();
m_map[k]=n;
m_vec.push_back(k);
 //return n; 
 return m_map[k];  // TODO just to return a ref lol?

}
const IdxTy&  operator[]( const StrTy & k  )const   {
auto x= m_map.find(k); if (x!=m_map.end()) return (*x).second;

 return badidx();

}
Fmap::iterator find(const StrTy & k )  { return m_map.find(k); } 
Fmap::const_iterator find(const StrTy & k ) const { return m_map.find(k); } 


typedef vector_iterator<StrTy,collections_traits::Tr > Vi;
typedef const_vector_iterator<StrTy,collections_traits::Tr > CVi;

//Fmap::iterator 
Vi find(const IdxTy & k )  
{ //MM_TRACE 
MM_ONCE(" this should not be called ",) 
return Vi(m_vec,k); // m_map.end();

 } 
//Fmap::const_iterator 
CVi find(const IdxTy & k ) const 
{//MM_TRACE  
MM_ONCE(" this should not be called ",) 

return CVi(m_vec,k); // m_map.end(); 

} 


const IdxTy & badidx() const { static IdxTy b=~0;  return b; }

template <class Td, class Ts>
void tokenize_vector(Td & d, Ts & s)
{
MM_LOOP(ii,s)
{
d.push_back((*this)(*ii));
}

}
// this is only for dumping which now should be a member func 

#if 0 
Fmap::const_iterator begin() const { return m_map.begin(); } 
Fmap::const_iterator end() const { return m_map.end(); } 
#endif

//typedef vector_iterator<std::vector<StrTy> > Vi;


CVi begin() const  { return CVi(m_vec); } 
CVi end() const  { return CVi(m_vec).end(); } 

Vi begin()   { return Vi(m_vec); } 
Vi end()   { return Vi(m_vec).end(); } 



// this iterates over the REVERSE map 
// ???? 

static StrTy  header(IdxTy & i ) 
{
switch (i)
{
case 0: { return StrTy("id"); } 
case 1: { return StrTy("value"); } 
case 2: { return StrTy(""); } 
default:
MM_ERR(" bad case need DMEL "<<MMPR(i))
}
//return h;
return StrTy("");
} 
template <class Ti, class Tf> 
const StrTy  field(const Ti & ii, const IdxTy & i, Tf & f )  const
{
switch (i)
{
//case 0: { return m_id; } 
case 0: { return f((*ii).first); } 
case 1: { return f((*ii).second); } 
//case 1: { return m_parent; } 
//case 2: { return m_rank; } 
//case 3: { return m_src; } 
default:
			MM_DMELF(f,"badheder ",i);  
//MM_ERR(" bad case need DMEL "<<MMPR(i))
}
return f(badidx());
} 
template <class Ti, class Tf> 
void  field(const Ti & _key, const IdxTy & i, Tf & f,const StrTy & w ) 
{
const IdxTy key=f.myatoi(_key);
switch (i)
{
case 0: { break; } 
case 1: { 
const IdxTy n=m_vec.size();
if (n!=key)
{
MM_DMELF(f,"loadingtokenizer",MMPR4(key,n,w,i))

}
m_map[w]=n;
m_vec.push_back(w);
 //return n; 
// (*find(key)).second= f.myatoi(w); 
//return f((*ii).first);
break; 
 } 
//case 1: { return f((*ii).second); } 
//case 1: { return m_parent; } 
//case 2: { return m_rank; } 
//case 3: { return m_src; } 
default:
			MM_DMELF(f,"badheder ",i);  
//MM_ERR(" bad case need DMEL "<<MMPR(i))
}
//return f(badidx());
} 
void clear() { m_map.clear(); m_vec.clear(); } 

private:
Fmap m_map;
Rmap m_vec;
StrTy BADINT;
}; // string_tokenizer


#endif





template <class Tobj, class Traits >
class mjm_collection 
{
typedef mjm_collection Myt;
protected:
typedef Traits Tr;
typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

typedef string_tokenizer St;

typedef typename To::const_iterator CItor;
typedef typename To::iterator Itor;

public:
mjm_collection(): m_dmel(0) {}
virtual ~mjm_collection() {}

// never used, see the fields( params 
void load_object(const StrTy & str) {} // { std::ifstream is (str); load_canon_map(is); }
void load_object(IsTy & is)
{
 //   const IdxTy base=0;
    IdxTy cnt=0;
    IdxTy skip=0;
    CommandInterpretter li(&is);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		load_line(li);
        //MM_ERR(" processing "<<li.dump())
//        if (sz<(base+3)) { ++skip; continue; }
//        const IdxTy level=myatoi(li.word(base));
//        const StrTy & t=li.word(base+1);
//        const StrTy & c=li.word(base+2);
//        m_fwd[level][t]=c;
        ++cnt;
    } // li 
    MM_ERR(" finish loading canon "<<MMPR2(cnt,skip))
} // load 
virtual void load_line(const CommandInterpretter & li ) {}
void save_object(OsTy & os)const  {}
void dump_object(const StrTy & label, const IdxTy & flags, OsTy & os) const 
{
MM_ERR(" dumping canon map ")
const StrTy sep=" ";
MM_LOOP(ii, m_obj)
{
dump_item(ii,os);
} // ii 
} // dump_object
virtual void dump_item( const CItor & ii, const IdxTy & flags, OsTy & os) const
{

} 


protected:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }

void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}


Dmel * m_dmel;
Ss m_ss;
To m_obj;

};  // mjm_collection

// turning into an R data frame with ad hoc crap lol 
class mjm_ragged_table
{
typedef collections_traits::Tr Tr;
typedef mjm_ragged_table Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
public:
typedef StrTy Word;
typedef std::vector<Word> Line;
typedef std::vector<Line> Lines;
typedef Lines::const_iterator LiCItor;
typedef Lines::iterator LiItor;
typedef std::map<StrTy, StrTy>  xlate_map;
typedef std::map<StrTy, D>  xlate_double_map;
typedef std::vector<IdxTy> LocVec;
typedef std::map<StrTy, LocVec> line_index_map;

typedef std::vector<IdxTy> sort_order;
mjm_ragged_table()
: m_min_words(1),m_max_words(0),m_sep(" "),m_splitter(1),m_ignore_hash(false),m_lip(0) {}

mjm_ragged_table(const IdxTy lines, const IdxTy words)
: m_min_words(1),m_max_words(0),m_sep(" "),m_splitter(1),m_ignore_hash(false),m_lip(0) {

 m_lines=Lines(lines);
const Line v= Line(words);
MM_LOOP(ii,m_lines) (*ii)=v;
}

~mjm_ragged_table() { delete m_lip; m_lip=0;  } 
LiCItor begin() const { return m_lines.begin(); } 
LiCItor end() const { return m_lines.end(); } 
LiItor begin()  { return m_lines.begin(); } 
LiItor end()  { return m_lines.end(); } 

const Line & operator[] (const IdxTy i )  const { return m_lines[i]; } 
Line & operator[] (const IdxTy i )   { return m_lines[i]; } 
// TODO FIXME need to make this more sane. do not call

IdxTy maxcol() const 
{
IdxTy mc=0;
MM_LOOP(ii,(*this)) { if ((*ii).size()>mc) mc=(*ii).size(); } 
return mc; 
}

void load(const StrTy & fn ,const bool debug=false)
{
//std::map<StrTy,IdxTy> m;
//std::istream & is =  std::cin;
std::ifstream  is(fn); //  =  std::cin;
load(is,debug);
}
void load_lines(const StrTy & fn ,const bool debug=false)
{
//std::map<StrTy,IdxTy> m;
//std::istream & is =  std::cin;
std::ifstream  is(fn); //  =  std::cin;
load_lines(is,debug);
}

sort_order ordering() const
{
sort_order x(size());
IdxTy i=0;
MM_LOOP(ii,x) {(*ii)=i; ++i; } 
return x;
}
// const Ty and Myt? 
//template <class Ty> void sort( sort_order & o, Ty & f) 
//{ std::sort(o.begin(),o.end(),f); }

template <class Ty> void sort( sort_order & o, Ty  f) 
{ std::sort(o.begin(),o.end(),f); }

template <class Ty> void sort(const IdxTy off,  sort_order & o, Ty  f) 
{ 
auto ii=o.begin();
for(IdxTy i=0; i<off; ++i) { ++ii; if (ii==o.end()) return; } 
std::sort(ii,o.end(),f); }


void pdf_columns()
{
m_sep=" ";
m_splitter=2;
}
void sep(const StrTy & s) { m_sep=s;} 
void splitter(const IdxTy & s) { m_splitter=s;} 
const Line & line(const IdxTy i ) const { return m_lines[i];}
Line & ncline(const IdxTy i )  { return m_lines[i];}
void clear()  {  m_lines.clear();}
const IdxTy  size() const { return m_lines.size(); } 
const void ignore_hash(const bool x )  { m_ignore_hash=x; } 
void push_back(const Line & li ) { m_lines.push_back(li); }  
void add(const Line & li ) { m_lines.push_back(li); }  
Myt& operator<<(const Line &li) {add(li); return (*this); }
///////////////////////////////////////////////////
//void Add( const Line & x)
//{
//	const IdxTy sz=x.size();
//	if (sz<m_min_words) m_min_words=sz;
//	if (sz>m_max_words) m_max_words=sz;
//	m_lines.push_back(x); 
//}
Myt & new_line() 
{
// the incremental approach wrecks min and max word counts per line,,,, 
m_lines.push_back(Line());
return *this; 
}
Myt & fresh_line() 
{
// the incremental approach wrecks min and max word counts per line,,,, 
if (m_lines.size()==0) return *this;
if (m_lines.back().size()==0) return *this; 
m_lines.push_back(Line());
return *this; 
}

StrTy sum(const IdxTy i, const IdxTy first, const IdxTy last,const StrTy sep="") const
{
Ss ss;
const Line & l=(*this)[i];
const IdxTy sz=l.size();
const IdxTy f=(sz>last)?last:sz;
if (f>first) ss<<l[first];
for(IdxTy k=(first+1); k<f; ++k) ss<<sep<<l[k];
return ss.str();
} // sum

Myt& operator<<(const StrTy & word) 
{

	const IdxTy szl=m_lines.size();
	if (szl==0) m_lines.push_back(Line());
	m_lines.back().push_back(word);
	const IdxTy sz=m_lines.back().size();
	if (sz<m_min_words) m_min_words=sz;
	if (sz>m_max_words) m_max_words=sz;

return *this; 
}
Myt& operator<<(const IdxTy & word) {Ss ss; ss<<word; return (*this)<<ss.str(); }
Myt& operator<<(const int & word) {Ss ss; ss<<word; return (*this)<<ss.str(); }
Myt& operator<<(const double & word) {Ss ss; ss<<word; return  (*this)<<ss.str(); }

//////////////////////////////////////////////////
void load_from_string(const StrTy & s, const bool debug=false)
{
Ss ss;
ss<<s;
CommandInterpretter li(&ss);
ConmfigureReader( li);
li.use_alt_eol(',',true);
Load(li,debug);
}

void load(CommandInterpretter & li,const bool debug=false )
{
    while (li.nextok())
    {
    if (!li.ok()) break;
    const IdxTy sz=li.size();
    if ( m_ignore_hash) if (sz>0) if (li.words()[0].c_str()[0]=='#') continue;
    if (m_min_words>m_max_words) m_min_words=sz;
    Add(li.words());
    if (debug) { MM_ERR(MMPR4(li.line(),sz,(*this).size(),li.ok())) }
    } // while 

} // load (CommandInterpretter .. 

void load(std::istream & is,const bool debug=false )
{
// IdxTy lines=0;
// IdxTy skip=0;
//const  IdxTy interval=100000;
//const IdxTy w=2;
CommandInterpretter li(&is);
ConmfigureReader(  li,debug);
    while (li.nextok())
    {
	// TODO FIXME  putting ok() in nextok() IS D ASDFASDFASD
		if (!li.ok()) break;
        const IdxTy sz=li.size();
		if ( m_ignore_hash) if (sz>0) if (li.words()[0].c_str()[0]=='#') continue;
		if (m_min_words>m_max_words) m_min_words=sz;
		Add(li.words());
		//const Line & w=li.words();
		//MM_SZ_LOOP(j,w,wsz)
	//	{
	//		MM_ERR(MMPR2(j,w[j]))
//		}
//		if (sz<m_min_words) m_min_words=sz;
//		if (sz>m_max_words) m_max_words=sz;
//		m_lines.push_back(li.words()); 
		if (debug) { MM_ERR(MMPR4(li.line(),sz,(*this).size(),li.ok())) }
//        if (sz<=w){ ++skip;  continue; }
//        ++m[li.word(w)];
//        ++lines;
//        if ((lines%interval)==0)
//            {MM_STATUS(MMPR4(m.size(),lines,li.word(w),m[li.word(w)])<<MMPR(skip))}
//      if (li.line().size()!=0) ts.push_back(li.word(1)); 
    } // nextok()

} // load

void load_lines(std::istream & is,const bool debug=false )
{
CommandInterpretter li(&is);
li.readline_ok(false);
ConmfigureReader(  li,debug);
MM_ONCE(" changed fudding load_lines to use configure but may be fuxekd "<<debug,)
//    li.set_split(1,' ');
    while (li.nextok())
    {
	// TODO FIXME AK putting ok() in nextok() IS ED ASDFASDFASD
		if (!li.ok()) break;
        const IdxTy sz=li.size();
		if (m_min_words>m_max_words) m_min_words=sz;
		typedef std::vector<StrTy> W;
		W w;
		w.push_back(li.line());
		const IdxTy fudd=li.line().size();
		Add(w);
		//if (debug) { MM_ERR(MMPR3(li.line(),sz,(*this).size())) }
		if (debug) { MM_ERR(MMPR(fudd)<<MMPR4(li.line(),sz,(*this).size(),li.ok())) }
    } // nextok()

} // load_lines


//void add(const Line & w) { Add(w); } 


void remove_ci()
{
delete m_lip;
m_lip=0;
}
CommandInterpretter *  add_ci()
{
remove_ci();
m_lip= new CommandInterpretter();
return m_lip;
}
CommandInterpretter * set_split(const IdxTy split, const char c)
{
MM_ONCE(" not setting m_splitter wtf ",)
	if (m_lip==0) add_ci();
    (*m_lip).set_split(split,c);
	return m_lip;
}

void add(const StrTy & line)
{
if (false) { MM_ERR(" adding "<<MMPR(line)) } 
if (m_lip!=0)
{

(*m_lip).nextline(line.c_str());
Add((*m_lip).words());
return; 
}
CommandInterpretter li;
//ConmfigureReader(  li);
li.nextline(line.c_str());
Add(li.words());

}



void write_file(const StrTy & fn, const IdxTy  flags=0)  const
{
std::ofstream ofs(fn);
ofs<<dump(flags);
}

void save(const StrTy & fn, const IdxTy  flags, const StrTy & sep)  const
{
std::ofstream ofs(fn);
ofs<<dump(flags,sep);
}


void append(const StrTy & fn,  const IdxTy flags,const StrTy & sep )const 
{
std::ofstream ofs(fn,std::ios_base::app);
ofs<<dump(flags,sep);
}




void append(const StrTy & fn,  const IdxTy flags=3)const 
{
std::ofstream ofs(fn,std::ios_base::app);
ofs<<dump(flags);
}


void save(const StrTy & fn,  const IdxTy flags=3)const 
{ write_file(fn,flags); }
//void Save( std::ostream & os, Ragged & r,  const IdxTy flags=3)
//{  os<<r.dump(flags); }


void dump_sep(std::ostream & os,  const IdxTy  flags, const StrTy & sep)  const
{
//std::ofstream ofs(fn);
os<<dump(flags,sep);

}



void dump_os(std::ostream & os,  const IdxTy  flags=0)  const
{
//std::ofstream ofs(fn);
os<<dump(flags);

}

void string(std::ostream & os, const IdxTy flags) const
{
MM_SZ_LOOP(i,(*this),sz)
{
const Line & l=(*this)[i];
MM_LOOP(jj,l) { os<<" "<<(*jj); } 
os<<"\\\\"; 

} // i 


} // string 


template <class Tm> void index_field
( Tm & d , const IdxTy ifield, const IdxTy flags, const IdxTy mark) const
{
const IdxTy nidx=ifield;
const IdxTy minlen=nidx+1;
IdxTy i=0;
//const bool flag_dup=Bit(flags,0);
MM_LOOP(ii,m_lines)
{
const Line & line=(*ii);
if (line.size()<minlen) continue;
const StrTy & nm=line[nidx];
//const StrTy & v=line[vidx];
if (mark!=0) d[nm]|=mark;
else d[nm]=i;
++i;
//rwm.set( nm,v,true);
} // ii 

} // index_field


void to_map(ReadWriteMap & rwm) const
{
const IdxTy nidx=0;
const IdxTy vidx=1;
const IdxTy minlen=((nidx>vidx)?nidx:vidx)+1;
MM_LOOP(ii,m_lines)
{
const Line & line=(*ii);
//MM_ERR(MMPR2(line.size(),minlen))
if (line.size()<minlen) continue;
const StrTy & nm=line[nidx];
const StrTy & v=line[vidx];
//MM_ERR(MMPR2(nm,v))

rwm.set( nm,v,true);
}
}
// note that with repeated keys this messes up
// this is designed to concatenate mail headers 
//void to_map(ReadWriteMap & rwm,const IdxTy key) //const
void to_map(ReadWriteMap & rwm,const IdxTy key) const
{
const IdxTy nidx=key;
const IdxTy vidx=(key==0)?1:0;
const IdxTy minlen=((nidx>vidx)?nidx:vidx)+1;
StrTy last="";
MM_LOOP(ii,m_lines)
{
const Line & line=(*ii);
//MM_ERR(MMPR2(line.size(),minlen))
if (line.size()<minlen) continue;
const StrTy&   nm=line[nidx];
StrTy  v="";
StrTy sep="";
for ( IdxTy j=0; j<line.size(); ++j) if ( j!=key){if (line[j]!="\n") if (line[j]!="\r" )   v=v+sep+ line[j]; sep=StrTy(" "); } 

const IdxTy sss=v.length();
char vx[v.length()+2];
IdxTy ptr=0;
for(IdxTy j=0; j<sss; ++j)
{
char c=v.c_str()[j];
if (c!='\r') if (c!='\n') { vx[ptr]=c; ++ptr;}
}
vx[ptr]=0;
if  ( nm.length()>0) rwm[nm]=vx;
else rwm[last]=rwm[last]+StrTy(" ")+vx;
last=nm;
//MM_ERR(" setting key "<<MMPR2(nm,v))
//rwm.set( nm,v,true);
}
}
typedef std::map<StrTy,IdxTy> KeyMap;
void index_locations(KeyMap & km, const IdxTy keypos=0) const
{
const IdxTy sz=size();
for(IdxTy i=0; i<sz; ++i)
{
const Line & l=(*this)[i];
if (l.size()<=keypos) continue;
if ( l[keypos]!="") km[l[keypos]]=i; // unique line names 
}
} // index_locations
bool filled_in(Myt * d, const Myt& unknown, const IdxTy keypos=0 )
{
KeyMap us,them;
if (d!=0) { *d=(*this); } 
bool rc=false;
index_locations(us,keypos);
unknown.index_locations(them,keypos);
MM_LOOP(ii,us)
{
const IdxTy loc=(*ii).second;
auto jj=them.find((*ii).first);
if (jj!=them.end()) 
{
const IdxTy locthem=(*jj).second;
bool filled= ((*this)[loc].size()<unknown[locthem].size()); 
if (d==0) return true;
if (filled ) { rc=true;  (*d)[loc]=unknown[locthem]; } 
} // end
} // us 
return rc;
} // filled_in


bool need_esc(const char c ) const  { return ((c==' ')||(c=='"')); }  
bool need_esc(const char c,const char cs ) const  { return ((c==cs)||(c=='"')); }  

// TODO FIXME add a latex output option 
StrTy dump( const IdxTy flags=0,const StrTy &  s=" ") const
{
Ss ss;
IdxTy i=0;
const bool use_space=((flags&1)!=0);
const bool add_seq=((flags&2)==0);
const bool add_quote=((flags&4)!=0);
const bool add_escapes=((flags&8)!=0);
const bool use_s=((flags&16)!=0);
const bool debug_parse=((flags&32)!=0);

StrTy seps="|";
if (use_space) seps=" ";
if (use_s) seps=s;
// added 2020-12-27 to bring arb sep out
if (use_s&&use_space) seps=m_sep;
StrTy quote="\"";

const StrTy sep=seps;
const char cs=sep.c_str()[0];
MM_LOOP(ii,m_lines)
{
IdxTy field=0;
const Line & line=(*ii);
if (debug_parse) ss<<MMPR(line.size()<<" ");
bool add_sep=false;
if (add_seq) {add_sep=true;  ss<<i; }
IdxTy debug_field_count=0;
MM_LOOP(jj,line)
{
if (debug_parse) { ss<<" *("<<debug_field_count<<")* "; } 
// TODO FIXME  this should quote or escape where needed
if ( add_sep) ss<<sep;
if (add_quote) ss<<quote;
if(add_escapes)
{
const IdxTy l1= (*jj).length();
const char * p=(*jj).c_str();
for(IdxTy i=0; i<l1; ++i ) {const char c=p[i];if ( need_esc(c,cs)) { ss<<'\\'; }   ss<<c; } 
}
else { ss<<(*jj); }  // TODO FIXME XXX need to escape 


if (add_quote) ss<<quote;
++field;
add_sep=true;
++debug_field_count;
}
ss<<CRLF;
++i;
} 

return ss.str();
}

StrTy fix_str(const StrTy & w, const bool add_escapes, const IdxTy flags, const 
StrTy & sep)
{
const IdxTy sz=w.length();
const char * p=w.c_str();
const IdxTy bs=3*sz+3;
char c[bs];
const char s=sep.c_str()[0];
//IdxTy pc=0;
IdxTy pout=0;
bool needs_q=false;
const bool qif=((flags&1)!=0);
const bool qalways=((flags&2)!=0);
if (qif) { 
for(IdxTy i=0; i<sz; ++i) { if (p[i]==s) { needs_q=true; break; }} 
}
const bool quote=qalways||(qif&&needs_q);
if (quote) { c[pout]='"'; ++ pout; } 
for(IdxTy i=0; i<sz; ++i) { 

c[pout]=p[i]; ++pout;
if (pout>=bs) { MM_ERR(" too big corrupt") } 
} // i  
if (quote) { c[pout]='"'; ++ pout; } 
c[pout]=0; ++pout;
if (pout>=bs) { MM_ERR(" too big corrupt") } 
return StrTy(c); 
} // fix_str
// huge overhead for any number of lines.. 
StrTy reassemble_value(const Line & line,  const IdxTy flags )
{
const bool use_space=((flags&1)!=0);
const bool add_seq=((flags&2)==0);
if (add_seq) { MM_ERR(" add_seq does nothing") } 
const bool add_quote=((flags&4)!=0);
const bool add_escapes=((flags&8)!=0);
const IdxTy qflag=(add_quote?2:0)+(add_escapes?1:0);
StrTy seps=" ";
if (use_space) seps=" ";
StrTy quote="\"";

const StrTy sep=seps;
Ss ss;
bool add_sep=false;
const IdxTy sz= line.size();
for(IdxTy i=1; i<sz; ++i)
{
// TODO FIXME  this should quote or escape where needed
if ( add_sep) ss<<sep;
//if (add_quote) ss<<quote;
ss<<fix_str(line[i],add_escapes,qflag,seps);
//ss<<line[i]; // TODO FIXME XXX need to escape 
//if (add_quote) ss<<quote;
//++field;
add_sep=true;
} // for 
return ss.str();
}


line_index_map line_index(const IdxTy key)
{
line_index_map m;
MM_SZ_LOOP(i,m_lines,sz)
{
const Line & line=m_lines[i];
const IdxTy sz=line.size();
if ( sz>key) { m[line[key]].push_back(i); }
}
return m;
}


// create map from field key to field val 
xlate_map vkmap(const IdxTy key)
{
xlate_map m;
const IdxTy mf=key; // (key>val)?key:val;
MM_SZ_LOOP(i,m_lines,sz)
{
const Line & line=m_lines[i];
const IdxTy sz=line.size();
if (sz<=mf) continue;
for(IdxTy val=0; val<sz; ++val)  { m[line[val]]=line[key]; }
}
return m;
}


// create map from field key to field val 
xlate_map xlate_field_map(const IdxTy key, const IdxTy val)
{
xlate_map m;
const IdxTy mf=(key>val)?key:val;
MM_SZ_LOOP(i,m_lines,sz)
{
const Line & line=m_lines[i];
const IdxTy sz=line.size();
if ( sz>mf) { m[line[key]]=line[val]; }
}
return m;
}

xlate_double_map xlate_field_max_map(const IdxTy key, const IdxTy val) const
{
xlate_double_map m;
const IdxTy mf=(key>val)?key:val;
MM_SZ_LOOP(i,m_lines,sz)
{
const Line & line=m_lines[i];
const IdxTy sz=line.size();
if ( sz>mf) { 
const D x=atof(line[val].c_str());
auto ii=m.find(line[key]);
//MM_ERR(" CK "<<MMPR4(key,val,x,line[key]))
if (ii==m.end()) m[line[key]]=x; // (*ii).second=x;
else
{
if (x>(*ii).second)  m[line[key]]=x; // (*ii).second=x;
}
//k	m[line[key]]=line[val]; 


}
}
return m;
}


// the given field is the canonical value all others are synonym
// make map to look up synonyms and return canon 
xlate_map synonym_map(const IdxTy key, const StrTy & def="")
{
xlate_map m;
//const IdxTy mf=(key>val)?key:val;
MM_SZ_LOOP(i,m_lines,sz)
{
const Line & line=m_lines[i];
const IdxTy szll=line.size();
const StrTy canon=(key<szll)?line[key]:def;
// and mapping canon to canon is desired. 
MM_SZ_LOOP(j,line,szl)
{
 m[line[j]]=canon; 
} // j 
} // i 
return m;
}

// these strings are terrible memory wastes
// this whole thing is mostly static and
// could just user a read buffer with a lot
// of holes if it is reall dynamic. 
void append(const IdxTy & line, const StrTy & v)
{
m_lines[line].push_back(v);
}



template <class Ty> 
void  append_list (const std::vector<Ty> & v)
{
MM_LOOP(ii,v)
{
Line w;
(*ii).line(w);
add(w);
}

} // append_list


Myt & operator+=(const Myt & that )
{
// this need to check words sizes item by item as the
// aggregates may be invalid although really just check for -1 or 
// something 
MM_LOOP(ii,that.m_lines)
{ Add(*ii);
// m_lines.push_back((*ii));
//const IdxTy sz=(*ii).size();
//		if (sz<m_min_words) m_min_words=sz;
//		if (sz>m_max_words) m_max_words=sz;
//		m_lines.push_back(li.words()); 
}
return *this; 
}
void min_line_size(const IdxTy m)
{
Lines x;
MM_LOOP(ii,m_lines) { if ((*ii).size()>=m) x.push_back((*ii)); } 
m_lines=x;
}

void transpose(Myt & d) const
{
d.clear();
const IdxTy sz=m_lines.size();
if (sz<1) return; 
// first line is header
const Line & header=m_lines[0];
const IdxTy fields=header.size();
d.m_lines=Lines(fields);
MM_ERR(" transpose "<<MMPR2(sz,fields))
for(IdxTy i=0; i<fields; ++i)
{
MM_ERR(" ASCK  "<<MMPR2(i,sz))
Line x(sz);
for(IdxTy j=0; j<sz; ++j){
if (false) { MM_ERR(" ASDCASDCASCASDCASDCASD  "<<MMPR2(j,m_lines[j].size())) }

 x[j]=m_lines[j][i]; 

}
d.m_lines[i]=x;
}


}




private:
bool Bit(const IdxTy f, const IdxTy b) const { return ((1<<b)&f)!=0;}
void ConmfigureReader( CommandInterpretter & li,const bool debug=false)
{
li.readline_ok(false);
if (m_sep.length()>0)
    li.set_split(m_splitter,m_sep.c_str()[0]);
    else li.set_split(m_splitter,' ');
//MM_ERR(" AFK "<<MMPR3(m_splitter, m_sep,m_sep.length()))
if (debug) { MM_ERR(" AFK "<<MMPR3(m_splitter, m_sep,m_sep.length())) } 
}
void Load( CommandInterpretter & li, const bool debug)
{
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if ( m_ignore_hash) 
			if (sz>0) if (li.words()[0].c_str()[0]=='#') continue;
		if (m_min_words>m_max_words) m_min_words=sz;
		Add(li.words());
		if (debug) { MM_ERR(MMPR3(li.line(),sz,(*this).size())) }
    } // nextok()
} // Load

void Add( const Line & x)
{
	const IdxTy sz=x.size();
	if (sz<m_min_words) m_min_words=sz;
	if (sz>m_max_words) m_max_words=sz;
	m_lines.push_back(x); 
}

Lines m_lines;
IdxTy m_min_words, m_max_words;
StrTy m_sep;
IdxTy m_splitter;
bool m_ignore_hash;
CommandInterpretter* m_lip;
}; // mjm_ragged_table




//////////////////////////////////////////////////////////////////////////
// map integers to strings similar to tokenizer but designed to be defined 
// numbers  and supports synonyms although string tokenizer could do that too. 
// Forward maps integers to a single string
// Reverse key is a string mapping to vector of integers.
// a proper string is in the fwd table, synonyums are not.
// fwd[rev[s][i]]==s means that s is the proper string for that node

class mjm_index_lut 
{
//typedef  zymo_traits::Tr  Tr;
typedef collections_traits::Tr Tr;
typedef mjm_index_lut Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef StrTy KeyTy;
typedef IdxTy Indx;
typedef std::map<Indx,KeyTy> Fwd;
public:
typedef Fwd::const_iterator Citor;
typedef Fwd::iterator Itor;
typedef std::vector<Indx> RevVec;
//typedef std::map<KeyTy,Indx > Rev;
typedef std::map<KeyTy,RevVec > Rev;
typedef Rev::const_iterator CRitor;
typedef Rev::iterator Ritor;
//typedef std::map<StrTy, IdxTy> Txmap;
public:
typedef RevVec reverse_type;
typedef mjm_index_lut mapped_type;

bool proper(const IdxTy node, const StrTy & name)
{
auto ii=m_fwd.find(node);
if (ii==m_fwd.end()) return false;
if ((*ii).first==node) return true; 
return false;
}

Citor find(const Indx & k ) const { return m_fwd.find(k); } 
CRitor find(const KeyTy & k ) const { return m_rev.find(k); } 

// this has to work on rev only 

Citor fend( ) const { return m_fwd.end(); } 
Citor fbegin( ) const { return m_fwd.begin(); } 
Itor fend( )  { return m_fwd.end(); } 
//Itor fbegin( )  { MM_ONCE(" rev will not be consistent doh",) ;  return m_fwd.begin(); } 
Itor fbegin( )  { return m_fwd.begin(); } 


// this is not unique due to same name in diffeernt
// hiearchy need to search the tree not the name
CRitor rend( ) const { return m_rev.end(); } 
CRitor rbegin( ) const { return m_rev.begin(); } 
Ritor rbegin( )  { return m_rev.begin(); } 
Ritor rend( )  { return m_rev.end(); } 

#if ! LUT_LOAD_FWD 
//#warning using the reverse lut 
// the lut is loaded from the rev table shich is better 
CRitor end( ) const { return rend(); } 
CRitor begin( ) const { return rbegin(); } 
Ritor begin( )  { return rbegin(); } 
Ritor end( )  { return rend(); } 
#else
Citor end( ) const { return fend(); } 
Citor begin( ) const { return fbegin(); } 
Itor begin( )  { return fbegin(); } 
Itor end( )  { return fend(); } 

#endif





mjm_index_lut(): BADIDX(~0), BADKEY("!")  {}
IdxTy size() const { return m_fwd.size(); } 
static const RevVec & NulRev()
{
static RevVec x;
return x; 
}
//const Indx & operator()(const KeyTy & s ) const 

// return vector but do not modify lut
const RevVec & operator()(const KeyTy & s ) const 
{ 
	CRitor ii=m_rev.find(s); 
	if ( ii==m_rev.end())  { return NulRev(); } // BADIDX; }
	return (*ii).second;
}
// return key but do not add to maps 
const KeyTy  & operator()(const Indx  & i ) const 
{ 
	Citor ii=m_fwd.find(i); 
	if ( ii==m_fwd.end())  { return BADKEY; }
	return (*ii).second;
//	const IdxTy sz=m_rev.size();
//	if ((i<0)||(i>=sz)) { return BADKEY; } 	
//	return m_rev[i];


}

// use like regular map except that nothing created 
const KeyTy  & operator[](const Indx  & i ) const 
{ 
return (*this)(i); 
//	const IdxTy sz=m_rev.size();
//	if ((i<0)||(i>=sz)) { return BADKEY; } 	
//	return m_rev[i];
}



// define a key value pair in forward dirction,
// loc now points to string s 
// do not add to reverse lookup . 
void add_at(const KeyTy & s, const IdxTy & loc)
{
	Citor ii=m_fwd.find(loc); 
	if ( ii!=m_fwd.end())
	{
		if ((*ii).second!=s)
		{
			MM_ERR(" need dmel maping error "<<MMPR3(s,loc,(*ii).second))
		}
	}
	m_fwd[loc]=s;
//	while ((m_rev.size())<=loc) m_rev.push_back(StrTy());
//if ( m_rev[loc]!="") { if (m_rev[loc]==s) return ; else
//{ MM_ERR(" mapping error "<<MMPR3(s,loc,m_rev[loc]) )}}
//if (false) 	m_rev[s]=loc;
if (false) 	m_rev[s].push_back(loc);

//	m_rev.push_back(s);
}
// define a key value pair location maps to string s and
// add location to reverse table too for consistency 
void add_definition(const KeyTy & s, const IdxTy & loc)
{
	Citor ii=m_fwd.find(loc); 
	if ( ii!=m_fwd.end())
	{
		if ((*ii).second!=s)
		{
			MM_ERR(" need dmel maping error "<<MMPR3(s,loc,(*ii).second))
		//	MM_DMELF(f,"badfield", " badfieldi "<<MMPR(i))
		}
	}
	m_fwd[loc]=s;
// the authoritative version is always first 
//if (false) 	m_rev[s].insertpush_back(loc);
auto & rv= m_rev[s];
rv.insert(rv.begin(),loc);
//	m_rev.push_back(s);
}


// a new string points to a possibly alread defined 
// forward location that could be empty 
void add_synonym(const KeyTy & s, const IdxTy & loc)
{
// TODO check consistency etc 
//m_rev[s]=loc;
m_rev[s].push_back(loc);
}
// add all the existing forward values to rev,
// note that unique is not checked creating possible
// osddities. 
void add_fwd()
{
//MM_LOOP(ii,m_fwd) m_rev[(*ii).second]=(*ii).first; 
// FIXME check for dups at least 
MM_LOOP(ii,m_fwd) m_rev[(*ii).second].push_back((*ii).first); 

}
// for loading, only the fwd direction is indexed with add_at.
// add all of those to reverse and prepare for cinnamons later
void fix_fwd()
{
m_rev.clear();
//IdxTy collisions=0;
for (auto ii=fbegin(); ii!=fend(); ++ii)
{
const auto&  k=(*ii).second;
//if (m_rev.find(k)!=m_rev.end()) 
//{
//++collisions;
//MM_ONCE("do not use for lookup search hierarchy  collision "<<k,)
//}
//m_rev[k]=(*ii).first;
m_rev[k].push_back((*ii).first);

}

//MM_ERR(" DONOT USE FOR LOOKUPcollision count "<<collisions)
}
// clean up fwd table but synonyms and propers can get mixed up 
void fix_rev()
{
m_fwd.clear();
IdxTy collisions=0;
for (auto ii=rbegin(); ii!=rend(); ++ii)
{
const StrTy & str=(*ii).first;
const RevVec &  k=(*ii).second;
MM_LOOP(jj,k)
{
if (m_fwd.find(*jj)!=m_fwd.end()) 
{ ++collisions; MM_ONCE(" first  collision "<<str,) }
//m_rev[k]=(*ii).first;
//m_rev[k].push_back((*ii).first);
m_fwd[*jj]=str;
}
} // fix_rev
MM_ERR(" FWD collision count "<<collisions)
}


// 
static StrTy  header(IdxTy & i ) 
{
switch (i)
{
case 0: { return StrTy("id"); } 
case 1: { return StrTy("value"); } 
case 2: { return StrTy(""); } 
default:
MM_ERR(" bad case need DMEL "<<MMPR(i))
}
//return h;
return StrTy("");
} 

template <class Ti , class Tf> 
const StrTy   field(const Ti & ii, const IdxTy & i, Tf & f )  const
{
switch (i)
{
case 0: { return f((*ii).first); } 
case 1: { 

std::vector<StrTy> footnotes;
const auto & v= (*ii).second;
const auto & s= (*ii).first;
MM_SZ_LOOP(i,v,sz) 
{StrTy note="";
// this should call proper()
auto jj=m_fwd.find(v[i]);
 if ( jj!=m_fwd.end()) if (((*jj).second)==s)  note="*"; footnotes.push_back(note); }  
return f(v,footnotes);
//return f((*ii).second);


 } 
//case 1: { return m_parent; } 
//case 2: { return m_rank; } 
//case 3: { return m_src; } 
default:
//MM_ERR(" bad case need DMEL "<<MMPR(i))
MM_DMELF(f,"badfield", " badfieldi "<<MMPR(i))
		//	MM_DMELF(f,"badheder ",i);  
}
return f(BADIDX);
} 

template <class Ti , class Tf> 
void  field(const Ti & key, const IdxTy & i, Tf & f, const StrTy & w ) 
{
switch (i)
{
case 0: break ; 
case 1: { 
#if  LUT_LOAD_FWD 
//m_fwd[key]=w;
f(m_fwd[key],w);
#else

std::vector<StrTy> footnotes;
//f(m_rev[key],w);
auto & v=m_rev[key]; 
f(v,footnotes, w);
MM_SZ_LOOP(j,footnotes,sz) { if (footnotes[j].length()!=0) m_fwd[v[j]]=key; } 
#endif

break;
}
//return f((*ii).first); } 
//case 1: { return f((*ii).second); } 
//case 1: { return m_parent; } 
//case 2: { return m_rank; } 
//case 3: { return m_src; } 
default:
//MM_ERR(" bad case need DMEL "<<MMPR(i))
MM_DMELF(f,"badfield", " badfieldi "<<MMPR(i))
		//	MM_DMELF(f,"badheder ",i);  
}
//return f(BADIDX);
} 

// move itor to new string value, get the stuff in m_rev too 
void move(const StrTy & vd,const Itor & ii)
{
const IdxTy loc=(*ii).first;
const StrTy vorg=(*ii).second;
if (vd==vorg) return; 
m_fwd[loc]=vd;

auto& dvec=m_rev[vd];
auto& svec=m_rev[vorg];
MM_LOOP(jj,dvec) { if ((*jj)==loc) { MM_ERR(" danger will robinson "<<MMPR3(vorg,vd,loc) )  break; }}
MM_LOOP(jj,svec) { if ((*jj)==loc) { dvec.push_back(loc); svec.erase(jj); break; }}

if (svec.size()==0) m_rev.erase(m_rev.find(vorg));


}


// TODO too low level access need to fix all of this 
void rmove(const StrTy & to, const StrTy & from)
{
//auto rto=m_rev.find(to);
auto rfrom=m_rev.find(from);
m_rev[to]=(*rfrom).second;
m_rev.erase(rfrom);

}

bool contains(const IdxTy n, const RevVec & v)
{
MM_LOOP(ii,v) { if ((*ii)==n) return true; }
return false; 
} // contains
void add_rev_unique(const StrTy & name, const RevVec & v)
{
MM_LOOP(ii,v)
{
const IdxTy & i=(*ii);
RevVec & vnew= m_rev[name];
if (!contains(i,vnew)) vnew.push_back(i);

}

}
 
// make a new rev table from this one but using split words
Myt split_words()
{
Myt x;

MM_LOOP(ii,m_rev)
{
const StrTy & name=(*ii).first;
const RevVec & v=(*ii).second;
Ss ss(name);
while (true) 
{ StrTy y; ss>>y; if (y.length()==0) break; 
//MM_ERR(MMPR3(name, y,v.size()))
x.add_rev_unique(y,v); }

} // ii 
return x; 
} // split_words

void clear()
{
m_fwd.clear();
m_rev.clear();
}


private:


Fwd m_fwd;
Rev m_rev;
Indx BADIDX;
KeyTy BADKEY;

}; // mjm_index_lut 


/////////////////////////////////////////////////////////////////////////
class collate_order_map
{
//typedef  zymo_traits::Tr  Tr;
typedef collections_traits::Tr Tr;
typedef collate_order_map Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef StrTy KeyTy;
typedef IdxTy Indx;
// REVERSE convention of the mjm_index_lut.
// Forward maps strings to integers 
// The integer is the order or location for the string.
// This is 1:1 although really that need not apply for formatting uses.
// Ints should be consecutive hence use of vector inverse
typedef std::map<KeyTy, Indx> Fwd;
typedef Fwd::const_iterator Citor;
typedef std::vector<KeyTy > Rev;
typedef std::map<StrTy, IdxTy> Txmap;
public:
static const IdxTy & bad() { static const IdxTy b=~0U; return b; } 
Citor find(const StrTy & k ) const { return m_fwd.find(k); } 
// these are in alpha order NOT value order .... doh 
Citor end( ) const { return m_fwd.end(); } 
Citor begin( ) const { return m_fwd.begin(); } 

collate_order_map(): BADIDX(~0), BADKEY("!")  {}
IdxTy size() const { return m_rev.size(); } 

const Indx & operator()(const KeyTy & s ) const 
{ 
	Citor ii=m_fwd.find(s); 
	if ( ii==m_fwd.end())  { return BADIDX; }
	return (*ii).second;
}

const KeyTy  & operator()(const Indx  & i ) const 
{ 
	const IdxTy sz=m_rev.size();
	if ((i<0)||(i>=sz)) { return BADKEY; } 	
	return m_rev[i];
}

const KeyTy  & operator[](const Indx  & i ) const 
{ 
	const IdxTy sz=m_rev.size();
	if ((i<0)||(i>=sz)) { return BADKEY; } 	
	return m_rev[i];
}



void add_at(const KeyTy & s, const IdxTy & loc)
{
	Citor ii=m_fwd.find(s); 
	if ( ii!=m_fwd.end())
	{
		if ((*ii).second!=loc)
		{
			MM_ERR(" need dmel maping error "<<MMPR3(s,loc,(*ii).second))
		}
	}
	m_fwd[s]=loc;
	while ((m_rev.size())<=loc) m_rev.push_back(StrTy());
if ( m_rev[loc]!="") { if (m_rev[loc]==s) return ; else
{ MM_ERR(" mapping error "<<MMPR3(s,loc,m_rev[loc]) )}}
m_rev[loc]=s;

//	m_rev.push_back(s);
}


Indx add( const KeyTy & s)
{
	Citor ii=m_fwd.find(s); 
	if ( ii==m_fwd.end())  { m_fwd[s]=m_fwd.size(); m_rev.push_back(s);  }
return (*this)(s);
} 
// this map is ordered but trash the order in this case 
void append( const Txmap & m)
{
MM_LOOP(ii,m)
{
add((*ii).first);
}

}
// the map is likely in sorted order, need to check 
void sort()
{
std::sort(m_rev.begin(),m_rev.end());
m_fwd.clear();
MM_SZ_LOOP(i,m_rev,sz) { m_fwd[m_rev[i]]=i; }

}
void dump( std::ostream & os, const StrTy & label, const IdxTy & flags) const 
{
const StrTy sep=" ";
const bool dump_fwd=true;
if (dump_fwd)
{
MM_LOOP(ii,m_fwd)
{
os<<label<<sep<<(*ii).first<<sep<<(*ii).second<<CRLF; 
} // ii 

} // dump_fwd
const bool dump_rev=false;
if (dump_rev)
{
IdxTy i=0; 
MM_LOOP(ii,m_rev)
{
os<<label<<sep<<i<<sep<<(*ii)<<CRLF; 
++ii;
}  //ii
} // dump_r
}


bool assert_consistent()
{
IdxTy errc=0;
auto & x=*this; 
MM_LOOP(ii,m_fwd) 
{
if (x((*ii).first)!=((*ii).second))
	{ ++errc; MM_ERR(" bad map "<<MMPR4(errc, (*ii).first,(*ii).second,x((*ii).first))) }

}
IdxTy i=0;
MM_LOOP(ii,m_rev) 
{
if (x((*ii))!=(i))
	{ ++errc; MM_ERR(" bad vector  "<<MMPR4(errc, i,(*ii),x((*ii)))) }
++i;
}


return ( errc==0) ; 
}
private:


Fwd m_fwd;
Rev m_rev;
Indx BADIDX;
KeyTy BADKEY;

}; // collate_order_map




#endif




// based on ncbi taxonomy data files 
class mjm_tax_struct
{
typedef mjm_tax_struct Myt;
typedef collections_traits::Tr Tr;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef std::vector<IdxTy> Children;
typedef std::vector<StrTy> Header;
typedef std::vector<IdxTy> Fields;

public:
enum { BAD=~0 };
//static const IdxTy & bad() { const static IdxTy b=~0U; return b; } 
mjm_tax_struct(): m_parent(BAD),m_id(BAD),m_kids(0),m_rank(BAD),m_src(BAD) {}
void set(const IdxTy & p, const IdxTy & id ) { m_parent=p; m_id=id; } 
void set(const IdxTy & p, const IdxTy & id, const IdxTy & rank ) { m_parent=p; m_id=id; m_rank=rank; } 
void set(const IdxTy & p, const IdxTy & id, const IdxTy & rank, const IdxTy & src ) { m_parent=p; m_id=id; m_rank=rank;  m_src=src;} 
//TODO FIXME  was adding invalid self to children need to fi x

void add_child(const IdxTy & c) {if (c!=bad())  m_kids.push_back(c); } 
// 
void add_child(const IdxTy & id, const IdxTy & c) {m_id=id; if (c!=bad())  m_kids.push_back(c); } 
bool had_child(const IdxTy & id ) const
{
	MM_SZ_LOOP(i,m_kids,sz) { if (m_kids[i]==id) return true; } 
	return false;  
}

static IdxTy BADU() {  return IdxTy(BAD); } 
bool is_valid() const { return m_id!=BADU(); } 
const IdxTy bad() const { return BADU(); } 
const IdxTy  size() const { return m_kids.size(); } 
const IdxTy&  id() const { return m_id; } 
const IdxTy&  parent() const { return m_parent; } 
const bool  had_parent() const { return m_parent!=BADU(); } 
const IdxTy&  rank() const { return m_rank; } 
const IdxTy&  src() const { return m_src; } 
static StrTy  header(IdxTy & i ) 
{
switch (i)
{
case 0: { return StrTy("id"); } 
case 1: { return StrTy("parent"); } 
case 2: { return StrTy("rank"); } 
case 3: { return StrTy("src"); } 
case 4: { return StrTy(""); } 
default:
MM_ERR(" bad case need DMEL "<<MMPR(i))
}
//return h;
return StrTy("");
} 



template < class Tf> 
void  vv( const IdxTy & key, const IdxTy & i, const IdxTy & v, Tf & f )  const
{
if (v==IdxTy(BAD)) {
MM_DMELF(f,"badfield", "bad field "<<MMPR3(i,key,v))
}
}


template < class Tf> 
const StrTy  field(const IdxTy&  key, const IdxTy & i, Tf & f )  const
{ // ??? 
//const IdxTy key=f.myatoi(_key);
switch (i)
{
// otherwise the unassigned as huge numbers and file size is huge 
case 0: { vv(key,i,m_id,f);  return f(int(m_id)); } 
case 1: { vv(key,i,m_parent,f); return f(int(m_parent)); } 
case 2: { vv(key,i,m_rank,f); return f(int(m_rank)); } 
case 3: { return f(int(m_src)); } 
default:
//MM_ERR(" bad case need DMEL "<<MMPR(i))
MM_DMELF(f,"badcase", "bad case "<<MMPR(i))
}
return f(BAD);
} 
template < class Tf> 
void   field( const IdxTy & i, Tf & f, const StrTy & w  )  
{
switch (i)
{
// otherwise the unassigned as huge numbers and file size is huge 
case 0: {  f(m_id,w); break;  } 
// TODO if there is no parent definition, this load will not
// create one as opposed to the ncbi loading. 
case 1: {  f(m_parent,w); break; } 
case 2: {  f(m_rank,w);break; } 
case 3: {  f(m_src,w); break; } 
default:
//MM_ERR(" bad case need DMEL "<<MMPR(i))
MM_DMELF(f,"badcase", "bad case "<<MMPR(i))
}
} 






template <class Tl > void dump( OsTy & os, const Tl & lut, const IdxTy flags=0) const
{
const StrTy sep=" ";
const bool dump_number=((flags&1)!=0);
MM_LOOP(ii,m_kids)
{
os<<sep;
const auto & i=(*ii);
auto jj=lut.find(i);
if (jj!=lut.fend()) os<<(*jj).second;  else os<<".";
if (dump_number) os<<"("<<i<<")"; 
}

}
// fudding subscripts are not fudding const fudd 
template <class Tl > StrTy tax_string( Tl & lut, const string_tokenizer & st, const IdxTy & flags=0 ) const
{
Ss ss;
//ss<< MMPR4(m_id,lut[m_id],m_parent,lut[m_parent])<<MMPR(int(m_rank)); 
ss<< MMPR4(m_id,lut[m_id],int(m_parent),lut[m_parent])<<MMPR(st(m_rank)); 
return ss.str();
} // tax_string 
const Children & kids() const { return m_kids; } 
IdxTy m_parent;
IdxTy m_id; // could be implicit store this in a map or something. 
Children m_kids;
IdxTy m_rank;
IdxTy m_src;
};  // mjm_tax_Struct


class mjm_tax_map : public std::map<collections_traits::Tr::IdxTy, mjm_tax_struct> 
 {

public:
template <class Tu> 
static StrTy  header(Tu & i )  { return mjm_tax_struct::header(i); } 

template <class Ti, class Tu,class Tf> 
const StrTy  field(const Ti & ii, const Tu & i, Tf & f )  const
{
return (*ii).second.field((*ii).first,i,f);
} 

template <class Ti, class Tu,class Tf> 
void  field(const Ti & _key, const Tu & i, Tf & f, const StrTy & s ) 
{
const unsigned int  & key=f.myatoi(_key); 
//return (*find(kkey)).second.field(i,f,s);
// TODO this needs to creaet parent and child pointers etc... 
// #error make family pointers 
if (i==1){ unsigned int  p=0; f(p,s);  (*this)[p].add_child(p,(*this)[key].id()); }
return (*this)[key].field(i,f,s);
}


}; // mjm_tax_map

//template <class Tobj, class Traits > class mjm_collection 
class mjm_tax_tree : public mjm_collection< mjm_tax_map, collections_traits::Tr>
{
//typedef collections_traits::Tr Tr;
typedef mjm_tax_tree Myt;
typedef mjm_tax_map MyMap;
typedef mjm_tax_struct MyStruct;
typedef mjm_collection< mjm_tax_map, collections_traits::Tr> Super;
typedef Super::Tr Tr; 
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef string_tokenizer St;

public:
//void clear() {(*this) = Myt(); } 
IdxTy size() const { return m_lut.size(); } 
template < class Ti, class Tf> 
const StrTy  field( const Ti & ii, const IdxTy & i, Tf & f )  const
{
return (*ii).second.field(i,f); 
}
template < class Ti, class Tf> 
void  field( const Ti & ii, const IdxTy & i, Tf & f, const StrTy & s )  
{
 (*ii).second.field(i,f,s); 
}


void read_ncbi_dmp(const StrTy & fn)
{
read_ncbi_nodes_dmp(fn+StrTy("/nodes.dmp"));
read_ncbi_names_dmp(fn+StrTy("/names.dmp"));

}
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245000/
void read_ncbi_nodes_dmp(const StrTy & fn) { std::ifstream  is(fn); read_ncbi_nodes_dmp( is); }
void read_ncbi_names_dmp(const StrTy & fn) { std::ifstream  is(fn); read_ncbi_names_dmp( is); }

/*
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt

nodes.dmp
---------

This file represents taxonomy nodes. The description for each node includes 
the following fields:

	tax_id					-- node id in GenBank taxonomy database
 	parent tax_id				-- parent node id in GenBank taxonomy database
 	rank					-- rank of this node (superkingdom, kingdom, ...) 
 	embl code				-- locus-name prefix; not unique
 	division id				-- see division.dmp file
 	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	genetic code id				-- see gencode.dmp file
 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	mitochondrial genetic code id		-- see gencode.dmp file
 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	comments				-- free-text comments and citations
*/

void read_ncbi_nodes_dmp(IsTy & is)
{
const IdxTy id_idx=0;
const IdxTy parent_idx=1;
const IdxTy rank_idx=2;
const IdxTy div_idx=4;
const IdxTy min_sz=5;

	CommandInterpretter li(&is);
	li.set_split(1,'|');
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if (sz<min_sz)
		{
			//DMel   
			//DMel("tooshort",m_ss<<MM_STR_LOC<<li.line());  
			MM_DMEL("tooshort",li.line());  
			continue;
		}
		const auto & w=li.words();
		const IdxTy div=myatoi(w[div_idx]);
		if (div!=0 ) continue;
		const IdxTy id=myatoi(Trim(w[id_idx]));
		const IdxTy parent=myatoi(Trim(w[parent_idx]));
		const StrTy & rank=Trim(w[rank_idx]);
		//m_obj[id].set(parent,id);
		m_obj[id].set(parent,id,m_st(rank));
		// if parent has not been define, the id is invliad
		m_obj[parent].add_child(parent,id);
    } // li
} // read_ncbi_nodes_dmp
/*

ames.dmp
---------
Taxonomy names file has these fields:

	tax_id					-- the id of node associated with this name
	name_txt				-- name itself
	unique name				-- the unique variant of this name if name not unique
	name class				-- (synonym, common name, ...)

*/

#if 0 

cat names.dmp  | sed -e 's/^[^|]*|[^|]*|[^|]*|//' | sed -e 's/|.*//g' | sed -e 's/[^a-z]//g' | sort | uniq -c | more
   1161 acronym
    303 anamorph
 407210 authority
    229 blastname
  14196 commonname
  24978 equivalentname
    486 genbankacronym
    107 genbankanamorph
  28135 genbankcommonname
   2959 genbanksynonym
  35988 includes
    628 inpart
   1386 misnomer
  35103 misspelling
1677153 scientificname
 167035 synonym
    179 teleomorph
 113616 typematerial
#endif



void read_ncbi_names_dmp(IsTy & is)
{
const IdxTy id_idx=0;
const IdxTy name_idx=1;
//const IdxTy rank_idx=2;
const IdxTy class_idx=3;
const IdxTy min_sz=4;
const bool adding_syn=true; // synonyms
	CommandInterpretter li(&is);
	li.set_split(1,'|');
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if (sz<min_sz)
		{
			if (sz==0) continue; 
			//DMel("tooshort",m_ss<<MM_STR_LOC<<li.line());  
			MM_DMEL("tooshort",li.line());  
			continue;
		}
		const auto & w=(li.words());
		const IdxTy id=myatoi(w[id_idx]);
		// this does look up by id and does not care if obj is valid
		if (m_obj.find(id)==m_obj.end() ) continue;
		//StrTy name=mjm_strings::to_lower(Trim(w[name_idx]));
		StrTy name=(Trim(w[name_idx]));
		const StrTy ty=Trim(w[class_idx]);
//		if ( ty=="scientific name") m_lut[id]=name;
		if (!adding_syn){  if ( ty=="scientific name") m_lut.add_at(name,id); } 
		else
		{
			name=mjm_strings::to_lower(name);
			name=BraceRemove(name);
			if ( ty=="scientific name") m_lut.add_definition(name,id);
			else if ( ty=="equivalent name") m_lut.add_synonym(name,id);
			else if ( ty=="common name") m_lut.add_synonym(name,id);
			else if ( ty=="synonym") m_lut.add_synonym(name,id);
			else if ( ty=="includes") m_lut.add_synonym(name,id);
			else if ( ty=="acronym") m_lut.add_synonym(name,id);
			else if ( ty=="genbank common name") m_lut.add_synonym(name,id);
			else if ( ty=="genbank synonym") m_lut.add_synonym(name,id);
			else if ( ty=="genbank acronym") m_lut.add_synonym(name,id);

		}
		//	else { MM_ERR(" class "<<MMPR(ty)) }
    } // li
//	m_lut[~0]="BAD";
	if (!adding_syn) m_lut.add_at("BAD",~0);
	else m_lut.add_definition("BAD",~0);
	//note that this is easy to get inconsistent and must
	// be fixed in this order, don't hand out non-const itor doh 
	// FIXME 
	fix_lut(adding_syn);
	if (!adding_syn) m_lut.fix_fwd();
} // read_ncbi_nodes_dmp


void dump(const StrTy & fn)
{
check_valid();
OsTy & os=std::cerr;
MM_LOOP(ii,m_obj)
{
const IdxTy id=(*ii).first;
const MyStruct & o=(*ii).second;
MM_ERR(MMPR4(id,m_lut[id],m_lut[o.id()],o.size()))
o.dump( os, m_lut);  os<<CRLF;

}

}

void check_valid()
{
IdxTy i=0;
MM_LOOP(ii,m_obj)
{
if (!(*ii).second.is_valid())
{
MM_ERR(" invalid "<<MMPR4(m_lut[(*ii).first], (*ii).first,(*ii).second.id(), (*ii).second.size()))
}
++i; 
}
MM_ERR(" total entries "<<i)
} //check_valid








void read_composite(const StrTy  & fn)
{
std::ifstream is(fn);
read_composite(is);
}

void read_composite(IsTy & is)
{
form_def<> f;
   CommandInterpretter li(&is);
f.setup(li);
if (!li.nextok()) return; 
read_map("tax-con",m_obj,f,li);
//write_map(os,(*this),f);
//St m_st; //mjm_index_lut  m_lut;
read_map("tax-lut",m_lut,f,li);
// note that the rev lut needs to be made after 
// fwd is const and then add synonyms. 

#if LUT_LOAD_FWD
m_lut.fix_fwd();
#else
//m_lut.fix_rev();
#endif

read_map("tax-tok",m_st,f,li);
MM_ERR(" done reading map composite") 
}

template <class Tmm, class Tf, class Tn > void read_map( const Tn & name,  Tmm & m , Tf & f, CommandInterpretter & li)
{
const StrTy sep=f.sep(); // " ";
typedef std::map<StrTy, IdxTy> Hm;
typedef std::map<StrTy, IdxTy> Om;
typedef std::map<IdxTy, IdxTy> Xm;
Hm h;
Om mh;
Xm xm;
IdxTy nfields=0;
bool debug=false; // (name==StrTy("tax-tok"));
MM_ERR(" starting "<<name)
while (true)
{
StrTy s= Tmm::mapped_type::header(nfields);
if (s.length()==0) break;
// doh, the aray lookup adds the fing element... 
h[s]=h.size()-1;
++nfields;
}
// should point to name line
if (li.line()!=name)
{
	MM_DMEL("maploadnamewronf",li.line()<<" not equal "<<name);  

}
if  (!(li.nextok()) ) { MM_ERR("blank data "<<name) return; } 
 MM_SZ_LOOP(i,li,wsz) {mh[li.word(i)]=i; }  
MM_LOOP(ii,h) { const IdxTy i=(*ii).second; const IdxTy j=mh[(*ii).first];  xm[j]=i; } 
const IdxTy fsz=nfields;
IdxTy cnt=0;
if (mh.size()!=fsz) 
	{ MM_ERR(" spurious elements maybe "<<MMPR2(mh.size(),fsz)) } 
if (debug) { MM_ERR(" fii ") } 
    while (li.nextok())
    {
	++cnt;
if (debug) { MM_ERR(" fii ") } 
		const IdxTy sz=li.size();
		if (fsz!=sz)
		{
			MM_ERR(" ending " <<name<<"  with "<<MMPR(cnt)<<MMPR4(fsz,sz,wsz,li.line()))
			 return;
		}
		const auto & w=li.words();
if (debug) { MM_ERR(" fii ") } 
		//const IdxTy key=myatoi(w[0]);
		const StrTy key=(w[0]);
		//`ZZauto & o=m[key];
//		const auto & ii=m.find(key);
		for(IdxTy i=0; i<nfields; ++i){  
if (debug)		MM_ERR(MMPR3(i,xm[i],w.size()))
			m.field(key,i,f, w[xm[i]]); }

	} // li 
MM_ERR(" hit end of file on  "<<name)

} // read_map

void write_composite(const StrTy & fn)
{
if (fn.length()==0) { write_composite(std::cout); return; } 
std::ofstream os(fn);
write_composite(os);

}

void add_synonym(const StrTy & s, const IdxTy & loc)
{m_lut.add_synonym(s,loc); } 
// see also fix_lut??? 
#if 0
void conform(const IdxTy & flags )
{
MM_ERR(" do not call this crap ") 
const bool to_lower=true;
const bool make_single=true;
// this is over the REVERSE or it would b
MM_LOOP(ii,m_lut)
{
    StrTy  r=(*ii).first ; // =s;
    if ( to_lower) r=mjm_strings::to_lower(r);
    if ( make_single) mjm_strings::make_single(r);
}  // m_lut

//m_lut.fix_fwd();


} // conform
#endif



void write_composite(OsTy & os)
{
form_def<> f;
write_map(os,"tax-con",m_obj,f);
//write_map(os,(*this),f);
//St m_st; //mjm_index_lut  m_lut;
write_map(os,"tax-lut",m_lut,f);
write_map(os,"tax-tok",m_st,f);
}

template <class Tmm, class Tf, class Tn > void write_map(OsTy & os, const Tn & name, const Tmm & m , Tf & f)
{
const StrTy sep=f.sep(); // " ";
IdxTy nfields=0;
os<<name<<CRLF;
while (true)
{
StrTy s= Tmm::mapped_type::header(nfields);
if (s.length()==0) break;
if (nfields!=0) os<<sep;
os<<s;
++nfields;
}
os<<CRLF;
MM_LOOP(ii,m)
{
	// TODO figure out what this should be ... 	
//	const auto & o=(*ii).second;
	//const typename Tmm::mapped_type & o=(*ii).second;
	if (nfields>0) { os<<m.field(ii,0,f); } 
	for(IdxTy i=1; i<nfields; ++i) os<<sep<<m.field(ii,i,f);
	os<<CRLF;
}
} // write_map


void write_single(const StrTy &  fn) 
	{if (fn.length()==0) { write_single(std::cout); return ; } std::ofstream os(fn); write_single(os); }
void write_single(OsTy & os)
{
MM_LOOP(ii,m_obj)
{
//template <class Tl > StrTy tax_string(const Tl & lut, const IdxTy & flags=0) const
//const IdxTy id=(*ii).first;
const MyStruct & o=(*ii).second;
os<<o.tax_string(m_lut,m_st,0)<<CRLF;
//MM_ERR(MMPR4(id,m_lut[id],m_lut[o.id()],o.size()))
} // ii

} // write_single


bool  basic_operations(const StrTy & cmd, const StrTy & p1)
{
Myt & tt=(*this);
do { 
if (cmd=="load-tax-nodes") { tt.read_ncbi_nodes_dmp(p1);  continue; }
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") {tt.clear();  tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }
MM_ERR("no command fits for tax tree "<<MMPR2(cmd,p1))
return false;
} while (false);
return true;
}




StrTy lineage( const IdxTy & taxid  ) const
{
std::map<IdxTy, IdxTy> m;
return lineage(taxid,m);
}

template <class Tm> 
StrTy lineage( const IdxTy & taxid, const Tm & m  ) const
{
	StrTy x; 
	const StrTy sep=":";
	IdxTy id=taxid;
	const bool check_map=(m.size()!=0);
	while (true)
	{
		CItor ii= m_obj.find(id);
		if (ii==m_obj.end()) return x;	
		const MyStruct & op=((*ii).second); 
		if (!op.is_valid()) return x; 
		const bool add= (( check_map&&(m.find(op.rank())!=m.end() ))||!check_map) ;
		if (add) { x=m_lut[op.id()]+sep+x; }
		id=(op).parent();
	}

}
//template <class Tm >
void lineage( std::vector<StrTy> & vd, const IdxTy & taxid  ) const
{
std::map<IdxTy,IdxTy > m;
lineage(vd,taxid,m);

}

void lineage( std::vector<IdxTy> & vd, const IdxTy & taxid  ) const
{
	IdxTy id=taxid;
	while (true)
	{
		CItor ii= m_obj.find(id);
		if (ii==m_obj.end()) return ;	
		const MyStruct & op=((*ii).second); 
		if (!op.is_valid()) return ; 
		vd.push_back(op.id()); 
		id=(op).parent();
	}

}

template <class Tm >
void lineage( std::vector<StrTy> & vd, const IdxTy & taxid, const Tm & m  ) const
{
	StrTy x; 
	IdxTy id=taxid;
	const bool check_map=(m.size()!=0);
	while (true)
	{
		CItor ii= m_obj.find(id);
		if (ii==m_obj.end()) return ;	
		const MyStruct & op=((*ii).second); 
		if (!op.is_valid()) return ; 
//		const bool add= (( check_map&&(m.find(op.rank())!=m.end() ))||!check_map) ;
//		if (add) { x=m_lut[op.id()]+sep+x; }
		const bool add= (( check_map&&(m.find(op.rank())!=m.end() ))||!check_map) ;
		if (add) vd.push_back(m_lut[op.id()]); 
		id=(op).parent();
	}
}
IdxTy depth(const IdxTy & taxid) const
{
	IdxTy d=0;
	IdxTy id=taxid;
	while (true)
	{
		CItor ii= m_obj.find(id);
		if (ii==m_obj.end()) return d;	
		// TODO FIXME this is not real consistent on BAD parent depth 
		const MyStruct & op=((*ii).second); 
		if (!op.had_parent()) return d;
		id=(op).parent();
		++d;
	}

return d;
}



IdxTy root(const IdxTy & taxid) const
{
	IdxTy id=taxid;
	while (true)
	{
		CItor ii= m_obj.find(id);
		if (ii==m_obj.end()) return id;	
		const MyStruct & op=((*ii).second); 
		if (!op.had_parent()) return id;
		id=(op).parent();
	}

return id;
}


IdxTy root() const
{ 
auto ii=m_obj.begin();
if (ii==m_obj.end()){ MM_ERR(" tax tree empty")  return  ~0; } 
IdxTy rkey=(*ii).second.id();
IdxTy rootid=root(rkey);
MM_ERR( " found root "<<MMPR4((*ii).first,rkey,rootid,m_lut(rootid)))
return rootid;
}

std::vector<IdxTy> roots(const IdxTy & flags=0) // const
{
const bool print_found=((flags&1)!=0);
std::vector<IdxTy> r;
//auto ii=m_obj.begin();
//if (ii==m_obj.end()){ MM_ERR(" tax tree empty")  return  r; } 
MM_LOOP(ii,m_obj)
{
IdxTy rkey=(*ii).second.id();
IdxTy rootid=root(rkey);
bool found=false;
MM_LOOP(jj,r) { if (rootid==(*jj)) {found=true; break; }}
if (!found){ r.push_back(rootid);
if (print_found){ MM_ERR( " found root "<<MMPR3(rootid,node_name(rootid),node_info(rootid,256))) }} 
}

if (print_found){ MM_ERR( " found roots "<<MMPR(r.size()))}

return r;
}


std::vector<IdxTy> roots_unify() 
{
std::vector<IdxTy> v;
v=roots();
if (v.size()==1) return v;
MM_ERR(" adding dummy node to unify ")
std::vector<StrTy> vnew;
vnew.push_back(StrTy("FakeRoot"));
vnew.push_back(StrTy("Dummy"));
MM_LOOP(ii,v)
{
MM_ERR(" exisitng root "<<MMPR((*ii)))
vnew[1]=node_name((*ii));
add_new_nodes(vnew, 0);
}
return roots();
}

/*
 310 class
      3 cohort
   8734 family
    494 forma
  84443 genus
     16 infraclass
    102 infraorder
      3 kingdom
 224865 norank
   1439 order
     10 parvorder
    220 phylum
1319049 species
    466 speciesgroup
    133 speciessubgroup
    138 subclass
   2829 subfamily
   1388 subgenus
      1 subkingdom
    346 suborder
     28 subphylum
  21614 subspecies
    488 subtribe
      4 superclass
    816 superfamily
      5 superkingdom
     49 superorder
      2 superphylum
   2009 tribe
   7149 varietas


      3 kingdom
    220 phylum
 	310 class
   1439 order
   8734 family
  84443 genus
 224865 norank
1319049 species

*/
template <class Tm>
void normal_taxes(Tm & m, St & st)
{
m[st("superkingdom")]=1;
m[st("phylum")]=1;
m[st("class")]=1;
m[st("order")]=1;
m[st("family")]=1;
m[st("genus")]=1;
m[st("norank")]=1;
m[st("species")]=1;
}
bool is_descendant(const IdxTy& higher, const IdxTy & lower) const
{
IdxTy node=lower;
while (node!=IdxTy(~0))
{
if (node==higher) return true; 
auto ii=m_obj.find(node);
if (ii==m_obj.end()) return false;
const MyStruct & o=(*ii).second; 
node=o.parent();

}

return false;

}

StrTy  node_info( const IdxTy n, const IdxTy flags=0)
{
Ss ss;
node_info(ss,n,flags);
return ss.str();
}

static IdxTy & bad()
{
static IdxTy b=~0U;
return b; 

}

void node_info(std::ostream & os, const IdxTy n, const IdxTy flags=0)
{
const IdxTy dump_flags=flags;
const IdxTy string_flags=flags;
const bool dump_full=((flags&(1<<8))!=0);
MM_ERR("node_info "<<MMPR2(flags,dump_full))
//typedef mjm_index_lut::reverse_type Rt;
St & st= m_st;
const auto & node=m_obj[n];
if (dump_full)
{
 auto node=m_obj[n];
StrTy name=m_lut(n);
IdxTy uprank=0;
while (node.parent()!=bad())
{
const StrTy name=m_lut(node.id());
os<<MMPR4(uprank,st(node.rank()),node.id(),name)<<CRLF;
node=m_obj[node.parent()];
++uprank;
}
node=m_obj[n];
const std::vector<StrTy> c=kids(n);
MM_SZ_LOOP(i,c,sz) {if (( i%3)==0){ if (i!=0)  os<<CRLF; }   os<<MMPR2(i,c[i])<<" \t"; }
os<<CRLF;
return;
}
//const IdxTy parent=node.parent();
os<<node.tax_string( m_lut, st, string_flags);
node.dump( os, m_lut,dump_flags);

}
IdxTy add_generations(const IdxTy node, const IdxTy nlevels, const StrTy & name, const IdxTy flags=0)
{
	IdxTy parent=node;
	IdxTy id=m_lut.size();
St & st= m_st;
	for(IdxTy i=0; i<nlevels; ++i)
{	
	while ( m_lut.find(id)!=m_lut.fend()) ++id; 
	m_lut.add_definition(name,id);
//	const IdxTy fakefudd=st(StrTy("fake"));
	// rank is not used fudd 
	m_obj[id].set(parent,id,0,st(StrTy("fake")));
	//MM_ERR(" ASSSCK SHCK "<<MMPR3(parent,id,(fakefudd))) 
	// if parent has not been define, the id is invliad
	m_obj[parent].add_child(parent,id);
	parent=id;
	}
return id;
}

IdxTy  add_new_nodes(const std::vector<StrTy> &  vin, const IdxTy flags=0)
{
typedef mjm_index_lut::reverse_type Rt;
St & st= m_st;
IdxTy parent= ~0;
// TODO danger will robinson defined later too lol 
//IdxTy rank=0;
typedef std::vector<StrTy> Ranks;
Ranks ranks;
IdxTy id= ~0;
MM_SZ_LOOP(i,vin,sz)
{
	if (i>=ranks.size())
	{
		Ss ss;
		ss<<"level"<<i;
		ranks.push_back(ss.str()); 
	}
	const StrTy & name=vin[i];
	const IdxTy rank=i;
	id= ~0;
	const Rt & locs=m_lut(name);
	// this needs to insure it is from the rigth PARENT or
	// else we have noway to know when to add.. 
	bool missing=true;
	MM_SZ_LOOP(j,locs,locsz)
	{ 
		if (m_obj[locs[j]].parent()==parent) 
			{ missing=false; id=locs[j]; break; } 
	} // j 
	if (missing)
	{
		id=m_lut.size();
		while ( m_lut.find(id)!=m_lut.fend()) ++id; 
		m_lut.add_definition(name,id);
		m_obj[id].set(parent,id,st(ranks[rank]));
		// if parent has not been define, the id is invliad
		m_obj[parent].add_child(parent,id);
	}
	parent=id;
} // i 
return id;
} // add_new_nodes

void descend(const std::vector<StrTy> &  vin, const IdxTy flags=0)
{
typedef mjm_index_lut::reverse_type Rt;
St & st= m_st;
std::map<IdxTy,IdxTy> mtax;
normal_taxes(mtax,st);
const IdxTy bitflags=0;
IdxTy bit=best_indexed_taxon(  vin, bitflags);
std::map<IdxTy, IdxTy> m;
MM_SZ_LOOP(i,vin,sz)
{
const StrTy & name=vin[i];
const Rt & locs=m_lut(name);
MM_SZ_LOOP(j,locs,szlocs)
{
const IdxTy node=locs[j];
const IdxTy parent=m_obj[node].parent();
const bool proper=m_lut.proper(node,name);
MM_MSG(MMPR2(i,sz)<<MMPR3(vin[i],node,int(parent))<<MMPR3(proper,m_lut[node], m_lut[parent])<<MMPR2(bit,m_lut[bit]))
} // j 
if (szlocs==0) 
{
StrTy kids="";
StrTy lineae="";
if (i!=0)
{
const Rt & locso=m_lut(vin[i-1]);
MM_SZ_LOOP(j,locso,szlocs)
{
auto & childs=m_obj[locso[j]].kids(); // m_kids
MM_SZ_LOOP(j,childs,ksz) { kids+=StrTy(" ")+m_lut(childs[j]); }
} // j
for (IdxTy pidx=0; pidx<=szlocs; ++pidx)
{
lineae+=" & ";
std::vector<StrTy> lncbi;
if ( pidx<szlocs) { lineage(lncbi,locso[pidx],mtax); } 
else { if ( szlocs>0) break; } 
const IdxTy sz1=lncbi.size();
const IdxTy szx=(sz1>sz)?sz1:sz;
// lineae=lineage(locso[0],m);
for (IdxTy k=0; k<szx;++k)
{
lineae+=":";
StrTy other="";
if (k<sz) { other=vin[k]; const IdxTy nncbi=m_lut(other).size(); lineae+=other; if (nncbi==0) lineae+="!"; } 
lineae+=",";
if (k<sz1) { StrTy ncbi=lncbi[sz1-1-k]; if ( ncbi==other) lineae+="="; else lineae+=ncbi; } 
} // k
} // szlocs
} // i!=0
MM_MSG( " no hits some choices  "<<MMPR2(i,sz)<<MMPR3(vin[i],lineae,kids))
}
} // i 
} // descend







template <class Tv> StrTy best_indexed_taxon_str( const Tv & vin, const IdxTy flags)
{
Ss ss;
IdxTy id=best_indexed_taxon(vin,flags);
StrTy name=m_lut(id);
ss<<name<<"_"<<id;

return ss.str();
}


template <class Tv> IdxTy best_indexed_taxon( const Tv & vin, const IdxTy flags)
const 
{
const bool notice_singles=!false;
const bool notice_family=!false;
const bool notice_dmel=false;

typedef mjm_index_lut::reverse_type Rt;

IdxTy best=~0;
const IdxTy sz=vin.size();
if (sz==0) return best;
for (IdxTy i=sz-1;i<sz; --i) 
{
const Rt & locs=m_lut(vin[i]);
// should check descendents etc but ok for now 
const IdxTy lsz=locs.size();
MM_ERR(" fuicking tree search "<<MMPR2(lsz,vin[i]))
if (lsz == 1) { 
if (notice_singles) { MM_ERR(" single "<<MMPR2(locs[0],m_lut(locs[0]))) } 

return locs[0]; } 
if ( lsz==0) continue;
// many choices, check tree... 
 // recursive check 
//for (IdxTy j=i-1;i<sz; --j) 
for (IdxTy j=i;j<sz; --j) 
{
const Rt & locj=m_lut(vin[j]);
const IdxTy jsz=locj.size();
if (jsz==0) continue;

for (IdxTy ii=0; ii<lsz; ++ii)
for (IdxTy jj=0; jj<jsz; ++jj)
if (is_descendant(locj[jj],locs[ii]))
{
if (notice_family ) { MM_ERR(" family "<<MMPR3(locs[ii],m_lut(locj[jj]),m_lut(locs[ii]))) } 
 return locs[ii];
}

} // j 
// have no idea what to do here
if (notice_dmel) { MM_ERR(" need DMEL multiple hits with no ancestsos "<<MMPR3(sz,vin[0],vin[sz-1])) } 
return best;
} // i 

return best; 
}

// Eventually m_lut should be int to int and use
// same string tokenizer as this. 
template <class Tv> IdxTy best_taxon( const Tv & vin, const IdxTy flags)
{
IdxTy best=~0;
St & st= m_st;
std::map<IdxTy,IdxTy> m;
normal_taxes(m,st);
std::vector<IdxTy> vint;
MM_LOOP(ii,vin) { vint.push_back(st(*ii)); } 
normal_taxes( m, m_st);
		//(m.find(op.rank())!=m.end()
const IdxTy rootid=root();
IdxTy node=rootid;
IdxTy vinpos=0;
IdxTy vinsz=vin.size();
best=subtree_search(vint,vinsz,vinpos,node,st);
if (best==IdxTy(~0) ) best=node;
MM_ERR(" best "<<MMPR3(best,m_obj[best].id(),  m_lut[m_obj[best].id()]))
return best; 
}

IdxTy subtree_search
	(std::vector<IdxTy> & vint,  const IdxTy & vinsz, IdxTy  vinpos,  IdxTy  node, St & st)
{// FIXME that should be fast lol  make the m_lut int to int 
// and fix the reverse mode 
//const MyStruct & o=m_obj[node]; 
//IdxTy idt=st(m_lut[m_obj[node].id()]);
//IdxTy idt=st(m_lut[o.id()]);

const MyStruct & o=m_obj[node]; 
IdxTy idt=st(m_lut[o.id()]);
if (idt==vint[vinpos])
{ // match, go down in both or finish
++vinpos;
MM_ERR(" match "<<MMPR3(node,idt, st(idt)))
if (vinpos==vinsz) {return node; } // ncbi taxon NOT st token 
// check children
MM_LOOP(ii,o.m_kids)
{
IdxTy tnode=(*ii);
IdxTy next=subtree_search(vint,vinsz,vinpos,tnode,st);
if (next!=IdxTy(~0)) return next;
} // ii 
MM_ERR(" match no kids through "<<MMPR3(node,idt, st(idt)))
return node; 
} // if
else
{ // look horizontally or vertically 


}

return ~0;
} // subtree_search



void traverse_full_tree(const StrTy & p, const IdxTy & flags )
{
// there should be at least one
auto ii=m_obj.begin();
if (ii==m_obj.end()){ MM_ERR(" tax tree empty")  return ; } 
IdxTy rkey=(*ii).second.id();
IdxTy rootid=root(rkey);
MM_ERR( " found root "<<MMPR4((*ii).first,rkey,rootid,m_lut(rootid)))
std::vector<IdxTy> vv;
traverse_tree(rootid,vv);
}

void traverse_tree(const IdxTy& rid, std::vector<IdxTy>  v)
{
auto ii=m_obj.find(rid);
if (ii==m_obj.end()){  return ; } 
const auto & o=(*ii).second;
v.push_back(o.id()); // this SHOULD equal rid... 
MM_LOOP(jj,o.m_kids) { traverse_tree(*jj,v); }

Ss ss;
ss<<MMPR(v.size());
MM_LOOP(jj,v)
{
ss<<":"<<m_lut((*jj));
}
MM_ERR(ss.str()) 

}
typedef ::tree_layout tree_layout;
#if 0
class tree_layout
{
typedef collections_traits::Tr Tr;  
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
class pos { 
public: pos(const IdxTy _id, const IdxTy _x, const IdxTy _y)
{ id=_id; x=_x; y=_y; }
 pos() {}
IdxTy id,x,y; 
};
public:
typedef pos position;
void operator()(const std::vector<IdxTy> & v)
{
const IdxTy sz=v.size();
const IdxTy cnt=m_level_sz[sz];
const IdxTy id=v[sz-1];
++m_level_sz[sz];
m_layout[id]=pos(id,sz,cnt);
// this needs to be done when the parent is positioned
// or work solely on id
//if (sz>1) m_connect.push_back(
if (sz>1) m_parent[id]=v[sz-2];
}
std::map<IdxTy, pos> m_layout;
std::map<IdxTy , IdxTy> m_level_sz;
std::map<IdxTy , IdxTy> m_parent;

}; // tree_layout
#endif

template <class Tf> 
void traverse_full_tree(Tf & f, const StrTy & p, const IdxTy & flags )
{
// there should be at least one
auto ii=m_obj.begin();
if (ii==m_obj.end()){ MM_ERR(" tax tree empty")  return ; } 
IdxTy rkey=(*ii).second.id();
IdxTy rootid=root(rkey);
MM_ERR( " found root "<<MMPR4((*ii).first,rkey,rootid,m_lut(rootid)))
std::vector<IdxTy> vv;
traverse_tree(f,rootid,vv);
}
template <class Tf> 
void traverse_tree(Tf & f, const IdxTy& rid, std::vector<IdxTy>  v)
{
auto ii=m_obj.find(rid);
if (ii==m_obj.end()){  return ; } 
const auto & o=(*ii).second;
v.push_back(o.id()); // this SHOULD equal rid... 
MM_LOOP(jj,o.m_kids) { traverse_tree(f,*jj,v); }
//MM_ERR(" K "<<MMPR(v.size()))
f(v);

/*
Ss ss;
ss<<MMPR(v.size());
MM_LOOP(jj,v)
{
ss<<":"<<m_lut((*jj));
}
MM_ERR(ss.str()) 
*/
}





void dump_lineages(const StrTy & fn, const IdxTy flags)
{
std::map<IdxTy, IdxTy> m;
if ((flags&1)!=0)  { normal_taxes(m,m_st);    }
dump_lineages(fn,m,flags);

}
template <class Tm> 
void dump_lineages(const StrTy & fn, Tm & m, const IdxTy & flags  )
{
MM_ERR("dumping lineages ")
const bool dump_all_nodes= ((flags&2)!=0);
const bool check_names= ((flags&4)!=0);
const IdxTy species_code=m_st(StrTy("species"));
//const IdxTy genus_code=m_st(StrTy("genus"));
if (check_names)
{
for (auto ii=m_lut.fbegin(); ii!=m_lut.fend(); ++ii) 
//	MM_LOOP(ii,m_lut)
	{
		const IdxTy i=(*ii).first;
		auto & o=m_obj[i];
		const StrTy & id=m_lut(o.id()); 
		const IdxTy wc=WordCount(id);
		const IdxTy nac=NonAlphaCount(id);
		if (( wc!=1)||(nac!=0)) 
		{ MM_ERR(MMPR3(wc,nac,id)<<MMPR3(i,m_st(o.rank()),lineage(i,m))) }
	} // ii 
}
else {
//MM_LOOP(ii,m_lut)
for (auto ii=m_lut.fbegin(); ii!=m_lut.fend(); ++ii) 
{
	const IdxTy i=(*ii).first;
	auto & o=m_obj[i];
//	StrTy&  v=(*ii).second;
	if (dump_all_nodes) { MM_MSG(MMPR3(i,m_st(o.rank()),lineage(i,m))) }
	else if (o.rank()==species_code) { MM_MSG(MMPR2(i,lineage(i,m))) }
} // ii 
}

MM_ERR("done dumping lineages ")
} // dump_lineages

// the name lut needs to be fixed. The species contain genus and wanted to make everything
// lower case for now...
// see also conform
//#warning  need to make this do reverse too
/*
 The Forward table maps integer node numbers to a taxon name.
Names are not unique and can be same for different nodes although presumably
shildren of the same parent must be unique. 
So, the reverse table can point to many node numbers. To reconstruct
the forward table from the reverse in the presence of synonyms, the 
prefered name has to be known. /

*/ 
void fix_lut(const bool have_syn=false)
{
const IdxTy species_code=m_st(StrTy("species"));
// obsolete TODO 
//const IdxTy genus_code=m_st(StrTy("genus"));
//MM_LOOP(ii,m_lut)
for(auto ii=m_lut.fbegin(); ii!=m_lut.fend(); ++ii)
{
	const IdxTy i=(*ii).first;
	auto & o=m_obj[i];
	// if we have synonyms, save the old one and move the 
	// reverse entry this is terribly inefficient and moving map
	// entries around would be easier during creation  
	StrTy  v=(*ii).second;
	const StrTy vorg=v;
	if( !have_syn) v=mjm_strings::to_lower(v);
	if (o.rank()==species_code)
	{
		// these APPEAR to be genus names but hard to be sure, no docs... 
		// this may only occur on subspecies, not that far yet lol. 
		if (!have_syn)	v=BraceRemove(v);
		const StrTy & genus=m_lut[m_obj[o.parent()].id()];
		v=NonRedundant(v,genus);	
		if (!have_syn) (*ii).second=v;
		else m_lut.move(v,ii);
	}
if (false) {	if (have_syn) if (v!=vorg)
	{
		// well that is great doh... 
		// find the original vorg entry in m_rev and move it to 
		// v. This may sound easy but memory mess inmap???
		auto rr=m_lut.find(vorg);
		// the reverse table could point to multiple
		// fwd's and already be moved but hopefully
		// all the same way. 
		if (rr==m_lut.rend())
		{
	//		MM_ERR(" DMEL rev missing "<<MMPR2(v,vorg))
		}
//		m_lut(	
	else 	m_lut.rmove(v,vorg);
//		m_lut(v)=(*rr).second; m_lut.m_rev.erase(rr); 
	}}
// remove genus from species 
} // ii 

if (!have_syn) return; 

//for(auto ii=m_lut.rbegin(); ii!=m_lut.rend(); ++ii)
{

//#warning r do the reverse table 

} 

} // fix_lut


IdxTy child( const IdxTy n,const IdxTy i) 
{
auto & node=  m_obj[n];
return node.m_kids[i];
}
IdxTy parent( const IdxTy n) 
{
auto & node=  m_obj[n];
return node.parent();
}

class name_order {
public:
template<class Tx,class Ty> bool operator()(const Tx & a, const Ty & b )
{
return (fudd->m_lut(a))<(fudd->m_lut(b));
}
Myt * fudd;
}; // name_order

void sort_for_ui()
{
name_order no;
no.fudd=this;
MM_LOOP(ii,m_obj)
{

MyStruct & o=(*ii).second; // m_obj[node]; 
std::sort(o.m_kids.begin(), o.m_kids.end( ), no);
//std::sort(o.m_kids.begin(), o.m_kids.end( ));
}

}

typedef mjm_index_lut::Citor  node_iterator;
typedef mjm_index_lut::CRitor  synonym_iterator;

typedef To::iterator  obj_iterator;

obj_iterator begin_obj() { return m_obj.begin(); } 
obj_iterator end_obj() { return m_obj.end(); } 

node_iterator begin_node() { return m_lut.fbegin(); } 
node_iterator end_node() { return m_lut.fend(); } 
synonym_iterator begin_synonym() { return m_lut.rbegin(); } 
synonym_iterator end_synonym() { return m_lut.rend(); } 
synonym_iterator  find_synonym(const StrTy & s) { return m_lut.find(s); } 
node_iterator  find_node(const IdxTy & n) { return m_lut.find(n); } 

mjm_index_lut split_words() { return m_lut.split_words(); } 

 std::vector<IdxTy>   node_numbers( const StrTy & x  ) { return m_lut(x); }
const StrTy & node_name( const IdxTy n ) const { return m_lut(n); }
//const IdxTy&  parent( const IdxTy n) { return m_obj[n].parent(); }  
const MyStruct & node(const IdxTy & n)const  { return (*(m_obj.find(n))).second; } 
const IdxTy  nkids ( const IdxTy n)  { return m_obj[n].m_kids.size(); } 
IdxTy child ( const IdxTy n,const IdxTy i )const  { return node(n).m_kids[i]; } 
std::vector<StrTy> kids ( const IdxTy n) 
{
std::vector<StrTy> kids ;
auto & node=  m_obj[n];
MM_LOOP(ii,node.m_kids)
{
Ss ss;
ss<<m_lut(*ii)<<"_"<<(*ii);
kids.push_back(ss.str()); 
}

return kids;

}
// terminal nodes only 
class tree_terminal_iterator{
typedef std::vector<IdxTy> Ptr;
public:
tree_terminal_iterator(Myt & t, const IdxTy r) : tt(t),m_r(r) ,m_bumps(0)
{
//m_node.push_back(m_r); m_desc.push_back(0);
push(m_r,0);
pop();
if (!valid()) return; 
descend_to_terminal();
/*
while (tt.nkids(node())!=0) 
{
//m_node.push_back(tt.child(node(),0)); m_desc.push_back(0);
push(tt.child(node(),0),0);
}
*/

}
void push(const IdxTy & no, const IdxTy & desc)
{
m_node.push_back(no);
m_desc.push_back(desc);
}

bool  pop()
{
if ( m_desc.back()<tt.nkids(m_node.back())) return false; 
//if ( m_desc.back()>=tt.nkids(m_node.back()))
{
m_desc.pop_back();
m_node.pop_back();
return true;
}
//return false; 
}

bool valid() const { return m_node.size()!=0; } 
void inc()
{
m_bumps=0;
do { 
if (!valid()) break; 
++m_desc.back(); ++m_bumps;
} while (pop()) ;
if (!valid()) return; 
descend_to_terminal();
/*
while (tt.nkids(node())!=0) 
{
//m_node.push_back(tt.child(node(),m_desc.back())); m_desc.push_back(0);
push(tt.child(node(),m_desc.back()),0);
}
*/
//m_node.back()=
}
IdxTy node() { return m_node[m_node.size()-1];}
IdxTy depth() { return m_node.size();}
IdxTy bumps() { return m_bumps;}


void descend_to_terminal()
{
while (tt.nkids(node())!=0) { push(tt.child(node(),m_desc.back()),0); }
}
// this precludes visiting the root... 
//IdxTy node() { return tt.child(m_node.back(),m_desc.back()); } // m_node[m_node.size()-1];}
Ptr  m_desc,m_node;
Myt & tt;
const IdxTy m_r;
IdxTy m_bumps;
}; // tree_terminal_iterator

class tree_hierarch_iterator{
typedef std::vector<IdxTy> Ptr;
public:
tree_hierarch_iterator(Myt & t, const IdxTy r) : tt(t),m_r(r) , m_bumps(0)
{
//m_node.push_back(m_r); m_desc.push_back(0);
// visit the parent first to print the hieratchy.. 
push(m_r,~0);
//pop(); if (!valid()) return; 
//while (tt.nkids(node())!=0) { push(tt.child(node(),0),0); }
//descend_to_terminal();


}
void push(const IdxTy & no, const IdxTy & desc) { m_node.push_back(no); m_desc.push_back(desc); }

bool  pop()
{
if ( m_desc.back()<tt.nkids(node())) return false; 
m_desc.pop_back(); m_node.pop_back();
return true;
}

bool valid() const { return m_node.size()!=0; } 
bool terminal() const { return 0==tt.nkids(node()); }  
void inc()
{
m_bumps=0;
// when all kids even zero are done pop otherwise descend 
do {
if (!valid()) break;
// just push the zero child like we would have done before .. 
//if (terminal())
//{ pop();
// only occurs when root is only node doh 
//if (!valid()) return; //   break;  
//} else{
//++m_desc.back();  ++m_bumps;
//if (pop())  break;
//}
++m_desc.back(); ++m_bumps;
} while (pop());
if (!valid()) return; 
push_child();
//descend_to_terminal();
//m_node.back()=
}

void descend_to_terminal()
{
while (tt.nkids(node())!=0) { push(tt.child(node(),m_desc.back()),0); }
while (tt.nkids(node())!=0) { push_child(); }
}

void push_child() { push(tt.child(node(),m_desc.back()),~0); }


IdxTy depth()const  { return m_node.size();}
IdxTy bumps()const  { return m_bumps;}

IdxTy node()const { return m_node[m_node.size()-1];}
StrTy dump() const
{
Ss ss;
ss<<MMPR4(m_r,m_bumps,m_desc.size(),m_node.size())<<MMPR(valid());
if (valid()) {ss<<MMPR3(tt.node_name(node()),tt.nkids(node()),m_desc.back());}  
return ss.str();
}
// this precludes visiting the root... 
//IdxTy node() { return tt.child(m_node.back(),m_desc.back()); } // m_node[m_node.size()-1];}
Ptr  m_desc,m_node;
Myt & tt;
const IdxTy m_r;
IdxTy m_bumps;
}; // tree_hierarch_iterator









StrTy l2s(const std::vector<StrTy> & lineage)
{
Ss ss; 
StrTy sep="";
MM_LOOP(ii,lineage) {   ss<<sep<<(*ii); sep="_"; } 
return ss.str(); 
}
const St & st() const { return m_st; } 
void clear() { m_lut.clear(); m_st.clear(); m_obj.clear(); } 

private:

StrTy Trim(const StrTy & s)
{
Ss ss;
ss<<s;
StrTy x;
while ( true) { StrTy y; ss>>y; if ( y=="") break; if (x!="") x=x+StrTy(" "); x=x+y; } 
return x; 
}
// BraceRemove 
StrTy BraceRemove(const StrTy & s)
{
IdxTy cnt=0;
const char * p=s.c_str();
const IdxTy sz=s.length();
char c[sz+1];
IdxTy level=0;
for (IdxTy i=0; i<sz; ++i) { const char & b=p[i];
if ( b=='[') { ++level; continue;}
if ( b==']') { --level; continue;}
if ( level==0) { c[cnt]=b; ++cnt; } 
}
c[cnt]=0;
return StrTy(c); 
}

StrTy NonRedundant(const StrTy & s, const StrTy & genus )
{
Ss ss;
ss<<s;
StrTy x;
while ( true) { StrTy y; ss>>y; if (y==genus) continue; if ( y=="") break; if (x!="") x=x+StrTy(" "); x=x+y; } 
return x; 
}

IdxTy WordCount(const StrTy & s)
{
Ss ss;
ss<<s;
IdxTy cnt=0;
while ( true) { StrTy y; ss>>y; if ( y=="") break; ++cnt; } 
return cnt; 
}

IdxTy NonAlphaCount(const StrTy & s)
{
IdxTy cnt=0;
const char * p=s.c_str();
const IdxTy sz=s.length();
for (IdxTy i=0; i<sz; ++i) { const char & c=p[i];
if ((c>='a')&&(c<='z')) continue;
if ((c>='A')&&(c<='Z')) continue;
++cnt;
}
return cnt; 
}




St m_st;
//std::map<IdxTy, StrTy> m_lut;
//collate_order_map  m_lut;
//St m_lut;
mjm_index_lut  m_lut;


}; //mjm_tax_tree


// #undef MM_DMEL

#endif

