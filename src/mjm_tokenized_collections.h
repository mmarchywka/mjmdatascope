#ifndef MJM_TOKENIZED_COLLECTIONS_H__
#define MJM_TOKENIZED_COLLECTIONS_H__

#include "mjm_globals.h"
#include "mjm_string_tokenizer.h"
//#include "mjm_collections.h"
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

#define LUT_LOAD_FWD  0

// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#ifndef MM_DMEL
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
#endif

#define MM_DMELF(fm,code,x) f.DMel(code, f.m_ss<<MM_STR_LOC<<" "<<x); 

namespace tokenized_collections_traits
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
typedef string_tokenizer St;
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // tokenized_collection_traits ns 






#if 0 


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

#endif

#if 1 
// turning into an R data frame with ad hoc crap lol 
class mjm_tokenized_ragged_table
{
typedef tokenized_collections_traits::Tr Tr;
typedef mjm_tokenized_ragged_table Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::St St;
typedef Tr::MyBlock  MyBlock;
public:
//typedef std::vector<StrTy> Line;
typedef IdxTy WordTy;
typedef std::vector<WordTy> Line;
typedef std::vector<Line> Lines;
typedef std::map<StrTy, StrTy>  xlate_map;
typedef std::map<StrTy, D>  xlate_double_map;
typedef std::vector<IdxTy> LocVec;
typedef std::map<StrTy, LocVec> line_index_map;
//typedef IdxTy WordTy;
//typedef std::vector<WordTy> Key;
typedef std::vector<StrTy> ExtLine;

typedef std::vector<IdxTy> sort_order;
mjm_tokenized_ragged_table()
: m_min_words(1),m_max_words(0),m_sep(" "),m_ignore_hash(false) {
NewSt(0);
}

mjm_tokenized_ragged_table(const IdxTy lines, const IdxTy words)
: m_min_words(1),m_max_words(0),m_sep(" "),m_ignore_hash(false) {
NewSt(0);
 m_lines=Lines(lines);
const Line v= Line(words);
MM_LOOP(ii,m_lines) (*ii)=v;
}

~mjm_tokenized_ragged_table()
{ if (m_delete_st) delete m_st; } 

St & st() { return *m_st; }
ExtLine st(const Line & k ) const  { return (*m_st)(k); }
Line st(const ExtLine & k )   { return (*m_st)(k); }
StrTy st(const IdxTy & k ) const  { return (*m_st)(k); }
IdxTy st(const StrTy & k ) { return (*m_st)(k); }


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



void sep(const StrTy & s) { m_sep=s;} 
const Line & line(const IdxTy i ) const { return m_lines[i];}
Line & ncline(const IdxTy i )  { return m_lines[i];}
void clear()  {  m_lines.clear();}
const IdxTy  size() const { return m_lines.size(); } 
const void ignore_hash(const bool x )  { m_ignore_hash=x; } 
void add(const Line & li ) { m_lines.push_back(li); }  
void add(const ExtLine & li ) { m_lines.push_back(st(li)); }  

void load_from_string(const StrTy & s, const bool debug=false)
{
Ss ss;
ss<<s;
CommandInterpretter li(&ss);
ConmfigureReader( li);
li.use_alt_eol(',',true);
Load(li,debug);
}

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
	// TODO FIXME ASSFUDD putting ok() in nextok() IS FUDDED ASDFASDFASD
		if (!li.ok()) break;
        const IdxTy sz=li.size();
		if ( m_ignore_hash) if (sz>0) if (li.words()[0].c_str()[0]=='#') continue;
		if (m_min_words>m_max_words) m_min_words=sz;
		Add(li.words());
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
ConmfigureReader(  li);
MM_ONCE(" changed fudding load_lines to use configure but may be fuxekd",)
//    li.set_split(1,' ');
    while (li.nextok())
    {
	// TODO FIXME ASSFUDD putting ok() in nextok() IS FUDDED ASDFASDFASD
		if (!li.ok()) break;
        const IdxTy sz=li.size();
		if (m_min_words>m_max_words) m_min_words=sz;
		typedef std::vector<StrTy> W;
		W w;
		w.push_back(li.line());
		Add(w);
		//if (debug) { MM_ERR(MMPR3(li.line(),sz,(*this).size())) }
		if (debug) { MM_ERR(MMPR4(li.line(),sz,(*this).size(),li.ok())) }
    } // nextok()

} // load_lines






void add(const StrTy & line)
{

CommandInterpretter li;
if (false) { MM_ERR(" adding "<<MMPR(line)) } 
li.nextline(line.c_str());
Add(li.words());

}
void write_file(const StrTy & fn, const IdxTy  flags=0)  const
{
std::ofstream ofs(fn);
ofs<<dump(flags);

}
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
const StrTy & nm=st(line[nidx]);
const StrTy & v=st(line[vidx]);
//MM_ERR(MMPR2(nm,v))

rwm.set( nm,v,true);
}

}

StrTy dump( const IdxTy flags=0) const
{
Ss ss;
IdxTy i=0;
const bool use_space=((flags&1)!=0);
const bool add_seq=((flags&2)==0);

StrTy seps="|";
if (use_space) seps=" ";

const StrTy sep=seps;
MM_LOOP(ii,m_lines)
{
IdxTy field=0;
const Line & line=(*ii);
bool add_sep=false;
if (add_seq) {add_sep=true;  ss<<i; }
MM_LOOP(jj,line)
{
// TODO FIXME  this should quote or escape where needed
if ( add_sep) ss<<sep;
ss<<st(*jj);
++field;
add_sep=true;
}
ss<<CRLF;
++i;
} 

return ss.str();
}

line_index_map line_index(const IdxTy key)
{
line_index_map m;
MM_SZ_LOOP(i,m_lines,sz)
{
const ExtLine  line=st(m_lines[i]);
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
const ExtLine  line=st(m_lines[i]);
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
const ExtLine  line=st(m_lines[i]);
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
const ExtLine  line=st(m_lines[i]);
const IdxTy sz=line.size();
if ( sz>mf) { 
const D x=atof(line[val].c_str());
auto ii=m.find(line[key]);
//MM_ERR(" FUDD "<<MMPR4(key,val,x,line[key]))
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
const ExtLine  line=st(m_lines[i]);
const IdxTy szll=line.size();
const StrTy canon=(key<szll)?line[key]:def;
// and mapping canon to canon is desired. 
MM_SZ_LOOP(j,line,szl)
{
 m[(line[j])]=(canon); 
} // j 
} // i 
return m;
}

void append(const IdxTy & line, const StrTy & v)
{
m_lines[line].push_back(st(v));
}

Myt & operator+=(const Myt & that )
{
MM_ERR(" danger will robinson needs to share the m_st lol")
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
MM_ERR(" ASSFUDD SHUIT FUDD  "<<MMPR2(i,sz))
Line x(sz);
for(IdxTy j=0; j<sz; ++j){
MM_ERR(" ASDCASDCASCASDCASDCASD  "<<MMPR2(j,m_lines[j].size()))

 x[j]=m_lines[j][i]; 

}
d.m_lines[i]=x;
}


}




private:
void NewSt(St * p) 
{
if (p==0) { m_st= new St(); m_delete_st=true; return; } 
m_st=p;
m_delete_st=false; 

}

void ConmfigureReader( CommandInterpretter & li,const bool debug=false)
{
li.readline_ok(false);
if (m_sep.length()>0)
    li.set_split(1,m_sep.c_str()[0]);
    else li.set_split(1,' ');
MM_ERR(" ASSFUDD "<<MMPR2(m_sep,m_sep.length()))
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

void Add( const ExtLine & x)
{
	const IdxTy sz=x.size();
	if (sz<m_min_words) m_min_words=sz;
	if (sz>m_max_words) m_max_words=sz;
	m_lines.push_back(st(x)); 
}

void Add( const Line & x)
{
	const IdxTy sz=x.size();
	if (sz<m_min_words) m_min_words=sz;
	if (sz>m_max_words) m_max_words=sz;
	m_lines.push_back((x)); 
}




St * m_st;
bool m_delete_st;
Lines m_lines;
IdxTy m_min_words, m_max_words;
StrTy m_sep;
bool m_ignore_hash;
}; // mjm_tokenized_ragged_table

#endif

// turning into an R data frame with ad hoc crap lol 
// this is a bad impl if the keys are highly reudnudant but 
// the tokenizer helps and there ar ways to fix it but the
// variable key lenght may be useful. 

class mjm_tokenized_vector_map
{
typedef tokenized_collections_traits::Tr Tr;
typedef mjm_tokenized_vector_map Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::St St;
typedef Tr::MyBlock  MyBlock;
public:
//typedef std::vector<StrTy> Line;
typedef IdxTy WordTy;
typedef std::vector<WordTy> Key;
typedef std::vector<StrTy> ExtKey;
typedef D Value;
typedef std::map<Key,Value> IdxMap;

mjm_tokenized_vector_map() { NewSt(0); }

~mjm_tokenized_vector_map()
{ if (m_delete_st) delete m_st; } 

St & st() { return *m_st; } 
ExtKey st(const Key & k ) const  { return (*m_st)(k); } 
Key st(const ExtKey & k ) { return (*m_st)(k); } 
void clear()  { m_map.clear(); }
const IdxTy  size() const { return m_map.size(); } 
void add(const ExtKey & k, const Value & v ) { m_map[st(k)]+=v; }  
void add(const Key & k, const Value & v ) { m_map[(k)]+=v; }  


void dump(OsTy & os=std::cout,  const IdxTy flags=0) const
{
//IdxTy i=0;
const StrTy sep=" ";
MM_LOOP(ii,m_map)
{
const ExtKey k=st((*ii).first);
const Value & v=(*ii).second;
const IdxTy sz=k.size();
if (sz!=0) os<<k[0];
for(IdxTy i=1; i<sz; ++i) os<<sep<<k[i];
os<<sep<<v<<CRLF;
//++i;
} 

}


StrTy to_string( const IdxTy flags=0) const
{
Ss ss;
IdxTy i=0;
MM_LOOP(ii,m_map)
{
++i;
} 
return ss.str();
}

private:
void NewSt( St *  p) 
{
if (p==0) { m_st= new St(); m_delete_st=true; return; } 
m_st=p;
m_delete_st=false; 

}


St * m_st;
bool m_delete_st;
IdxMap m_map;
}; // mjm_tokenized_ragged_table



class mjm_group_sample_map
{
typedef tokenized_collections_traits::Tr Tr;
typedef mjm_group_sample_map Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::St St;
typedef Tr::MyBlock  MyBlock;
public:
//typedef std::vector<StrTy> Line;
typedef IdxTy Otu;
typedef IdxTy AnGroup;
typedef IdxTy Sample;
typedef IdxTy Group;
typedef StrTy ExtTy;
//typedef std::vector<WordTy> Key;
//typedef std::vector<StrTy> ExtKey;
typedef D Value;
typedef std::map<Otu,Value> IdxMap;
typedef std::map<Sample,IdxMap> SampleMap;
typedef std::map<Group,SampleMap> GroupMap;
typedef std::map<IdxTy,IdxTy> ListMap;

mjm_group_sample_map() { NewSt(0);Init(); }
mjm_group_sample_map(St * p ) { NewSt(p);Init(); }

~mjm_group_sample_map()
{ if (m_delete_st) delete m_st; } 

St & st() { return *m_st; } 
ExtTy st(const Group & k ) const  { return (*m_st)(k); } 
Group st(const ExtTy & k ) { return (*m_st)(k); } 
void clear()  { m_map.clear(); }
const IdxTy  size() const { return m_map.size(); } 
//void add(const ExtKey & k, const Value & v ) { m_map[st(k)]+=v; }  
//void add(const Key & k, const Value & v ) { m_map[(k)]+=v; }  

void add(const Group & g, const Sample & s, const Otu & o ,  const Value & v ) 
{ m_map[g][s][o]+=v; 
 m_map[g][s][m_otu_total]+=v; 
++m_groups[g];
++m_samples[s];
++m_otus[o];
}  
void add(const Group & g, const Sample & s, const AnGroup & a, const Otu & o ,  const Value & v ) 
{ m_map[g][s][o]+=v; 
 m_map[g][s][m_otu_total]+=v; 
++m_groups[g];
++m_samples[s];
++m_otus[o];
++m_angroups[a];
auto ii=m_sample_to_angroup.find(s);
if (ii==m_sample_to_angroup.end()) m_sample_to_angroup[s]=a;
else
{
const IdxTy oldan=(*ii).second;
if (oldan!=a)
{
MM_ERR(" inconsistent groups "<<MMPR4(st(g),st(s),st(a),st(o))<<MMPR(st(oldan)))
}

}
}  

void add(const ExtTy & g, const ExtTy & s, const ExtTy & o ,  const Value & v ) 
{ 
	// m_map[st(g)][st(s)][st(o)]+=v; 
	add(st(g),st(s),st(o),v); 

}  

void add(const ExtTy & g, const ExtTy & s, const ExtTy & an, const ExtTy & o ,  const Value & v ) 
{ add(st(g),st(s),st(an),st(o),v); }  

void dump(OsTy & os=std::cout,  const IdxTy flags=0) const
{
//IdxTy i=0;
const StrTy sep=" ";
MM_LOOP(ii,m_map)
{
//const ExtKey k=st((*ii).first);
//const Value & v=(*ii).second;
//const IdxTy sz=k.size();
//if (sz!=0) os<<k[0];
//for(IdxTy i=1; i<sz; ++i) os<<sep<<k[i];
//os<<sep<<v<<CRLF;
//++i;
} 

}
void calculate()
{
MM_LOOP(gg,m_groups)
{
const IdxTy gi=(*gg).first;
const StrTy group=st(gi);
const auto gmi=m_map.find(gi);
if (gmi==m_map.end()) continue;
const auto gm=(*gmi).second;
MM_LOOP(ss,m_samples)
{
const IdxTy si=(*ss).first;
const StrTy sample=st(si);
const auto smi=gm.find(si);
if (smi==gm.end()) continue;
auto sm=(*smi).second;
IdxTy nzed=0;
D s=0;
D s2=0;
D sl=0;
D sl2=0;
D vmin=1;
D vmax=0;
D h=0;
IdxTy n= m_otus.size()-m_ndummy;
D norm=1.0/sm[m_otu_total];
MM_LOOP(oo,m_otus)
{
const IdxTy oi=(*oo).first;
const StrTy otu=st(oi);
const auto omi=sm.find(oi);
if (omi==sm.end()) {++nzed;   continue; } 
const D v=(*omi).second*norm;
if (v==0) { ++nzed; continue; }
if (v<vmin) vmin=v;
if (v>vmax) vmax=v;
s=s+v;
s2=s2+v*v;
const D vl=log(v);
sl=sl+vl;
sl2=sl2+vl*vl;
h=h-v*vl;
} // oo
MM_MSG(MMPR4(group,sample,nzed,n)<<MMPR4(vmin,vmax,s,s2)<<MMPR4(sl,sl2,exp(sl/n),h))

} // ss
} // gg


} // calculate

void calculate_otu()
{
MM_LOOP(gg,m_groups)
{
const IdxTy gi=(*gg).first;
const StrTy group=st(gi);
const auto gmi=m_map.find(gi);
if (gmi==m_map.end()) continue;
const auto gm=(*gmi).second;
MM_LOOP(oo,m_otus)
{
const IdxTy oi=(*oo).first;
const StrTy otu=st(oi);

IdxTy nzed=0;
D s=0;
D s2=0;
D sl=0;
D sl2=0;
D vmin=1;
D vmax=0;
D h=0;
IdxTy n=0; //  m_otus.size()-m_ndummy;
std::map<D,IdxTy> p;
MM_LOOP(ss,m_samples)
{
const IdxTy si=(*ss).first;
//const IdxTy ansi=m_sample_to_angroup[si];
//const StrTy angroup=st(ansi);
//const StrTy sample=st(si);
const auto smi=gm.find(si);
// TODO could iterate over the group submap? 
if (smi==gm.end()) continue;
++n;
auto sm=(*smi).second;
const auto omi=sm.find(oi);
if (omi==sm.end()) {++nzed;   continue; } 
D norm=1.0/sm[m_otu_total];

const D v=(*omi).second*norm;
if (v==0) { ++nzed; continue; }
if (v<vmin) vmin=v;
if (v>vmax) vmax=v;
s=s+v;
s2=s2+v*v;
const D vl=log(v);
++p[vl];
sl=sl+vl;
sl2=sl2+vl*vl;
h=h-v*vl;
} // ss
D fnzed=1;
if (n!=0){ fnzed=1.0*nzed/n;  s=s/n; sl=exp(sl/n); } 
IdxTy c=nzed;
auto ii= p.begin();
while (ii!=p.end()){if (c>=(n>>1)) break; c+=(*ii).second;  ++ii; }
h=(*ii).first;
MM_MSG(MMPR4(group,otu,n,fnzed)<<MMPR4(vmin,vmax,s,s2)<<MMPR4(sl,sl2,sl,h))

} // oo
} // gg


} // calculate_otu

// needs a header ( otu followed by group names ) 
template <class Ragged> // FUDD 
void calculate_otu(Ragged & rout)
{
typedef typename Ragged::Line Fock;
Fock line;
if (false){
line.push_back("otu");
MM_LOOP(gg,m_groups)
{
line.push_back(st((*gg).first)); 
}}
else
{
line.push_back("hierarchy");
MM_LOOP(oo,m_otus)
{
line.push_back(st((*oo).first)); 
}}

rout.add(line);

MM_LOOP(gg,m_groups)
{
const IdxTy gi=(*gg).first;
const StrTy group=st(gi);
const auto gmi=m_map.find(gi);
if (gmi==m_map.end()) continue;
const auto gm=(*gmi).second;
MM_LOOP(ss,m_samples)
{
const IdxTy si=(*ss).first;
const IdxTy ansi=m_sample_to_angroup[si];
const StrTy angroup=st(ansi);
const StrTy sample=st(si);
const auto smi=gm.find(si);
if (smi==gm.end()) continue;
auto sm=(*smi).second;
Fock line;
line.push_back(angroup);
D total=sm[m_otu_total];
D norm=1.0/sm[m_otu_total];
if (total==0)
{
MM_ERR(" total is zero, "<<MMPR3(group,sample,angroup))
norm=1;

}
MM_LOOP(oo,m_otus)
{
const IdxTy oi=(*oo).first;
const StrTy otu=st(oi);
D v=0;
const auto omi=sm.find(oi);
if (omi!=sm.end()) { v=(*omi).second*norm; }
//if (v==0) { ++nzed; continue; }
Ss ss;
ss<<std::setprecision(10)<<v;
line.push_back(ss.str());
} // oo
//MM_MSG(MMPR4(group,sample,nzed,n)<<MMPR4(vmin,vmax,s,s2)<<MMPR4(sl,sl2,exp(sl/n),h))
rout.add(line);
} // ss
} // gg




#if 0 
MM_LOOP(gg,m_groups)
{
const IdxTy gi=(*gg).first;
const StrTy group=st(gi);
const auto gmi=m_map.find(gi);
if (gmi==m_map.end()) continue;
const auto gm=(*gmi).second;
MM_LOOP(oo,m_otus)
{
const IdxTy oi=(*oo).first;
const StrTy otu=st(oi);

IdxTy nzed=0;
D s=0;
D s2=0;
D sl=0;
D sl2=0;
D vmin=1;
D vmax=0;
D h=0;
IdxTy n=0; //  m_otus.size()-m_ndummy;
std::map<D,IdxTy> p;
MM_LOOP(ss,m_samples)
{
const IdxTy si=(*ss).first;
//const StrTy sample=st(si);
const auto smi=gm.find(si);
// TODO could iterate over the group submap? 
if (smi==gm.end()) continue;
++n;
auto sm=(*smi).second;
const auto omi=sm.find(oi);
if (omi==sm.end()) {++nzed;   continue; } 
D norm=1.0/sm[m_otu_total];

const D v=(*omi).second*norm;
if (v==0) { ++nzed; continue; }
if (v<vmin) vmin=v;
if (v>vmax) vmax=v;
s=s+v;
s2=s2+v*v;
const D vl=log(v);
++p[vl];
sl=sl+vl;
sl2=sl2+vl*vl;
h=h-v*vl;
} // ss
D fnzed=1;
if (n!=0){ fnzed=1.0*nzed/n;  s=s/n; sl=exp(sl/n); } 
IdxTy c=nzed;
auto ii= p.begin();
while (ii!=p.end()){if (c>=(n>>1)) break; c+=(*ii).second;  ++ii; }
h=(*ii).first;
MM_MSG(MMPR4(group,otu,n,fnzed)<<MMPR4(vmin,vmax,s,s2)<<MMPR4(sl,sl2,sl,h))

} // oo
} // gg
#endif


} // calculate_otu









StrTy to_string( const IdxTy flags=0) const
{
Ss ss;
IdxTy i=0;
MM_LOOP(ii,m_map)
{
++i;
} 
return ss.str();
}

private:
void NewSt( St *  p) 
{
if (p==0) { m_st= new St(); m_delete_st=true; return; } 
m_st=p;
m_delete_st=false; 

}
void Init()
{
m_otu_total_name="#TOTAL";
m_otu_total=st(m_otu_total_name);
m_ndummy=1;
}

St * m_st;
bool m_delete_st;
StrTy m_otu_total_name;
IdxTy m_otu_total;
IdxTy m_ndummy;
GroupMap m_map;
ListMap m_groups,m_samples, m_sample_to_angroup,m_angroups, m_otus;
}; // mjm_tokenized_ragged_table



#undef MM_DMEL

#endif

