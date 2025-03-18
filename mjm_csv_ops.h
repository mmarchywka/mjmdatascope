#ifndef MJM_CSV_OPS_H__
#define MJM_CSV_OPS_H__

#include "mjm_text_data.h"
#include <stdlib.h>
// #include <tcl.h>
//#include "mjm_tcl_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

namespace mjm_csv 
{
typedef unsigned int IdxTy;
typedef  int IntTy;
typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::stringstream SsTy;


typedef std::vector<StrTy > LineBaseTy;

class mjm_csv_line : public LineBaseTy
{
typedef mjm_csv_line Myt;
typedef LineBaseTy Super;
typedef Super::value_type FieldTy;

public:
// this will mess up m_ok unless push_back and clear are derived.
mjm_csv_line(): m_is_comment(false),m_ok(false) {}
bool is_comment() const { if ( !m_ok) Check(); return m_is_comment; } 
bool is_blank() const { return size()==0; } 
// can keep a copy of original csv LINE literal when space is
// not an issue. Needs to be able to say "is_comment" etc.

private:

void Check() const
{
if ( size()==0) m_is_comment=false;
else 
{
const FieldTy  & x=(*this)[0];
// zero length means this is zero... 
m_is_comment=(x.c_str()[0]=='#');

} // else
m_ok=true;
} //Check

mutable bool m_is_comment,m_ok;

};
typedef mjm_csv_line LineTy;



// note this is SIGNED 
typedef std::vector<IntTy > IntLineTy;
typedef std::vector<LineTy> FileBaseTy;


class mjm_csv_file : public FileBaseTy
{
public:
typedef std::vector<StrTy> ListTy;
typedef std::vector<IntTy> IntListTy;
typedef std::map<StrTy, IdxTy> ColIndex;
enum { BAD=~0U};
void add_column(const StrTy & nm) { Add(m_column_names,nm); } 
// comment line as first is header names
void set_columns(const LineTy & names) { SetColumns(m_column_names,names); } 
template <class Ty> void set_column(const Ty & nm,const IdxTy n) 
{ Set(m_column_names,StrTy(nm),n);  } 
template < class Ty> void set_columns(const Ty & c )  { Set(m_column_names,c); } 
StrTy name(const IdxTy i ) const
{
if ( i>=m_column_names.size()) return "";
return m_column_names[i];

}
// column name location index
void make_index() { Index(); }

IdxTy column(const char *  nm) const { return column(StrTy(nm)); } 
IdxTy column(const StrTy & nm) const
{ auto ii= m_col_index.find(nm);
if ( ii!=m_col_index.end()) { return (*ii).second; }
return BAD; } 



StrTy dump_columns(const StrTy & sep ) const 
{ SsTy ss; for (IdxTy i=0; i<m_column_names.size(); ++i) 
{ if ( i!=0) ss<<sep; ss<<m_column_names[i]; }
return ss.str(); }
// some utility until we come up with data types etc

static void atoi( IntListTy & d, const ListTy & s) 
{
for ( auto i=s.begin(); i!=s.end(); ++i)
{ d.push_back(::atoi((*i).c_str())); }  
}



private:
void Add ( ListTy & l, const StrTy & nm) { l.push_back(nm); } 
void Set(ListTy & l , const StrTy & nm,const IdxTy n) 
{ if ( l.size()<=n) l.resize(n+1); l[n]=nm; } 
template < class Ty> void Set(ListTy & l, const Ty & c )  { l=c; } 

void SetColumns(ListTy & l, const LineTy & names) 
{ 
l.clear();
for ( IdxTy i=0; i<names.size(); ++i) l.push_back(names[i]);
if (names.is_comment()) if(l[0].length()>0) l[0]=StrTy(names[0].c_str()+1);
} 

void Index()
{
m_col_index.clear();
for ( IdxTy i=0; i<m_column_names.size(); ++i) 
m_col_index[m_column_names[i]]=i; 

}
ListTy m_column_names;
ListTy m_column_types; // type if available 
// these really nee line numbers etc. 
// do not use for comments in the file... 
ListTy m_comments; //  free form when reading comments put here.  

ColIndex m_col_index;

}; // mjm_csv_file

typedef mjm_csv_file FileTy;
// this should be line to line realy but key can be compoisted to a string for now 
typedef std::map<StrTy,LineTy> IndexMapTy;
typedef std::map<StrTy,LineTy> BinarySiteMapTy;
typedef std::map<LineTy,LineTy> LineIndexMapTy;
typedef std::map<LineTy,FileTy> GeneralLineIndexMapTy;
typedef std::map<LineTy,FileTy> InverseLineIndexMapTy;
typedef std::map<StrTy,StrTy> StrIndexMapTy;


class mjm_csv_format
{

public:
enum { USE_STRING_RULES=0 // whatever the >> operator does 
,QUOTES=1 // quotes around each field
,QUOTE_WHEN_SEP=2 // put quote aruond field when it contains a sep
,ESCAPE_SEP=4 // put backslash before included sep char
,ESCAPE_QUOTES=8 // put backslash before included quote, really needed be def 
,OUTPUT_HEADER=256
,COMMENT_HEADER=512 // put comment on header line 


};

mjm_csv_format () : m_branch(USE_STRING_RULES),m_sep(StrTy(" ")){}
IdxTy branch() const { return m_branch; } 
const StrTy & sep() const  { return m_sep; } 
bool write_header() const { return true; }
bool comment_header() const { return true; }

private:

IdxTy m_branch;
StrTy m_sep;



}; // format

typedef mjm_csv_format FormatRules;

// this really contains element stats and a dictionary
// but the two are seaprate. 
// fix this 
class mjm_element_stats
{
class element_info
{
typedef element_info Myt;
public:
element_info(): m_count(0) {}
Myt& operator++() { ++m_count; return *this;}
private:

IdxTy m_count;

}; // info 
typedef element_info ElemInfo; 
typedef StrTy ElemTy;
typedef std::map<ElemTy,ElemInfo> StatMap;
public:
mjm_element_stats() {}
void load(const FileTy & f) { Load(f); } 

private:
void Load(const FileTy & f)
{
for ( auto i=f.begin(); i!=f.end(); ++i)
{
const LineTy & line=(*i);
for ( auto j=line.begin(); j!=line.end(); ++j)
{
// just increment gross count across all fields 
++m_map[(*j)];

} // j 

} //i 


}


StatMap m_map;


}; // mjm_element_stats

typedef  mjm_element_stats ElementStats;

// a dictionary maps strings in csv file into codes
// these can be fixed or adpative and seek to either
// compress, abbreviate human readable, anonymize,
// etc
// really wantg regex here but maybe best done with sed. 
class mjm_csv_dictionary
{
//typedef element_info ElemInfo; 
public:
typedef StrTy ElemTy;
typedef StrTy ShortTy;

void add(const ElemTy & e, const ShortTy & v) {m_map[e]=v;}
void add_decode(const ElemTy & e, const ShortTy & v) {m_decode_map[v]=e;}
// this could fail due to not existing
// maps can either adapt or fail 
const ShortTy & encode(const ElemTy &e ) { return Encode(e); } 
// missing here is an error, could check fwd map lol 
const ElemTy & decode(const ShortTy &v ) { return Decode(v); } 

// arrgh make dict from stats in ex and rules in ds
template <class Tx, class Ty> void make (const Tx & ex, const Ty & ds)
{

}

private:
typedef std::map<ElemTy , ShortTy> DictMap;
typedef std::map<ShortTy,ElemTy > DecodeDictMap;
const ShortTy & Encode(const ElemTy & e ) { return (e); } 
const ElemTy & Decode(const ShortTy &v ) { return v; } 
DictMap m_map;
DecodeDictMap m_decode_map;


}; 
typedef mjm_csv_dictionary Dictionary;



class mjm_csv_dictionary_spec
{
public:

private:

}; 
typedef mjm_csv_dictionary_spec DictionarySpec;




class mjm_csv_ops
{
//
typedef LineTy::const_iterator LI;
typedef FileTy::const_iterator FI;
typedef IndexMapTy::const_iterator II;
typedef LineIndexMapTy::const_iterator LII;
typedef InverseLineIndexMapTy::const_iterator ILII;
public:
template <class Ty> static void dump_generic(OsTy & os,const Ty & x, const StrTy & label="dump_generic")
{
const StrTy sep=" ";
const IdxTy sz=x.size();
IdxTy idx=0;
for ( IdxTy i=0; i<sz; ++i)
{
const auto & line=x[i];
const IdxTy szl=line.size();
os<<label;
for ( IdxTy j=0; j<szl; ++j)
{
os<<sep<<line[j];

} //j 
os<<CRLF;
++idx;
} // i 



}


static void dump_line(OsTy & os, const LineTy & line, const StrTy & sep,const IdxTy leadsep=0)
{
LI li=line.begin(); 
if ( leadsep==1) { if (li!=line.end()) {os<<(*li); ++li; } }
while (li!=line.end()) 
{ os<<sep<<(*li); ++li; }

}
static void dump_file(OsTy & os, const FileTy & file, const StrTy & sep,const IdxTy leadsep=0)
{
FI fi=file.begin();
while ( fi!=file.end())
{
dump_line(os,(*fi),sep,1);
if ( leadsep==0) os<<CRLF;
if ( leadsep==1) os<<sep;
++fi;}



}

static void dump_map(OsTy & os, const InverseLineIndexMapTy & s, const StrTy & sep)
{
ILII ii=s.begin();
while ( ii!=s.end())
{
dump_line(os,(*ii).first,sep,1);
os<<sep<<"|"; 
const FileTy & f=(*ii).second;
dump_file(os,f,sep,1);

os<<CRLF;
++ii;
}


}

static void dump_map(OsTy & os, const IndexMapTy & s, const StrTy & sep)
{
II i=s.begin();
while ( i!=s.end())
{
os<<((*i).first)<<sep<<"|";
const LineTy & line= (*i).second;
dump_line(os,line,sep); //LI li=line.begin(); while (li!=line.end()) { os<<sep<<(*li); ++li; }
os<<CRLF;
++i; 
}

}

static StrTy dump( const LineTy & l)
{
SsTy ss;
const StrTy sep1=" ";
const StrTy sep2=",";
 
for ( IdxTy i=0; i<l.size(); ++i) { if ( i!=0 ) ss<<sep2; ss<<i<<sep1<<l[i]; } 
return ss.str();
}

static void convert(IntLineTy & il, const LineTy & l)
{

LI li=l.begin(); 
while ( li!=l.end())
{
il.push_back(::atoi((*li).c_str())); 
++li; }
}

// include the forward map for easer, 
static void make_binary_site_map(BinarySiteMapTy & bm, const InverseLineIndexMapTy & im,const LineIndexMapTy & idx)
{
//first go through each line-key and compile a list of all sites
typedef std::vector< ILII> IV;
typedef IV::const_iterator IVII;
typedef std::map<IntTy,IV >  LoMap;
typedef LoMap::const_iterator LMI;
// this is for something stupid....
const IdxTy seqs=idx.size();
IdxTy io=0;
IdxTy order[seqs];
for ( IdxTy i=0; i<seqs; ++i) order[i]=~0;  // this shouuld be somethign invald
LoMap lomap;
ILII ilii=im.begin();
//MM_MSG("")
while ( ilii!=im.end())
{
const LineTy & line=(*ilii).first;
LI li = line.begin();
while ( li!=line.end())
{
IntTy val=::atoi((*li).c_str()); 
lomap[val].push_back(ilii); // this relies on map being stable and copy of value.. 
++li; }
// as long as we are going through making the keys again, wth just make
// the ordering based on the already existin lists too LOL
const FileTy & fi=(*ilii).second;
FI fii=fi.begin();
while ( fii!=fi.end())
{
// terribly specialized
const IntTy seq=::atoi((*fii)[0].c_str());
//MM_MSG(" seq is "<<seq<<" and sv is "<<sv);
order[io]=seq;
++io;
++fii; }



 
++ilii; }
std::vector<IntTy> site_values;
LMI lmi=lomap.begin();
while ( lmi!=lomap.end()) { site_values.push_back((*lmi).first); ++lmi; } 
std::sort(site_values.begin(),site_values.end());
MM_MSG("")

const IdxTy sites=lomap.size();
typedef char CodeTy;
CodeTy codes[seqs][sites];
::memset(codes,'0' ,seqs*sites*sizeof(CodeTy));
//LII lii=idx.begin(); while ( lii!=idx.end()) { const LineTy & key =(*lii).first; ++lii; } 
std::vector<IntTy>::const_iterator  svi=site_values.begin();
//MM_MSG("")
IdxTy siteno=0;

while ( svi!=site_values.end())
{
MM_MSG(" site "<<siteno<<" value "<<(*svi))
++svi;++siteno;
}
svi=site_values.begin();
siteno=0;
// this is dumb, the allowed keyw are in teh map originally do
// not need to do exhasutive search again
while ( svi!=site_values.end())
{
//MM_MSG(" processing site "<<(*svi))
const IdxTy sv=(*svi);
IVII ivii=lomap[sv].begin();
IVII ivie=lomap[sv].end();
while (ivii!=ivie)
{ // the first part of this is the key, which contains sv at some point and thenthesecond part
// is a list of sequences that need serial numbers. 
const FileTy & fi = (*(*ivii)).second;
// this is specialized to a single numeric field as the key will fail else...
FI fii=fi.begin();
while ( fii!=fi.end())
{
// terribly specialized
const IntTy seq=::atoi((*fii)[0].c_str());
//MM_MSG(" seq is "<<seq<<" and sv is "<<sv);

codes[seq][siteno]='1';
++fii; }

++ivii; }
//MM_MSG("")

++siteno;
++svi; }

// printing code, the marix needs its own object
StrTy sep=" ";
OsTy & os=std::cout;
if ( (seqs*sites)<1) return; 
CodeTy * plast =codes[order[0]];
std::vector<IdxTy> sames;
for ( IdxTy i=0; i<seqs; ++i)
{
const IdxTy oi=order[i];
if ( oi==~0U ){  MM_MSG(" invalid at "<<i<<" with "<<((int)oi)) continue; } 
CodeTy * p =codes[oi];
IdxTy diffs=0,ones=0;;
for ( IdxTy j=0; j<sites; ++j){ if ( p[j]!='0') ++ones; if ( plast[j]!=p[j]) ++diffs; } 
const bool print_this=(diffs>0)||(i==0);
if ( print_this) {
os<<sep<<sames.size()<<CRLF;
if ( sames.size()!=0)
{
os<<MM_MARK<<" same as ";
for ( IdxTy ss=0; ss<sames.size(); ++ss) os<<sep<<sames[ss];
sames.clear();
os<<CRLF;

 } //sames

os<<i<<sep<<order[i]<<sep<<diffs<<sep<<ones<<"\t";
for ( IdxTy j=0; j<sites; ++j)
{ 

if ( p==plast) os<<p[j]; 
else {
const IdxTy pcode=((p[j]=='0')?0:1)+((plast[j]=='0')?0:2);
// if (plast[j]==p[j]) os<<" "; else os<<p[j];
switch ( pcode) {
case 0: os<<'0'; break;
case 1: os<<'+'; break;
case 2: os<<'-'; break;
case 3: os<<'1'; break;
}; 


 } 
}
// this is no delayed until sames sizeis abaiable but this is dumb
//os<<CRLF;

} else sames.push_back(order[i]);
plast=p;
}

if ( sames.size()!=0)
{
os<<sep<<sames.size()<<CRLF;
os<<MM_MARK<<" same as ";
for ( IdxTy ss=0; ss<sames.size(); ++ss) os<<sep<<sames[ss];
sames.clear();
os<<CRLF;

 } //sames
else os<<sep<<0<<CRLF;

//MM_MSG("")


// this is memory problem for large number of sites and sequences but usualyh not a big deal
// for hundred sites and thousand sequences.  ( about 1 Meg or so )

// these need to be numerical, and sorted.



}
static void invert_map(InverseLineIndexMapTy & d, const LineIndexMapTy & s)
{
typedef LineIndexMapTy::const_iterator II;
II i=s.begin();
while ( i!=s.end())
{
d[(*i).second].push_back((*i).first);
++i; 
}

}

static void invert_map(IndexMapTy & d, const StrIndexMapTy & s)
{
typedef StrIndexMapTy::const_iterator II;
II i=s.begin();
while ( i!=s.end())
{
d[(*i).second].push_back((*i).first);
++i; 
}

}



// turn line into string except leaving out field ex "rest of line"
static StrTy rol(const LineTy & line, const IdxTy ex,const StrTy & sep)
{
StrTy x;
const IdxTy len=line.size();
for  ( IdxTy i=0; i<len; ++i) if ( i!=ex) {if ( x.length()!=0) x+=sep; x+=line[i]; } 
return x; 
}
// make line from old line  leaving out field ex "rest of line"
static LineTy rolv(const LineTy & line, const IdxTy ex,const StrTy& sep=StrTy())
{
LineTy x;
const IdxTy len=line.size();
for  ( IdxTy i=0; i<len; ++i) if ( i!=ex) { x.push_back(line[i]); } 
return x; 
}


// given a list of key fields, make a key and value and then dump by value
static void index_by_key(StrIndexMapTy & d, const FileTy & s, const IdxTy key, const StrTy & sep)
{
const IdxTy sz=s.size();
const StrTy def=StrTy("");
for ( IdxTy i=0; i<sz; ++i)
{
const LineTy & line=s[i];
const IdxTy len=line.size();
const StrTy & k=(key<len)?line[key]:def;
d[k]+=rol(line,key,sep);

} //i 

}

// sep has not use here, easier to just default lol 
static void index_by_key(LineIndexMapTy & d, const FileTy & s, const IdxTy key, const StrTy & sep=StrTy(""))
{
const IdxTy sz=s.size();
const StrTy def=StrTy("");
for ( IdxTy i=0; i<sz; ++i)
{
const LineTy & line=s[i];
const IdxTy len=line.size();
const StrTy & k=(key<len)?line[key]:def;
LineTy lkey;
// in theory the key can be a series of fields... 
lkey.push_back(k);
d[lkey]=(rolv(line,key));

} //i 

}


static void index_by_key(GeneralLineIndexMapTy & d, const FileTy & s, const IdxTy key, const bool retain_key)
{
const IdxTy sz=s.size();
const StrTy def=StrTy("");
const IdxTy ekey=(retain_key)?(~0U):key;
for ( IdxTy i=0; i<sz; ++i)
{
const LineTy & line=s[i];
const IdxTy len=line.size();
const StrTy & k=(key<len)?line[key]:def;
LineTy lkey;
// in theory the key can be a series of fields... 
lkey.push_back(k);
d[lkey].push_back(rolv(line,ekey));

} //i 

}



// I suppose the keys should be a member of the file but this
// is not always the case 
// this could probably be done with a perl has more easily.. 
// append field fs from s to line in d with matching keys
// depending on dups flags, handle duplicates etc
static void append_matching_keys( FileTy & d, const IdxTy kd, const FileTy & s, const IdxTy ks, const IdxTy fs, const IdxTy dups)
{
enum { APPEND_MISSING=1, GLOB_EXTRA=2, PICK_ONE=4};
const bool append_missing=((dups&APPEND_MISSING)!=0);
const bool glob=((dups&GLOB_EXTRA)!=0);
const bool pick_one=((dups&PICK_ONE)!=0);
const StrTy missing="missing";
const StrTy sep=" "; // for globs 
//LineIndexMapTy di,si;
GeneralLineIndexMapTy si;
// this should be checked for unique keys... 
index_by_key(si,s,ks,true);
d.add_column(s.name(fs));
//index_by_key(di,d,kd);
for (auto i=d.begin(); i!=d.end(); ++i)
{
LineTy & ld=(*i);
// skip anything without a key field 
if (ld.size()<=kd) continue;
if (ld.is_comment())
{
// just add the names
if ( i==d.begin()) {ld.push_back(s.name(fs)); continue;}
}
const StrTy & keys=ld[kd];
LineTy key; key.push_back(keys);
const FileTy & sie=si[key];
const IdxTy sie_sz=sie.size();
MM_MSG(" with key "<<keys<<" have hits of "<<sie_sz)
if ( sie_sz==1)
{
// eh, this fs is now wrong.... key field was removed...
// see retain_key lol
ld.push_back(sie[0][fs]);

}else
{
if (sie_sz==0) { if (append_missing) ld.push_back(missing); }
else // more than one
{
StrTy v="";
if ( glob) {
// this needs to take all of the keys and make a new string 
v=sie[0][fs];
for ( IdxTy j=1; j<sie_sz; ++j) v+=sep+sie[j][fs];
ld.push_back(v);

} 
if ( pick_one) ld.push_back(sie[0][fs]);
} // multiple 
} // unique 



} // 

}



template < class Ty> static bool ok(const Ty & is) { return is.good()&&!is.eof(); } 
// load with extra junk at end, 
static IdxTy open_input( IsTy * & is, bool & dodel , const StrTy & nm) 
{
dodel=(nm!=StrTy("-"));
 is = dodel?(new std::ifstream(nm.c_str())):(&std::cin);
return 0;
}

static IdxTy open_output( OsTy * & is, bool & dodel , const StrTy & nm) 
{
dodel=(nm!=StrTy("-"));
 is = dodel?(new std::ofstream(nm.c_str())):(&std::cout);
return 0;
}



// this is really 2 steps, find stats, then make a dictionary from the stats
// using the rules applied to stats then finally substituting key and vale

static void make_dictionary_from_file(Dictionary & di, const FileTy & f, const DictionarySpec & ds)
{
// make the stats collection
ElementStats es;
es.load(f); // this should count all the crap
// probably return a dictionary created here from a factory?
// bit more likely to expand dictionaries than make new 
di.make(es,ds); // I guess the dictionary should be able to make itself`k
// the use the ds rules and elements to make a dictionary 

}

// this can be used for compression eventually but now just use
// for making concise names in fasta
// names to csv file, make dictionary, encode, extract back
// into fasta file 
// usually not encode whole set of columns, just some 
static void encode_file(FileTy & d, const FileTy & s, const Dictionary & di)
{



}




// dumb again, move switch out of line by line crap as always bulk
static void make_string_from_fields(StrTy & d, const LineTy &   s, const FormatRules & f)
{


switch ( f.branch())
{



// just put together with spaces 
default: {
const StrTy sep=f.sep(); // " ";
auto i=s.begin();
if ( i!=s.end()) { d+=(*i); ++i; }
for (; i!=s.end(); ++i) {d+=sep+(*i);}

} // default

} // switch 

}
// this should use quotes and escapes etc 
// this is reallyh dumb, the decision should be made on rules for the file
// not re checked after each line
static void parse_line_into_fields(LineTy & d, const ChTy *   s, const FormatRules & f)
{
// eventuyally this needs to handle commas, comments, and quotes etc
// and backslash lines etc
// this needs to concaternate whilein quote and accomodate escaped quotes
// or other conventions, perhaps also read "#" for comments or header
//if ( ok(ss)) x.push_back(y);
//SsTy ss((line));
switch ( f.branch())
{

default:
//SsTy ss(s); while ( ok(ss)) {StrTy y; ss>>y; if ( ok(ss)) d.push_back(y); }
// trailing spaces are making a mess... 
SsTy ss(s); while ( ok(ss)) {StrTy y=""; ss>>y; if ( true) d.push_back(y); }
}; // switch


}
// this is not named well, it refers to inner loop, see new methods to
// load and save LineTy as a text line with escapes and quotes
static bool load_line_from_string(LineTy & x, IsTy * is, ChTy * line, const IdxTy max)
{
is->getline(line,max);
if ( !ok(*is)) return false;
FormatRules fr;
parse_line_into_fields(x,line,fr);

return true;
}
// all fields after "fields" concatenated into a single field to
// make uniform 
static void squash_right(LineTy & x, const IdxTy fields, const StrTy & sep)
{
const IdxTy sz=x.size();
if ( sz<fields) return;
StrTy s=x[fields];
for ( IdxTy i=fields+1; i<sz; ++i) s=s+sep+x[i];
x[fields-1]=x[fields-1]+sep+s;
x.resize(fields);

}
// load a file with variable columns doing things to make it consistent
// either concateninating or dropping etc
static IdxTy load_ragged(FileTy & f, const StrTy & nm, const IdxTy fields, const IdxTy flags)
{
const IdxTy max=1<<16;
const StrTy sep=" ";
bool dodel=false; // =(nm!=StrTy("-"));
IsTy* is= 0; 
open_input(is,dodel,nm); // IsTy * is = dodel?(new std::ifstream(nm.c_str())):(&std::cin);
const IdxTy mask=(1<<16)-1;
const  IdxTy f1=flags&(mask);
const  IdxTy f2=(flags>>16)&(mask);
ChTy * line = new ChTy[max+2];
while (ok(*is))
{
LineTy x;
if ( !load_line_from_string(x,is,line,max)) break;
if ( x.is_comment())
{
if ( f.size()==0) f.set_columns(x); 
// to continue or not? May as well process for now... 
}

switch (f1)
{
case 0: { if ( x.size()>fields) squash_right(x,fields,sep); break ;}
case 1: { if ( x.size()>fields) x.resize(fields);  break ;}
case 2: { for ( IdxTy i=fields; i<x.size(); ++i)  
	MM_MSG(" too long "<<dump(x) )  break ;}
// too big squach too small filll
case 5: { if ( x.size()>fields) squash_right(x,fields,sep); 
case (8+5): { 
if ( x.size()!=fields) MM_MSG(" size is wrong fixing "<<dump(x));
if ( x.size()>fields) squash_right(x,fields,sep); 
 x.resize(fields);  break ;}

break ;}

default: if ( x.size()>fields) squash_right(x,fields,sep);

} // switch 

switch (f2)
{

case 1: { while ( x.size()<fields) x.push_back(LineTy::value_type());  break ;}


} // f2



f.push_back(x);
}
// mem leak on throw? 
delete [] line;
if ( dodel) delete is;
return 0;


}

static IdxTy write(FileTy & f, const StrTy & nm)
{

bool dodel=false; // delete when done, should be handled by object 
OsTy* os= 0; 

MM_MSG(" doing wirte now try to open  ")
MM_MSG(" doing wirte now try to open  "<<nm)

open_output(os,dodel,nm); 
FormatRules fr;
MM_MSG(" doing wirte now "<<nm)
if ( fr.write_header())
{
if ( fr.comment_header()) (*os)<<"#";
(*os)<<f.dump_columns(fr.sep());
(*os)<<CRLF;
}
IdxTy no=0;
for ( auto i=f.begin(); i!=f.end(); ++i)
{
// there is a function make_string_from_fields but it takes a format thing
const LineTy & line= (*i);
if ( line.is_comment()) 
	if ( no==0) {if (fr.write_header()) { ++no;  continue;}}


MM_MSG(" size ="<<line.size())
StrTy x="";
make_string_from_fields(x,line,fr);
(*os)<<x; (*os)<<CRLF;
++no;
} // i 


if ( dodel) delete os;
return 0; 


} // write

static IdxTy load(FileTy & f, const StrTy & nm)
{
const IdxTy max=1<<16;
//const bool dodel=(nm!=StrTy("-"));
// IsTy * is = dodel?(new std::ifstream(nm.c_str())):(&std::cin);
bool dodel=false; // =(nm!=StrTy("-"));
IsTy* is= 0; 
open_input(is,dodel,nm); // IsTy * is = dodel?(new std::ifstream(nm.c_str())):(&std::cin);


ChTy * line = new ChTy[max+2];
while (ok(*is))
{

// crap this was killing every ohter lines ... 
//is->getline(line,max);
//if ( !ok(*is)) break;
LineTy x;
if ( !load_line_from_string(x,is,line,max)) break;
if ( x.is_comment())
{
if ( f.size()==0) f.set_columns(x); 
// to continue or not? May as well process for now... 
}
//SsTy ss((line));
//LineTy x;
//while ( ok(ss)) { StrTy y; ss>>y; if ( ok(ss)) x.push_back(y); }
f.push_back(x);
}
// mem leak on throw? 
delete [] line;
if ( dodel) delete is;
return 0; 
}
 


}; // mjm_csv_ops

}; // mjm_csv ns

#endif


