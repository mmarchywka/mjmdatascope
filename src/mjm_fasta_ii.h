#ifndef MJM_FASTA_II_H__
#define MJM_FASTA_II_H__

#include "mjm_text_data.h"
#include <stdlib.h>
// remove tcl as not needed I hope 2018-02-21
//#include <tcl.h>
//#include "mjm_tcl_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
#include <mjm_csv_ops.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

class fasta_ii_typedefs
{
public:
typedef unsigned int IdxTy;
typedef  int IntTy;
typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::stringstream SsTy;
};

namespace mjm_fasta 
{
typedef fasta_ii_typedefs Tr;

typedef Tr::IdxTy IdxTy;
typedef  Tr::IntTy IntTy;
typedef Tr::ChTy ChTy;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;


typedef StrTy NameTy;
//  the sequence should contain a signature...... 
// bad idea for memory management however.

typedef StrTy SeqTy;

typedef std::vector<IdxTy> IndexEntry;
typedef std::map<ChTy,IndexEntry > StrIndex;

// this crap should be in text data 
namespace text
{

// yet another white space finder
bool is_white(const ChTy c)
{ // do no make 0 white, needs to find eol lol 
if ( c==' ') return true;
if ( c=='\t') return true;
if ( c=='\r') return true;
if ( c=='\n') return true;

return false; 
}

// count the number of white chars from p.
IdxTy whites(const ChTy * p)
{
IdxTy i=0;
// needs to terminate on hitting zed 
while ( is_white(p[i])) ++i;  
return i;

}

// make a name conform to some criteria
StrTy nick_names(const StrTy & name, const IdxTy flags)
{
StrTy n;

return n;
}

class dummy
{
public:

static void index(StrIndex & d, const StrTy & s,const IdxTy start=0, const IdxTy len=0)
{
//IndexEntry;
const ChTy * p=s.c_str();
const IdxTy sz=((len==0))?s.length():len;
for ( IdxTy i=0; i<sz; ++i)
{ d[p[i]].push_back(i); } // i

}

// actually sequenial may be faster as thsi could thrash miserbly 
template <class Ty> static void xdensity(Ty & d, 
const IndexEntry & r, const IndexEntry& s )
{
const IdxTy szr=r.size();
const IdxTy szs=s.size();
for (IdxTy i=0; i<szr; ++i)
for (IdxTy j=0; j<szs; ++j)
d[r[i]][s[j]]=1;

}

static  void xdensity(OsTy & d, 
const IndexEntry & r, const IndexEntry& s )
{
const StrTy sep=" ";
const IdxTy szr=r.size();
const IdxTy szs=s.size();
for (IdxTy i=0; i<szs; ++i)
for (IdxTy j=0; j<szr; ++j)
d<<s[i]<<sep<<r[j]<<sep<<1<<CRLF;

}

static  void xdensity(OsTy & d, 
const IndexEntry & r, const IndexEntry& s,const IdxTy band )
{
const StrTy sep=" ";
const IdxTy szr=r.size();
const IdxTy szs=s.size();
//MM_MSG(" xdensity sizes "<<szr<<" "<<sze)
for (IdxTy i=0; i<szs; ++i)
for (IdxTy j=0; j<szr; ++j)
{
// ASSFUDD this is UN FUDDING SIGNED FUDDING AGAIN 
IntTy dd=IntTy(s[i]-r[j]);
//const IdxTy assfudd=dd;
// for same sequences this is a waste
if ( dd<0) dd=-dd;
const bool in_band=(dd<IntTy(band));
//MM_MSG(" FUDDER "<<dd<<" "<<assfudd<<" "<<in_band<<" "<<s[i]<<" "<<r[j])
if (in_band) d<<s[i]<<sep<<r[j]<<sep<<1<<CRLF;
}

}





template <class Ty> static void xdensity(Ty & d, const StrTy & r, const StrTy & s, 
const IdxTy ir, const IdxTy is, const IdxTy dr, const IdxTy ds)
{
StrIndex ri,si;
index(ri,r,ir,dr);
index(si,s,is,ds);
for (auto i=ri.begin(); i!=ri.end(); ++i)
 xdensity(d,(*i).second,si[(*i).first]);

} // xdensity

template <class Ty> static void xdensity(Ty & d, const StrTy & r, const StrTy & s, 
const IdxTy band)
{
StrIndex ri,si;
index(ri,r);
index(si,s);
for (auto i=ri.begin(); i!=ri.end(); ++i)
 xdensity(d,(*i).second,si[(*i).first],band);

} // xdensity






}; // dummy
}; // text  


class mjm_sequence_sign : public std::map<ChTy,IdxTy> 
{
typedef mjm_sequence_sign Myt;
typedef std::map<ChTy,IdxTy> Super;
typedef Super::const_iterator CI;
public:

StrTy dump(const StrTy & label="") const 
{
SsTy ss;
const StrTy sep=" ";
const StrTy sep2=" ";
const bool do_crlf=false; 
//const IdxTy sz=size();
IdxTy i=0;
for (auto ii=begin(); ii!=end(); ++ii)
{
if ( ii!=begin()) ss<<sep2;
ss<<label<<sep<<(*ii).first<<sep<<(*ii).second;
if ( do_crlf) ss<<CRLF;
++i;
//++ii;
}// i

return ss.str();
}

StrTy dump_more(const IdxTy & flags, const StrTy & label="") const 
{
SsTy ss;
const StrTy sep=" ";
const StrTy sep2=" ";
const bool do_crlf=false; 
//const IdxTy sz=size();
IdxTy i=0;
double pln=0;
double plnx=0;
const bool h0flag=true;
const bool hxflag=true;
double cnt=0;
double cntx=0;
const double scale=log(2.0);
for (auto ii=begin(); ii!=end(); ++ii)
{
if ( ii!=begin()) ss<<sep2;
const ChTy & cm=(*ii).first;
const IdxTy & n=(*ii).second; 
ss<<label<<sep<<cm<<sep<<n;
cnt+=n;
if ( cm!='.') if ( cm!='-') cntx+=n;
if ( do_crlf) ss<<CRLF;
++i;
//++ii;
}// i
for (auto ii=begin(); ii!=end(); ++ii)
{
const ChTy & cm=(*ii).first;
const IdxTy & n=(*ii).second; 
if ( cm!='.') if ( cm!='-')if (cntx!=0)  { const double px=1.0*n/cntx; plnx-=px*log(px)/scale; } 
if (cnt!=0)  { const double px=1.0*n/cnt; pln-=px*log(px)/scale; } 

}



if (h0flag) ss<<sep<< pln;
if (hxflag) ss<<sep<< plnx;

return ss.str();
}



// "hadamard"? LOL 
// should this be signed?
IntTy operator-(const Myt & that ) const
{
IntTy n=0;
for (auto ii=begin(); ii!=end(); ++ii)
{
const IdxTy x=(*ii).second;
IdxTy y=0;
ChTy key=(*ii).first;
auto ix=that.find(key);
if ( ix!=that.end()) y=(*ix).second;
n+=(x>y)?(x-y):(y-x);
}

return n;

} // minus

};
typedef mjm_sequence_sign Signature; 



typedef std::vector<Signature> Signatures;
// while maybe a good idea for OO, it is bad for performance
// as you want signatures in a sequential area for quick tests
// without memory thrashing and even get rid of large sequences
// from memory. 
class mjm_annointed_sequence
{

//void add(const NameTy & nm, const SeqTy & seq, const Signature &sig

}; // annointed

class fasta_file
{
typedef fasta_file Myt;
typedef std::vector<NameTy> Names;
typedef std::vector<SeqTy> Seqs;
typedef std::stringstream Ss; 
class _streaming_state_buffer
{
public:
_streaming_state_buffer() {is=0; Init(); }
_streaming_state_buffer(IsTy * s) {  is=s; Init(); }
~_streaming_state_buffer() { delete []  mbuf; }
bool ok() const { return valid_fasta; } 
void Init()
{
valid_fasta=false;  
have_line=false;  
maxlen=(1<<20)-3; 
mbuf=new char[maxlen+3]; 

} // Init
bool valid_fasta;
StrTy name,seq;
IdxTy maxlen;
IsTy * is;
char * mbuf;
bool have_line;
StrTy line;

}; // _streaming_state_buffer



// signature of a sequence, counts of AA or other eqauality
// result
// count amino acides, or something else for check sum 

public:
typedef _streaming_state_buffer streaming_state_buffer;
typedef streaming_state_buffer Ssb;
// the compiler is complaining here with non const name(fudd) 
void name(const StrTy & nm) { m_name=nm; }
const StrTy&  name()const { return m_name; }
const StrTy&  preamble()const { return m_preamble; }
// TODO FIXME  WTF the compiler keeps calling itslef with this ..... 
// oh it is the fudding const lol 
void write(const StrTy &nm ) const {write(nm.c_str()); } 
void write(const ChTy * f) const // it woulc not call the overload once it was const doh!!!
{
OsTy * os=new std::ofstream(f);
write(*os,60,0);
delete os;
}

void write ( OsTy & os, const IdxTy linesize, const IdxTy flags)  const
{ Write(os,linesize,flags); } 


void write_interleaved(const ChTy * f)
{
OsTy * os=new std::ofstream(f);
write_interleaved(*os,60,0);
delete os;
}
void write_interleaved ( OsTy & os, const IdxTy linesize, const IdxTy flags) 
{ WriteInterleaved(os,linesize,flags); } 
void next(streaming_state_buffer & ssb) { Next(ssb); } 
void load(const StrTy &  f) { load(f.c_str()); }
void load(const ChTy * f)
{ IsTy *is= new std::ifstream(f); load(is); delete is; }
void load_interleaved(const ChTy * f)
{ IsTy *is= new std::ifstream(f); load_interleaved(is); delete is; }


void load ( IsTy * is)  {Load(is); } 
void load_interleaved( IsTy * is)  {LoadInterleaved(is); } 
void check(OsTy *os) {Check(os); } 
void add(const NameTy & nm, const SeqTy & seq) { Add(nm, seq); }
void add(const NameTy & nm, const SeqTy & seq, const Signature &sig
 ) { Add(nm, seq,sig); }
void add(const Myt & s, const IdxTy i) // push file i from s
{ Add(s,i); } 
IdxTy size() const { return m_seqs.size();}
const SeqTy & seq(const IdxTy &i ) const { return m_seqs[i];}
const SeqTy & seq(const IdxTy &i,const StrTy & s )  {Set(i,s);  return m_seqs[i];}
// this is confusing with non-const see other name overload 
const StrTy & name(const IdxTy &i ) const { return m_names[i];}
// compiler warning wtf
const StrTy & name(const IdxTy &i ) { return m_names[i];}
const Signature & sig(const IdxTy &i ) const { return m_sigs[i];}
StrTy dump(const IdxTy flags=0) const
{ Ss ss; 
{ ss<<MMPR(size()); }
if (flags&(1))  return ss.str(); 
MM_SZ_LOOP(i,(*this),szf)
{
ss<<MMPR3(i,name(i),seq(i).length())<<";";
} // i 
return ss.str();
} // dump 
IdxTy remove_gt() 
{
IdxTy i=0;
MM_SZ_LOOP(j,m_names,szn)
{
const bool isgt=m_names[j].c_str()[0]=='>';
if (isgt){ m_names[j]=m_names[j].substr(1); ++i; } 
} // j 
return i;
} // remove_gt
private:
// stupid line readers.. 
enum {  MAXLEN=(1<<20)-2} ;
StrTy WriteName(const StrTy & n,const IdxTy flags=0)  const 
{ 
switch (flags)
{
// read the first "word" dropping the rest 
case 1: { StrTy x=""; SsTy s(n.c_str()); s>> x; return x; } 
// remove whte space ...
case 2: { StrTy x=""; SsTy s(n.c_str()); s>> x;
while (ok(&s)) { StrTy y=""; s>>y;if ( y.length()!=0) x=x+StrTy("_")+y;
 return x; } 
// remove whte space ...
}

} // switch 

return n; 
} 
void Write ( OsTy & os, const IdxTy linesize, const IdxTy flags)  const 
{
IdxTy sz= m_names.size();
for (IdxTy i=0; i<sz; ++i)
{
os<<">"<<WriteName(m_names[i])<<CRLF;
const ChTy *   s=m_seqs[i].c_str();
const IdxTy ss=m_seqs[i].length();
IdxTy pc=0;
IdxTy col=0;
while ( pc<ss)
{
// argh so slow....
os<<s[pc]; ++pc; ++col;
// could make two loops lol 
if ( linesize!=0) if ( col>=linesize) { os<<CRLF; col=0; }
}
if ( col!=0) os<<CRLF;
}
}


void WriteInterleaved ( OsTy & os, const IdxTy linesize, const IdxTy flags) 
{
const IdxTy sz= m_names.size();
IdxTy pc=0; 
const ChTy * p[sz];
IdxTy len[sz];
IdxTy minus[sz];
IdxTy maxlen=0;
for ( IdxTy i=0; i<sz; ++i) 
{p[i]=m_seqs[i].c_str(); 
len[i]=m_seqs[i].length();
if ( len[i]>maxlen) maxlen=len[i]; 
minus[i]=0; // dashed
}
while (pc<maxlen)
{
const IdxTy eol=pc+linesize;
for ( IdxTy i=0; i<sz; ++i) 
{
os<<WriteName(m_names[i],1); // this should only be done once as it takes time
os<<"\t";
const IdxTy eeol=(eol<len[i])?eol:len[i];
const ChTy * pi=p[i];
for ( IdxTy j=pc; j<eeol; ++j) { if ( pi[j]=='-') ++minus[i]; os<<pi[j]; } 
os<<" "<<(eeol-minus[i]);
os<<CRLF;

}

os<<CRLF;
pc=eol;
} //pc


} // writeinterleaved

// in process 


void LoadInterleaved ( IsTy * is) 
{
typedef std::map<StrTy,SeqTy> Seqs;
typedef std::map<StrTy,IdxTy> Counts;
Seqs smap;
Counts cmap;
{ // scoping  
char bf[MAXLEN+2];
StrTy s;
while ( ok(is))
{
bf[0]=0; (*is).getline(bf,MAXLEN); if ( !ok(is)) break;
s=StrTy(bf);
// format contains "junk" lines and is generally the clustalw format
// with id sequence count followed by matches and blank lines 
const ChTy * c=s.c_str();
IdxTy len =::strlen(c);
IdxTy endcount=0;
bool is_skip=(len==0);
SsTy ss(c);
StrTy nm="";
StrTy seq="";
StrTy count="";
StrTy crap="";
ss>>nm>>seq>>count;
if ( ok(&ss)) ss>>crap;
// there are nonsequence lines, they hopefully start with bad char 
// note this relies on naming conventions too as upper case names
// will fool it. 
if ( !is_skip)
{
// this is a fudd with upper case names, 
// eh, this is not supposed to be done on names rather word 2... 
const ChTy * cs=seq.c_str();
// note that the fasta could be "marked up" with other crap later
// right now this is just clusalw output
  is_skip=((( cs[0]<'A')||(cs[0]>'Z'))&&(cs[0]!='-')) ;
///if ( crap.length()
// count should make sense if it is there, otherwise name and sequence
// would go to bad location anyway but still a fudd
//if ( count!=0) 
if ( count.length()!=0) 
{
endcount=::atoi(count.c_str());
if ( endcount==0){
 is_skip=true;
MM_MSG(" skipping interleaved fasta cuz count is wrong "<<s)
}
else { if ( endcount<=cmap[nm])
{ MM_MSG(" count is out of order "<<endcount <<" vs "<< cmap[nm]<<" for "<<s)}
else cmap[nm]=endcount;
}

} // validate count
// try a space count, this should be more tolerant of
// consecutive spaces etc fudd shot fudd  
//IdxTy words=0; for ( IdxTy i=0; i<len; ++i) if((c[i]==' ')||(c[i]=='\t')) {++words} 
} // is_skip
// this should be id delim sequence dlim maybecount
if (!is_skip) { smap[nm]+=seq; }
} // ok 
} // scoping


// this is a big memory waste on large files but wtf.
// if you have aligned genomes, well lol
for(auto i=smap.begin(); i!=smap.end(); ++i)
{ Add((*i).first,(*i).second); }



}
static  bool GetLine ( StrTy & d, char * bf, IsTy * is, const IdxTy max)
{
bf[0]=0; (*is).getline(bf,MAXLEN); 
d=StrTy(bf);
return ( ok(is)) ;
}

void LoadGenbank ( IsTy * is) 
{
SeqTy seq;
char bf[MAXLEN+2];
while ( ok(is))
{
StrTy s;
// 2022-08-15 
//if ( !GetLine(s,bf,is,MAXLEN)) break; 
if ( !GetLine(s,bf,is,MAXLEN)) {} //  break; 
// eventually add a name, seq 

} // ok is

} //loadgenbank
void Next(streaming_state_buffer & ssb) { 
StrTy & seq=ssb.seq;
StrTy & name=ssb.name;
ssb.valid_fasta=false;
seq="";
name="";
bool have_name=false;
IsTy* is=ssb.is;
char * bf=ssb.mbuf;
MM_ERR(MMPR2(ok(is),ssb.have_line))
while ( ok(is)||ssb.have_line)
{
StrTy s;
// this may leave old result here is this is bad.. 
if (ssb.have_line) { s=ssb.line; ssb.have_line=false; } 
//else { if ( !GetLine(s,bf,is,ssb.maxlen)) break;  } 
// 2022-08-15
else { if ( !GetLine(s,bf,is,ssb.maxlen)) {}  } 
MM_ERR(MMPR(s))
// the zero char is always there, it may be 0 for eos... 
if ( s.c_str()[0]=='>')
{
if (have_name) 
{
ssb.have_line=true;
ssb.line=s;
return; 
} // have_name
// for now assume a fasta can be blank... 
else { name=StrTy(s.c_str()+1);  have_name=true; ssb.valid_fasta=true; } 
}
else if (have_name) seq+=s;
} // ok



}  // Next

/*
zymo preamble, 
::::::::::::::
./data/2022-02-18/zr5958_sv.seq.fna
::::::::::::::

*/

void Load( IsTy * is) 
{
//typedef std::vector<StrTy> LI;
//typedef std::vector<IdxTy> II;
//LI x,nm; II pos; 
SeqTy seq;
StrTy _preamble;
bool have_name=false;
// if we write these out without line endings, this will be a 
// mess reading back lol. 
char bf[MAXLEN+2];

//while ( (*is).good()&&!(*is).eof())
while ( ok(is))
{
StrTy s;
// this may leave old result here is this is bad.. 
// if the line is too long and can't tell from 
// crlf and could contain > then this is ambiguous
// but generally ok as long as name len less than MAXLEN
// 2022-08-15 this will delete the last line .... 
if ( !GetLine(s,bf,is,MAXLEN)) {} // break; 
//bf[0]=0; (*is).getline(bf,MAXLEN); if ( !ok(is)) break;
//s=StrTy(bf);
//MM_MSG("ASSFUDD "<<s)
// the zero char is always there, it may be 0 for eos... 
if ( s.c_str()[0]=='>')
{
//if (seq.size()!=0)
{
//  m_seqs.push_back(seq); 
if ( seq.length()){  Add(seq); have_name=false; } 
seq="";
} // push old one 
//SsTy ss(s.c_str()); //StrTy strim; //ss>>strim;
m_names.push_back(s);
// should have a blank sequence 
if (have_name) { Add(StrTy()); } 
have_name=true;
}
//else seq+=s;
// ignore leading crap 
else if (have_name) seq+=s;
else _preamble+=s;
} // ok

// this can allow name and seqs to be different size 
//if (seq.size()!=0)  m_seqs.push_back(seq);
if (have_name|| seq.length())
{
//MM_MSG(" FUDDER "<<seq)
 Add(seq);
}
m_preamble+=_preamble; // may load multiepls... 

} // load


void Add(const NameTy & nm, const SeqTy & seq)
{
m_names.push_back(nm);
Add(seq);
}
void Set(const IdxTy i, const StrTy & seq)
{
Signature s;
Sign(s,seq);
m_sigs[i]=s;
m_seqs[i]=seq;
}
void Add( const SeqTy & seq)
{
Signature s;
Sign(s,seq);
m_sigs.push_back(s);
 m_seqs.push_back(seq);

}
// the sequence triple needs to be an object 
void Add(const NameTy & nm, const SeqTy & seq, const Signature &sig)
{
m_names.push_back(nm);
m_sigs.push_back(sig);
 m_seqs.push_back(seq);
}
// there is the fasta editor 
void Add(const Myt & s, const IdxTy i) // push file i from s
{ Add(s.name(i),s.seq(i),s.sig(i)); }

bool Check(OsTy * os)
{
bool isok=true;
const IdxTy szn= m_names.size();
const IdxTy szseq= m_seqs.size();
const IdxTy szsig= m_sigs.size();
if (( szn!=szseq)||(szn!=szsig))
{
isok=false;
if ( os!=0) (*os)<<MM_MARK<<" sizes wrong naems="<<szn<<" seqs="<<szseq<<" sigs="<<szsig<<CRLF;


}


return isok;

}

 
static bool ok(const IsTy * is) { return  ( (*is).good()&&!(*is).eof()) ; }
static void Sign( Signature s, const SeqTy & seq)
{
const ChTy * p=seq.c_str();
const IdxTy sz=seq.length();
for ( IdxTy i=0; i!=sz; ++i) ++s[p[i]];

}

// MEMBERS

StrTy m_name;
StrTy m_preamble;
Names m_names;
Seqs m_seqs;
Signatures m_sigs;

}; // file 


typedef fasta_file  FileTy;

// iterator over a tgriangle of 2 files 
class triangle_iterator
{
typedef triangle_iterator Myt;
public:
triangle_iterator(const FileTy & si, const FileTy & sj )
:m_si(si),m_sj(sj),m_sz_i(si.size()),m_sz_j(sj.size()),
	m_i(1),m_j(0),m_good(m_sz_i+m_sz_j) {}

IdxTy i() const { return m_i; } 
IdxTy j() const { return m_j; } 
// this may require the compielr to optimize more ...`
const Signature & sigi()const { return m_si.sig(m_i);}
const Signature & sigj()const { return m_sj.sig(m_j);}

const SeqTy & seqi()const { return m_si.seq(m_j);}
const SeqTy & seqj()const { return m_sj.seq(m_j);}

const StrTy & namei() const { return m_si.name(m_i);}
const StrTy & namej() const { return m_sj.name(m_j);}
void clear() { Clear(); } 

operator bool() const { return m_good; } 
Myt & operator++()  { Inc(); return *this;} 



private:

void Clear()
{
m_i=1;
m_j=0;
m_good=((m_sz_i+m_sz_j)!=0);
}
void Inc() {++m_i; if ( m_i>=m_sz_i) { ++m_j;  m_i=m_j+1;
if((m_j>=m_sz_j)||(m_i>=m_sz_i)) m_good=false; }
}

const FileTy & m_si;
const FileTy & m_sj;
const IdxTy m_sz_i, m_sz_j;
// normally these go on stack,  could go on heap
IdxTy m_i,m_j;
bool m_good;
};


// this is a baseclass for an iterator over fasta
// it provides the traversal logic, the deriver provides
// interface. Note that some efficiency issues may come up.
class pair_iterator_members
{
typedef FileTy F;
protected:
pair_iterator_members(const F & si, const F& sj, const IdxTy i, const
IdxTy j, const bool good)
: m_si(si),m_sj(sj),m_sz_i(si.size()),m_sz_j(sj.size()),m_i(i),m_j(j),m_good(good) 
{}

bool is_null() const { return (m_i|m_j)==0; } 
//bool inc_i() { ++m_i; if (m_i<m_sz_i) return false; m_i=0; return true; } 
// let user re init. 
bool inc_i() { ++m_i; if (m_i<m_sz_i) return false;  return true; } 
bool in_range() const {return ((m_j<m_sz_j)&&(m_i<m_sz_i)) ;}
const FileTy & m_si;
const FileTy & m_sj;
const IdxTy m_sz_i, m_sz_j;
// normally these go on stack,  could go on heap
IdxTy m_i,m_j;
bool m_good;





};

class pair_iterator_logic_rectangle : protected pair_iterator_members
{
typedef pair_iterator_members Super;
public:
//pair_iterator_logic_rectangle(const FileTy & si, const FileTy & sj )
//:m_si(si),m_sj(sj),m_sz_i(si.size()),m_sz_j(sj.size()),
//	m_i(0),m_j(0),m_good(m_sz_i+m_sz_j) {}

pair_iterator_logic_rectangle(const FileTy & si, const FileTy & sj )
// lol, the compiler cuaght this. 
//:Super(si,sj,0,0,0!=(m_sz_i+m_sz_j)){}
:Super(si,sj,0,0,0!=(si.size()+sj.size())){}

protected:

void Clear() { Super::m_i=1; Super::m_j=0; Super::m_good=!Super::is_null(); }
void Inc() {if( Super::inc_i())  { ++Super::m_j;
//Super:: m_i=Super::m_j+1; // in this case, we waste an assignment to i. 
Super:: m_i=0; // in this case, we waste an assignment to i. 
 m_good=Super::in_range(); }}

} ;// rectangle


// if uyou want efficiency code the normal fing loop
template <class Ty> class fasta_pair_iterator : public Ty
{
// these names are wrong 
 typedef  Ty Parent;
 typedef  pair_iterator_members Super;
typedef fasta_pair_iterator<Ty> Myt;

public:
fasta_pair_iterator<Ty>(const FileTy & si, const FileTy & sj ) : Parent(si,sj){}

IdxTy i() const { return Super::m_i; } 
IdxTy j() const { return Super::m_j; } 
// this may require the compielr to optimize more ...`
const Signature & sigi()const { return Super::m_si.sig(i());}
const Signature & sigj()const { return Super::m_sj.sig(j());}

const SeqTy & seqi()const { return Super::m_si.seq(i());}
const SeqTy & seqj()const { return Super::m_sj.seq(j());}

const ChTy *  seqip()const { return Super::m_si.seq(i()).c_str();}
const ChTy *  seqjp()const { return Super::m_sj.seq(j()).c_str();}

const StrTy & namei() const { return Super::m_si.name(i());}
const StrTy & namej() const { return Super::m_sj.name(j());}
void clear() { Parent::Clear(); } 

operator bool() const { return Super::m_good; } 
Myt & operator++()  { Parent::Inc(); return *this;} 

private:


};// fasta pair iterator



typedef fasta_pair_iterator<pair_iterator_logic_rectangle> RectangleItor;




class fasta_ii 
{
//
//typedef LineTy::const_iterator LI;
//typedef FileTy::const_iterator FI;
//typedef IndexMapTy::const_iterator II;
//typedef LineIndexMapTy::const_iterator LII;
//typedef InverseLineIndexMapTy::const_iterator ILII;

typedef std::vector<IdxTy> IdxList;
typedef std::map<IdxTy,IdxList> OneToMany;

typedef StrIndex SeqIndex;

public:

// dump some subset of cross density or rc density between two seqs
// thatcould be identical 
static void xdensity()
{
//const SeqTy & ssi=s.seq(i);
//const StrTy & name=s.name(i);
//for (IdxTy j=0; j<szs; ++j)
}

static void rna_complement(SeqTy & d, const SeqTy & s)
{
const IdxTy sz=s.size();
const ChTy * p =s.c_str();
ChTy c[sz+1]; c[sz]=0;
for ( IdxTy i=0; i<sz; ++i)
{
ChTy b=p[i];
switch (b)
{
case 'A': { c[i]='T'; break; }
case 'T': { c[i]='A'; break; }
case 'C': { c[i]='G'; break; }
case 'G':  { c[i]='C'; break; }


};

} // i 

d= StrTy(c);

} //rna_complement



// test the rectnagle itor 
template <class Ty,class Tx > static void all_distances(Ty & d, const FileTy & r, const FileTy &s) 
{
// really we can just pass this, 
RectangleItor ii(r,s);
// really this should have an int as well as bool 
// although an index func is ok too 
while (ii)
{

const ChTy * pi=ii.seqip();
const ChTy * pj=ii.seqjp();
const IdxTy szi=ii.seqi().length();
const IdxTy szj=ii.seqj().length();
const IdxTy ds=szi-szj;
const IdxTy sz=(ds<0)?szj:szi;
IdxTy distance=(ds<0)?(-ds):ds;
for ( IdxTy k=0; k<sz; ++k) if ( pi[k]!=pj[k]) ++distance;
Tx tx;
tx.push_back(ii.i());
tx.push_back(ii.j());
tx.push_back(distance);

d.push_back(tx);
//d.push_back(Tx(ii.i(),ii.j(),distance));

++ii;

}
}

 static void stream_all_diffs(OsTy& os, const FileTy & r, const FileTy &s,const StrTy & label) 
{
// really we can just pass this, 
const StrTy sep=" ";
const StrTy blank="-";
RectangleItor ii(r,s);
// really this should have an int as well as bool 
// although an index func is ok too 
while (ii)
{

const ChTy * pi=ii.seqip();
const ChTy * pj=ii.seqjp();
const IdxTy szi=ii.seqi().length();
const IdxTy szj=ii.seqj().length();
const IdxTy ds=szi-szj;
// this should be taking the minimum
const IdxTy minsz=(ds<0)?szi:szj;
const IdxTy maxsz=(ds<0)?szj:szi;
IdxTy distance=(ds<0)?(-ds):ds;
// this counts the missing ones for distnace but ignores them for
// sites
IdxTy serial=0;
for ( IdxTy k=0; k<minsz; ++k) if ( pi[k]!=pj[k]) ++distance;
for ( IdxTy k=0; k<minsz; ++k) {
if ( pi[k]!=pj[k]){
os<<label<<sep<<(ii.i())<<sep<< (ii.j())<< sep<<(distance); 
os<<sep<<serial<<sep<<k<<sep<<pi[k]<<sep<<pj[k];
os<<CRLF;
++serial;
}}
for ( IdxTy k=minsz; k<maxsz; ++k) {

os<<label<<sep<<(ii.i())<<sep<< (ii.j())<< sep<<(distance); 
if (szi<szj)
os<<sep<<serial<<sep<<k<<sep<<blank<<sep<<pj[k];
else os<<sep<<serial<<sep<<k<<sep<<pi[k]<<sep<<blank;

os<<CRLF;
++serial;


}


++ii;

}
}



// find each element in s equal to 4
// put into some kind of list 
// this should just find "distances" based on signature diffs
// or actual diffs 
template <class Ty> static void list_duplicates(Ty & d, const FileTy & r, const FileTy &s) 
{

IdxTy szs=s.size();
IdxTy szr=r.size();
for (IdxTy i=0; i<szr; ++i)
{
// this triple could be an object but again we normally want sigs in one plae
// in memory 
const Signature & si=s.sig(i);
const SeqTy & ssi=s.seq(i);
const StrTy & name=s.name(i);
for (IdxTy j=0; j<szs; ++j)
{
const Signature & sj=r.sig(j);
const SeqTy & ssj=r.seq(j);
const StrTy & namej=r.name(j);
if ( si==sj) // this will be cheap and rarely true 
{
if ( ssi==ssj) 
{

} // really eq 

} // == 
} //j 
} //i

}


// put uniques into d from s, may want todelete "s" 
// on the fly to save space or make an index of dups. 
static void remove_duplicates(FileTy & d, const FileTy &s) 
{
IdxTy sz=s.size();
IdxTy szassert=sz; // s.m_seqs.size(); fudd 
if ( sz!=szassert) 
	MM_MSG(" fasta "<<s.name()<<" sigs="<<szassert<<" seqs="<<sz)
// doh, forgot to keep zero 
if ( sz>0)  d.add(s.name(0),s.seq(0),s.sig(0)); 
for (IdxTy i=1; i<sz; ++i)
{
const Signature & si=s.sig(i);
const SeqTy & ssi=s.seq(i);
const StrTy & name=s.name(i);
bool keep=true;
IdxTy j=0;
for (; j<i; ++j)
{
// the sequence and signature should be paired and define an eq operator
const Signature & sj=s.sig(j);
if ( si==sj)
{
	const SeqTy & ssj=s.seq(j);
// skip this one has already been used. 
if ( ssj==ssi) { keep=false; break; }
} 
} // j
if ( keep)
{ d.add(name,ssi,si); }
else // keep
{
MM_MSG(" dropping "<<i<<" "<<name<<" for equality  "<<j)
}



} // i



}

// there should be a state based range sleection thing etc 
// should be a class, vector of .pos and .state
template <class Ty> static StrTy make_name( const Ty & stt)
{
StrTy x=" range=";
SsTy s;
const IdxTy sz=stt.size();
IdxTy i=0;
while ( i<sz)
{
if ( stt[i]==0) { ++i; continue; }
IdxTy j=i+1;
while (( j<sz)&&(stt[j]==1)) ++j;
if (j>(i+1) )
//{ s<<i<<"-"<<(j-1)<<","; } else {s<<i; if (j<sz) s<<",";}
{ s<<i<<"-"<<(j-1); } else {s<<i; }
i=j;
if (j<sz) s<<",";
} // i 

return x+s.str(); 
}
// make a transition table from a string of items and 
// ranges
template <class Ty> static bool  make_stt(Ty & stt,const StrTy & s)
{
// for now take zero as 1
const IdxTy base_base=0; // we could number bases from 1... 
const ChTy * p=s.c_str();
ChTy t[s.length()+1];
IdxTy i=0;
IdxTy maxnum=0;
std::vector<IdxTy> ostack;
std::vector<ChTy> opstack;

while ( p[i]!=0)
{
IdxTy j=0;
while (( p[i]>='0')&&(p[i]<='9')) { t[j]=p[i]; ++j; ++i;}
t[j]=0;
if ( j==0) { opstack.push_back(p[i]); ++i; } 
else
{const IdxTy num= ::atoi(t)-base_base; 
if ( num>maxnum) maxnum=num; ostack.push_back(num); 
} 

} // i 
// this has too many complications, drunk and tired just enumerate

if ( stt.size()<=maxnum) stt.resize(maxnum+1); // still zero based 
const IdxTy stt_sz=stt.size();
const IdxTy DEF=0;
for ( IdxTy i=0; i<stt_sz; ++i) stt[i]=DEF;
IdxTy optr=0;
IdxTy opptr=0;
// go through each operand and mark spaces
const IdxTy osize=ostack.size();
const IdxTy opsize=opstack.size();

// just consume all operands, if operators remain then
// get end range filled in aka trailing hyphen
while ( optr<osize)
{
//const 
// may adjust this later
IdxTy n=ostack[optr]; // one number 
++optr;
// if there is no operator, just set this location

//if ( opptr>=opsize)  stt[n]=1;
if ( opptr>=opsize)  stt[n]=2;
else{
const ChTy opx=opstack[opptr];
++opptr;
 if ( opx==',') stt[n]=2;
else if ( opx=='-' ) 
{
if ( optr<osize)
{
const IdxTy m=ostack[optr]; // one number 
// optional spacing between segs 
if ( n<=m) {stt[n]=2; ++n; } 
for ( IdxTy i=n; i<=m; ++i) stt[i]=1;
++optr;
// here should be a comma next
if ( opptr<opsize) if ( opstack[opptr]!=',')
{
MM_MSG(" missing comma after range "<<s<<" at "<<m)
return false; 
}
++opptr;
} else
{ MM_MSG(" no operand , problem parsing "<<s<<" with "<<opx<<" and "<<n)
return false; 
} // osize

}
else
{ MM_MSG(" bad operator , parsing error for "<<s << " with opeartor "<< opx<<" on "<<n) 
return false; 
}
} //if else 


} // optr 


return true; 
}
// extract a list of ranges from fasta, stt is list of transition points
// names canbe appended with rnage info
template <class Ty> static void extract_ranges(FileTy &d, const FileTy &s,const Ty & stt, const bool update_names)
{

IdxTy sz=s.size();
IdxTy statesz=stt.size();
//IdxTy sc=0; // state counter
//IdxTy xt=~0U; // transition point
//if ( statesz>0) xt=stt[0].pos;
const StrTy sfx=make_name(stt);
for (IdxTy i=0; i<sz; ++i)
{
const SeqTy & si = s.seq(i);
// thismust be FIRST to escape space crap
 StrTy  nm = sfx+StrTy(":")+s.name(i);//+StrTy(":")+sfx;
const ChTy * p =si.c_str();
const IdxTy slen= si.length();
IdxTy pc=0;
// we added a choice of filler to separate segments
ChTy temp[2*slen+1];
for ( IdxTy j=0; j<slen; ++j)
{
//if ( j==xt)
//{ state=stt[sc].state; // if not, we do not hit this one again anyway... 
//++sc; if ( statesz>sc) xt=stt[sc].pos; }
const IdxTy state=(j<statesz)?stt[j]:0;
if ( state==1) {temp[pc]=p[j]; ++pc; } 
if ( state==2) {if ( pc!=0) {temp[pc]='.'; ++pc; } temp[pc]=p[j]; ++pc; } 

}
temp[pc]=0; 
d.add(nm,StrTy(temp));


}



}

// read the sequences and pick majority at each location
// create "d" with this sequence and maybe others. 

static void find_consensus(FileTy &d, const FileTy &s)
{
SsTy comments;
// this is the number of files, 
IdxTy sz=s.size();
SeqTy sd,notes;
typedef std::vector<Signature> Distro;
Distro distro; 
MM_MSG(" size is "<<sz)
for (IdxTy i=0; i<sz; ++i)
{
const SeqTy & si = s.seq(i);
const ChTy * p =si.c_str();
const IdxTy slen= si.length();
if ( slen>distro.size()){
MM_MSG(" resize tgo "<<slen)

 distro.resize(slen); 
}

for ( IdxTy j=0; j<slen; ++j) ++distro[j][p[j]];
} // i 
//MM_MSG("sht fycj")
sz=distro.size();
for (IdxTy i=0; i<sz; ++i)
{
//MM_MSG(" fucccckkkkk "<<i) 
const Signature & sig=distro[i];
if ( sig.size()==1)
{
const ChTy c=(*(sig.begin())).first;
sd+=c;
notes+='*';

} // size==1
else if ( sig.size()==0) {
//comments<<MM_MARK<<" zero map should not happen "<<s.name()<<" i="<<i<< CRLF;
MM_MSG(" zero map should not happen "<<s.name()<<" i="<<i) 
// << CRLF;
 sd+='?'; notes+='?';}
else
{
//MM_MSG("zzzzzzzzzzzzzzzzzzzz")
// find the max, record the number of entries.  Or
// fraction having tha majority?
ChTy c='x';
IdxTy max=0;
IdxTy total=0; // this should be the number of sequences
for ( auto j=sig.begin(); j!=sig.end(); ++j)
{
total+=(*j).second;
if( (*j).second>=max) { c=(*j).first; max=(*j).second; } 

} // j
sd+=c;
if ( false)
{
// number of acids, 
IdxTy cn=sig.size();
if ( cn>9) notes+='+';
else notes+=ChTy('0'+cn);
}
else
{
// fraction haing this cide, must be less than 1 
double pct=(10*max/total);
notes+=ChTy('0'+pct);

}

//comments<<MM_MARK<<" variation "<<s.name()<<" i="<<i<<" "<<sig.dump()<< CRLF;
comments<<MM_MARK<<" variation "<<s.name()<<" i="<<i<<" "<<sig.dump_more(3)<< CRLF;


} // size!=1

}


//MM_MSG("asssssssssssssssss")
d.add(StrTy("consensus ")+s.name(),sd);
MM_MSG(" not adding notes to consensus file for now ")
//d.add(StrTy("notes ")+s.name(),notes);

MM_MSG((comments.str()))

}


}; //  fasta_ii

}; // namespace  mjm_fasta


#endif

