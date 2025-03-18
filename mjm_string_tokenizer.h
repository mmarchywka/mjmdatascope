#ifndef MJM_STRING_TOKENIZER_H__
#define MJM_STRING_TOKENIZER_H__

#include "mjm_globals.h"
//#include "mjm_gen_viz.h"
//#include "mjm_data_model_error_log.h"
#include "mjm_generic_iterators.h"
//#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
//#include "mjm_calendar.h"
//#include "mjm_strings.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>

#define LUT_LOAD_FWD  0

// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
#define MM_DMELF(fm,code,x) f.DMel(code, f.m_ss<<MM_STR_LOC<<" "<<x); 

namespace string_tokenizer_traits
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
//typedef mjm_block_matrix<D> MyBlock;
//typedef  data_model_error_log Dmel;
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // biom_hdf5_traits

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
typedef string_tokenizer_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
//typedef Tr::MyBlock  MyBlock;
typedef StrTy _public_type;
typedef IdxTy _coded_type;

typedef std::map<_public_type, _coded_type> Fmap;
typedef std::vector<_public_type> Rmap;
public:
enum {BAD=~0};
typedef Myt mapped_type;
string_tokenizer():BADINT("invalid") {}
typedef _public_type public_type;
typedef _coded_type coded_type;
// decode the int n to the string if it is a  valid idx

const public_type&  operator()( const coded_type &  n ) const {if ( n<m_vec.size())  return m_vec[n]; return BADINT;  }
// note this makes a copy of the vector element into the other
// string s instead of just returning a ptr, 
const IdxTy  operator()( public_type & s, const coded_type &  n ) const {if ( n>=m_vec.size())  return BAD; s=m_vec[n]; return 0;  }
// note that ptrs are not stable if vector is modieied.. 
void operator()( const public_type * & p, const coded_type &  n ) const {if ( n<m_vec.size())  p=&m_vec[n];  }



const coded_type&  operator()( const char *  k  )  
{ return (*this)(public_type(k)); } 
const coded_type&  operator()( const public_type & k  )  {
auto x= m_map.find(k); if (x!=m_map.end()) return (*x).second;
const IdxTy n=m_vec.size();
m_map[k]=n;
m_vec.push_back(k);
 //return n; 
 return m_map[k];  // TODO just to return a ref lol?

}
// return true if this word has been coded or isin vocabulary
bool have(const public_type &k) const{return m_map.find(k)!=m_map.end();} 

typedef std::vector<coded_type> coded_vector;
typedef std::vector<public_type> public_vector;

coded_vector operator()( const public_vector &  keys)
{
coded_vector rv;
MM_LOOP(ii,keys) {rv.push_back((*this)((*ii))); }
return rv;
}

public_vector operator()( const coded_vector&  ikeys) const
{
public_vector rv;
MM_LOOP(ii,ikeys) {rv.push_back((*this)((*ii))); }
return rv;
}

// great so what about string valued map?????? 
template<class Tv> 
std::map<coded_type,Tv> operator()( const std::map<public_type,Tv>&  keys)
{
std::map<coded_type,Tv> rv;
MM_LOOP(ii,keys) {rv[(*this)((*ii).first)]=(*ii).second; }
return rv;
}




const coded_type&  operator[]( const public_type & k  )const   {
auto x= m_map.find(k); if (x!=m_map.end()) return (*x).second;

 return badidx();

}
Fmap::iterator find(const StrTy & k )  { return m_map.find(k); } 
Fmap::const_iterator find(const StrTy & k ) const { return m_map.find(k); } 


typedef vector_iterator<public_type,string_tokenizer_traits::Tr > Vi;
typedef const_vector_iterator<public_type,string_tokenizer_traits::Tr > CVi;

//Fmap::iterator 
Vi find(const coded_type & k )  
{ //MM_TRACE 
MM_ONCE(" this should not be called ",) 
return Vi(m_vec,k); // m_map.end();

 } 
//Fmap::const_iterator 
CVi find(const coded_type & k ) const 
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
// MEMBERS 
Fmap m_map;
Rmap m_vec;
StrTy BADINT;
}; // string_tokenizer


#undef MM_DMEL
#undef MM_DMELF

#endif

