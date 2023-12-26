#ifndef MJM_CPP_REGEX_H__
#define MJM_CPP_REGEX_H__
// TODO FIXME no idea why this was here may be needed for
//something else 2021-01-25 working on toobib
//#include "mjm_text_data.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include  "mjm_conflicts.h"
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>
#include <algorithm>
#include <regex>

#include "mjm_constants.h"


#define WORK_REGFUK (std::regex_constants::egrep | std::regex_constants::icase)
#define WORK_REGFUK_CASE (std::regex_constants::egrep)
// m_r(s.c_str(),REGFUK),m_nm(nm),m_string(s) {}

class regex_typedefs
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


namespace mjm_cpp_regex
{

typedef regex_typedefs Tr;
typedef Tr::IdxTy IdxTy;
typedef  Tr::IntTy IntTy;
typedef Tr::ChTy ChTy;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;


class stuff
{
public:
typedef std::regex Rx;

static Rx regex(const StrTy & x ) { return Rx(x,WORK_REGFUK); }
static Rx regex_case(const StrTy & x ) { return Rx(x,WORK_REGFUK_CASE); }
// this should allow a way to be compiled.  Make
// another class with non-statics for compiled or use regex struct
// is this greedy or not etc? Too few options. 
static void replace(StrTy & d, const StrTy & s, const StrTy& rgex, const StrTy & r, const IdxTy flags)
{
SsTy x;
    Rx w(rgex.c_str());
x<<std::regex_replace(s,w,r);
d=x.str();
}

 static bool has( const StrTy &  p , const Rx & regex, const IdxTy flags)
{ return has(p.c_str(),regex,flags); } 
 static bool has( const char * p , const Rx & regex, const IdxTy flags)
{
const char * sp=p;
//const char * e=sp+strlen(sp);
StrTy fudd=sp;
    auto words_begin = 
        std::sregex_iterator(fudd.begin(),fudd.end(), regex);
    auto words_end = std::sregex_iterator();
//MM_ERR(MMPR3(fudd,std::distance(words_begin,words_end), regex))
return (words_begin!=words_end);

}

 static StrTy  first( const char * p , const Rx & regex, const IdxTy flags)
{
const char * sp=p;
//const char * e=sp+strlen(sp);
StrTy fudd=sp;
    auto words_begin = 
        std::sregex_iterator(fudd.begin(),fudd.end(), regex);
    auto words_end = std::sregex_iterator();
//MM_ERR(MMPR3(fudd,std::distance(words_begin,words_end), regex))
const bool had= (words_begin!=words_end);
if (had)
{
//auto shot=*words_begin;
//MM_ERR(MMPR4(shot.position(),shot.length(),shot.str(),p))

}
return had?((*words_begin).str()):StrTy(""); 
}


template <class Ty> static void find(Ty & hits, const StrTy & s, const StrTy & regex, const IdxTy flags)
{

// http://en.cppreference.com/w/cpp/regex/regex_iterator
    Rx words_regex(regex.c_str());
find(hits,s,words_regex,flags); 
}
template <class Ty> static void find(Ty & hits, const StrTy & s, Rx & words_regex, const IdxTy flags)
{

    auto words_begin = 
        std::sregex_iterator(s.begin(), s.end(), words_regex);
    auto words_end = std::sregex_iterator();
 
// std::cout << "Found " << std::distance(words_begin, words_end) << " words:\n";
 
    for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
        std::smatch match = *i;                                                 
        std::string match_str = match.str(); 
const IdxTy pos=match.position();
        //std::cout << match_str << '\n';
       // std::cerr << match_str << '\n';
       MM_ERR(MMPR( match_str))

//typename typedef Ty::value_type Tv;
typedef typename Ty::value_type Tv;
hits.push_back(Tv(pos,match_str)); 

   }   

}





}; // stuff




}; // ns mjm_cpp_Regex

#endif



