#ifndef MJM_FORMAT_TEXT_H__
#define MJM_FORMAT_TEXT_H__

#include "mjm_globals.h"
// expanded functions 
// this makes include order scred up 
//#include "mjm_misc_parse.h"

#include "mjm_read_buffer.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*

g++ -DTEST_SHOOTING__ -Wall  -std=gnu++11 -gdwarf-3 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_shooting3.h

g++ -DTEST_FORMAT_TEXT__ -Wall  -gdwarf-3 -O0 -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -Wno-unused-function -I.. -x c++ mjm_format_text.h


*/

class mjm_format_text
{
typedef mjm_format_text Myt;
public:
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;


}; // Tr

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::StrTy StrTy;
typedef Tr::SsTy SsTy;
static StrTy sep() { return " "; } 
static StrTy npad() { return " "; } 
// right justify into field size s with precision p 
static StrTy format_float( const D & v, const IdxTy p, const IdxTy s)
{
SsTy ss;
ss.precision(p);
ss<<v;
StrTy sx;
ss>>sx;
IdxTy slen=sx.length();
SsTy ssf;
const int  pad=s-slen;
for (int i=0; i<pad; ++i) ssf<<npad();
ssf<<sx;
return ssf.str();
}
template <class Ti, class Te> static IdxTy max_sz_map(Ti  ii, const Te & ff)
{
IdxTy max=0;
while (ii!=ff) { if ((*ii).first.length()>max) max=(*ii).first.length(); ++ii; } 
return max;
}
static StrTy left(const StrTy & s, IdxTy sz)
{
StrTy x=s;
while (x.length()<sz) { x=x+" "; } 
return x; 
}
static StrTy align_on( const StrTy & s, const StrTy & pattern, const IdxTy col, const IdxTy size, const IdxTy flags)
{
StrTy x=s;
const char * xc=x.c_str();
const char * xp=pattern.c_str();
const IdxTy len=x.length();
IdxTy i=0;
if (pattern.length()>0)
{
for (i=0; i<len; ++i)if (xc[i]==xp[0]) break; 
while (i<(size-col)){ ++i;  x=" "+x; } 
while (x.length()<size) x=x+" ";

}

return x;
}


static StrTy prepend_lable(const StrTy & label, const StrTy & src,const IdxTy flags)
{
mjm_read_buffer<> rdb;


// this caps and makes string from base 
return rdb.string(); 
}


// eventually want to generate latex output inyo
// strings as it must be of finite length for tex lol
// need not stream 
// https://en.wikibooks.org/wiki/LaTeX/Tables



}; //mjm_format_text


#ifdef  TEST_FORMAT_TEXT__
// actual system
typedef double D;
typedef unsigned int IdxTy;


int main(int argc,char **args)
{

typedef mjm_format_text Myt;

typedef unsigned int IdxTy;
typedef double D;
Myt x;


return 0;

} //c main


#endif // TEST_SHOOTING__

#endif

